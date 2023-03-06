import numpy as np
from pydrake.all import *
from free_ilqr import free_iLQR
import time
from pydrake.geometry import (
    MeshcatVisualizer,
    MeshcatVisualizerParams,
    Role,
    StartMeshcat,
)

####################################
# Parameters
####################################

T = 2.0      # total simulation time (S)
dt = 1e-3     # simulation timestep

# Solver method
# must be "ilqr" or "sqp"
method = "ilqr"

# Initial state
q0 = np.array([0.0, 0.1, 0, -1.2, 0, 0.8, 0])
x0 = np.hstack((q0, 0*q0))

# Target state
x_nom = np.copy( x0 )
x_nom[0] = 1.0
x_nom[1] = 1.0

# Quadratic cost int_{0^T} (x'Qx + u'Ru) + x_T*Qf*x_T
Q = np.diag([50.0] * 7 + [2.] * 7) * dt
R = np.diag([0.0] * 7) * dt
Qf = np.diag([5000] * 7 + [10.] * 7) * dt

####################################
# Tools for system setup
####################################

def create_system_model(plant):
    # sdf = FindResourceOrThrow("drake/examples/acrobot/Acrobot.urdf")
    sdf = FindResourceOrThrow("drake/manipulation/models/iiwa_description/sdf/iiwa14_no_collision.sdf")
    robot = Parser(plant=plant).AddModelFromFile(sdf)
    plant.WeldFrames(plant.world_frame(), plant.GetFrameByName("iiwa_link_0", robot))
    plant.Finalize()
    return plant, robot

####################################
# Create system diagram
####################################
builder = DiagramBuilder()
plant, scene_graph = AddMultibodyPlantSceneGraph(builder, dt)
plant, _ = create_system_model(plant)
assert plant.geometry_source_is_registered()

controller = builder.AddSystem(ConstantVectorSource(np.zeros(7)))
builder.Connect(
    controller.get_output_port(),
    plant.get_actuation_input_port()
    )

# Start the visualizer.
meshcat = StartMeshcat()
meshcat.SetProperty("/Background", "visible", False)
visualizer = MeshcatVisualizer.AddToBuilder(builder, scene_graph, meshcat)

diagram = builder.Build()
diagram_context = diagram.CreateDefaultContext()
plant_context = diagram.GetMutableSubsystemContext(plant, diagram_context)

#####################################
# Solve Trajectory Optimization
#####################################
    
# Create system model for controller
plant_ = MultibodyPlant(dt)
plant_, _ = create_system_model(plant_)
input_port_index = plant_.get_actuation_input_port().get_index()

#-----------------------------------------
# DDP method
#-----------------------------------------

if method == "ilqr":
    # Set up the optimizer
    num_steps = int(T/dt)
    ilqr = free_iLQR(
        plant_, 
        num_steps, 
        input_port_index=input_port_index,
        beta=0.9
    )

    # Define initial and target states
    ilqr.SetInitialState(x0)
    ilqr.SetTargetState(x_nom)

    # Define cost function
    ilqr.SetRunningCost(Q, R)
    ilqr.SetTerminalCost(Qf)

    # Set initial guess
    plant.SetPositionsAndVelocities(plant_context, x0)
    tau_g = -plant.CalcGravityGeneralizedForces(plant_context)
    S = plant.MakeActuationMatrix().T
    u_gravity_comp = 1.0 *  S @ np.repeat(tau_g[np.newaxis].T, num_steps-1, axis=1)
    u_guess = u_gravity_comp
    ilqr.SetInitialGuess(u_guess)

    states, inputs, solve_time, optimal_cost = ilqr.Solve()
    print(f"Solved in {solve_time} seconds using iLQR")
    print(f"Optimal cost: {optimal_cost}")
    timesteps = np.arange(0.0, T, dt)

#-----------------------------------------
# Direct Transcription method
#-----------------------------------------

elif method == "sqp":
    context_ = plant_.CreateDefaultContext()
    # Set up the solver object
    trajopt = DirectTranscription(
            plant_, 
            context_, 
            input_port_index=input_port_index,
            num_time_samples=int(T/(dt))
        )
    
    # Add constraints
    x = trajopt.state()
    u = trajopt.input()

    x_init = trajopt.initial_state()
    x_final = trajopt.final_state()
    
    trajopt.prog().AddConstraint(eq( x_init, x0 ))
    trajopt.prog().AddConstraint(eq( x_final, x_nom ))
    x_err = x - x_nom
    trajopt.AddRunningCost(u.T @ R @ u)
    # trajopt.AddFinalCost(x_err.T @ Qf @ x_err)
    
    # Solve the optimization problem
    st = time.time()
    res = Solve(trajopt.prog())
    solve_time = time.time() - st
    assert res.is_success(), "trajectory optimizer failed"
    solver_name = res.get_solver_id().name()
    optimal_cost = res.get_optimal_cost()
    print(f"Solved in {solve_time} seconds using SQP via {solver_name}")
    print(f"Optimal cost: {optimal_cost}")
    
    # Extract the solution
    timesteps = trajopt.GetSampleTimes(res)
    states = trajopt.GetStateSamples(res)
    inputs = trajopt.GetInputSamples(res)

else:
    raise ValueError(f"Unknown method {method}")

#####################################
# Playback
#####################################

# Just keep playing back the trajectory
visualizer_context = visualizer.GetMyContextFromRoot(diagram_context)
# visualizer.StartRecording(False)
import matplotlib.pyplot as plt
for i in range(2):
    plt.plot(timesteps, states[i,:])
plt.show()

while True:
    for i in range(len(timesteps)):
        t = timesteps[i]
        x = states[:,i]
        
        diagram_context.SetTime(t)
        plant.SetPositionsAndVelocities(plant_context, x)
        visualizer.ForcedPublish(visualizer_context)

        time.sleep(5e-3)
        
    time.sleep(3.0)



