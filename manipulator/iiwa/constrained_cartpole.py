from pydrake.common import FindResourceOrThrow, temp_directory
from pydrake.geometry import (
    MeshcatVisualizer,
    MeshcatVisualizerParams,
    Role,
    StartMeshcat,
)

from pydrake.all import (DrakeVisualizer, 
                        ConnectContactResultsToDrakeVisualizer,
                        MultibodyPlant)

from pydrake.multibody.parsing import Parser
from pydrake.multibody.plant import AddMultibodyPlantSceneGraph

from pydrake.systems.primitives import ConstantVectorSource

from pydrake.systems.analysis import Simulator
from pydrake.systems.framework import DiagramBuilder
from pydrake.math import RigidTransform, RollPitchYaw
import numpy as np

# basic settings
dt = 0.01


####################################
# system setup
####################################

def add_robot_model_into_plant(plant):
    sdf = FindResourceOrThrow("drake/examples/multibody/cart_pole/cart_pole.sdf")
    parser = Parser(plant)
    robot_model = parser.AddModelFromFile(sdf)
    plant.Finalize()
    return plant, robot_model

####################################
# Create system diagram
####################################
builder = DiagramBuilder()
plant, scene_graph = AddMultibodyPlantSceneGraph(builder, time_step=dt)
plant, _ = add_robot_model_into_plant(plant)
assert plant.geometry_source_is_registered(), 'Geometry object is not registered'

controller = builder.AddSystem(ConstantVectorSource(np.zeros(1)))
builder.Connect(
    controller.get_output_port(),
    plant.get_actuation_input_port()
)

DrakeVisualizer().AddToBuilder(builder, scene_graph)
ConnectContactResultsToDrakeVisualizer(builder, plant, scene_graph)

diagram = builder.Build()
diagram.set_name("cartpole")
diagram_context = diagram.CreateDefaultContext()
plant_context = plant.GetMyMutableContextFromRoot(diagram_context)
print(plant.CalcTotalMass(plant_context))
print(plant_context)

# Create system model for controller
plant_ = MultibodyPlant(time_step=dt)
plant_, _ = add_robot_model_into_plant(plant_)
input_port_index = plant_.get_actuation_input_port().get_index()
# print(type(plant_))


method = 'ilqr'

if method == "ilqr":
    # Set up the optimizer
    from ilqr import cons_iLQR
    num_steps = 200
    ilqr = cons_iLQR(
        plant_, 
        num_steps, 
        input_port_index=input_port_index,
        beta=0.9
    )

    # Define initial and target states
    # Initial state
    x0 = np.array([0,np.pi-0.1,0,0])

    # Target state
    x_nom = np.array([0,np.pi,0,0])

    ilqr.SetInitialState(x0)
    ilqr.SetTargetState(x_nom)

    # Define cost function
    # Quadratic cost int_{0^T} (x'Qx + u'Ru) + x_T*Qf*x_T
    Q = np.diag([10,10,0.1,0.1])
    R = 0.001*np.eye(1)
    Qf = np.diag([100,100,10,10])
    ilqr.SetRunningCost(dt*Q, dt*R)
    ilqr.SetTerminalCost(Qf)

    # Set initial guess
    u_guess = np.zeros((1,num_steps-1))
    ilqr.SetInitialGuess(u_guess)

    states, inputs, solve_time, optimal_cost = ilqr.Solve()
    print(f"Solved in {solve_time} seconds using iLQR")
    print(f"Optimal cost: {optimal_cost}")
    timesteps = np.arange(0, num_steps) * dt
