from pydrake.all import *
import numpy as np
import os
import time

class cons_iLQR:
    """constrained iLQR solver for Trajectory Optimization within the Drake framewprk"""
    def __init__(self, system, input_port_index:int, N:int, beta:float) -> None:
        """
        Args:
            system:             Drake Multibody Plant describing discrete dynamics
                                with x_{t+1} = f(x_t,u_t).
            N:                  Number of timesteps (Time Horizon)
            input_port_index:   InputPortIndex for the control input u_t. 
                                Default is to use the first port. 
            
        """
        # system check
        assert system.IsDifferenceEquationSystem()[0],  "Only discrete-time system is supported" 
        
        self.system = system
        self.context = system.CreateDefaultContext()
        self.input_port = system.get_input_port(input_port_index)

        # Autodiff version of the system 
        # for computing dynamics gradients
        self.system_ad = system.ToAutoDiffXd()
        self.context_ad = self.system_ad.CreateDefaultContext()
        self.input_port_ad = self.system_ad.get_input_port(input_port_index)

        # Horizon Length
        self.N = N

        # Define state and input sizes
        self.n = self.context.get_discrete_state_vector().size()
        self.m = self.input_port.size()

        # Initial and target states
        self.x0 = np.zeros(self.n)
        self.x_xom = np.zeros(self.n)

        # reference trajectory
        self.xref = None
        self.uref = None

        # Quadratic cost terms
        self.Q = np.eye(self.n)
        self.R = np.eye(self.m)
        self.Qf = np.eye(self.n)

        # <X, U> trajectory
        self.x_bar = np.zeros((self.n, self.N))
        self.u_bar = np.zeros((self.m, self.N-1))

        # Arrays to store dynamics gradients
        self.fx = np.zeros((self.n,self.n,self.N-1))
        self.fu = np.zeros((self.n,self.m,self.N-1))

        # Local feedback controller
        # u = u_bar - eps*k_t - K_t@(x-x_bar)
        self.k = np.zeros((self.m, self.N-1))
        self.K = np.zeros((self.m, self.n, self.N-1))

        # Coefficents Qu'*Quu^{-1}*Qu for computing the expected 
        # reduction in cost dV = sum_t eps*(1-eps/2)*Qu'*Quu^{-1}*Qu
        self.dV_coeff = np.zeros(self.N-1)

    def SetInitialState(self, x0):
        """set the initial state"""
        self.x0 = x0

    def SetTargetState(self, x_nom):
        """set the target state"""
        self.x_nom = x_nom

    def SetReferenceTraj(self, xref, uref):
        """set the reference trajectory"""
        self.xref = xref
        self.uref = uref

    def AddCost(self, cost_func, k:int):
        """add cost at the k-th time"""
        pass

    def SetRunningCost(self, Q, R):
        """
            set the quadratic running cost: 
            (x-x_nom)' @ Q @ (x-x_nom) + u' @ R @ u
        """
        assert Q.shape == (self.n,self.n)
        assert R.shape == (self.m,self.m)

        self.Q = Q
        self.R = R

    def SetTerminalCost(self, Qf):
        """
            set the terminal cost: 
            (x-x_nom)' @ Qf @ (x-x_nom)
        """
        assert Qf.shape == (self.n, self.n)
        self.Qf = Qf
    
    def SetInitialGuess(self, u_guess):
        """set the initial guess of control"""
        assert u_guess.shape == (self.m, self.N-1)
        self.u_bar = u_guess

    def SetControlLimits(self, u_min, u_max):
        pass

    def _running_cost_partials(self, x, u, step):
        """
        get the partial derivatives of the running cost
            l = (x-x_ref)' @ Q @ (x-x_ref) + u' @ R @ u
        for the given state and input values.
        Args:
            x:  state
            u:  control
            step: time step
        Returns:
            lx:     1st order partial w.r.t. x
            lu:     1st order partial w.r.t. u
            lxx:    2nd order partial w.r.t. x
            luu:    2nd order partial w.r.t. u
            lux:    2nd order partial w.r.t. u and x
        """
        lx = 2 * self.Q @ (x - self.x_nom)
        lu = 2 * self.R @ u
        lxx = 2 * self.Q
        luu = 2 * self.R
        lux = np.zeros((self.m,self.n))

        return (lx, lu, lxx, luu, lux)

    def _terminal_cost_partials(self, x):
        """
        Return the partial derivatives of the (quadratic) terminal cost
            lf = x' @ Qf @x
        Args:
            x: numpy array representing state
        Returns:
            lf_x:   gradient of terminal cost
            lf_xx:  hessian of terminal cost
        """
        lf_x = 2*self.Qf @ x - 2*self.x_nom.T @ self.Qf
        lf_xx = 2*self.Qf

        return (lf_x, lf_xx)

    def _calc_dynamics(self, x, u):
        """
        numerical-value type: compute the next state 
            x_next = f(x,u)
        Args:   
            x:  An (n,) numpy array representing the state
            u:  An (m,) numpy array representing the control input
        Returns:
            x_next: An (n,) numpy array representing the next state
        """
        # Set input and state variables in our stored model accordingly
        self.context.SetDiscreteState(x)
        self.input_port.FixValue(self.context, u)

        # Compute the forward dynamics x_next = f(x,u)
        state = self._context.get_discrete_state()
        self._system.CalcDiscreteVariableUpdates(self.context, state)
        x_next = state.get_vector().value().flatten()
        return x_next

    def _calc_dynamics_partials(self, x, u):
        """
        compute the first-order partial derivitives of the dynamics
        
        Args:   
            x: (n,) numpy array representing the state
            u: (m,) numpy array representing the control input
        Returns:
            fx: (n,n) numpy array representing the partial derivative 
                of f with respect to x.
            fu: (n,m) numpy array representing the partial derivative 
                of f with respect to u.
        """
        # Create autodiff versions of x and u
        xu = np.hstack([x,u])
        xu_ad = InitializeAutoDiff(xu)
        x_ad = xu_ad[:self.n]
        u_ad = xu_ad[self.n:]

        # Set input and state variables in stored model accordingly
        self.context_ad.SetDiscreteState(x_ad)
        self.input_port_ad.FixValue(self.context_ad, u_ad)

        # Compute the forward dynamics x_next = f(x,u)
        state = self.context_ad.get_discrete_state()
        self.system_ad.CalcDiscreteVariableUpdates(self.context_ad, state)
        x_next = state.get_vector().CopyToVector()
       
        # Compute partial derivatives
        G = ExtractGradient(x_next)
        fx = G[:,:self.n]
        fu = G[:,self.n:]

        return (fx, fu)