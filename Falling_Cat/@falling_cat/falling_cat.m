classdef falling_cat < handle
    %Falling Cat Model
    
    properties
        Name = 'Falling Cat',
        g = 9.81,
        k = 6000,
        m1 = 2.5,
        m2 = 0.35,
        l1 = 0.15,
        l2 = 0.3,
        l2E = 0.3,
        l2R = 0.15,
        w = 0.2,
        d = 0.2,
        h = 0.1,
        I1 = [],
        I2 = zeros(3),
        
    end
    
    methods
        function obj = falling_cat()
            % Construct an instance of this class
            disp('[INFO]: Creating Falling Cat Model.');
            obj.I1 = 1/12*obj.m1*[ obj.d^2 + obj.h^2,  0,                   0;
                                   0,                  obj.w^2 + obj.h^2,   0;
                                   0,                  0,                   obj.w^2 + obj.d^2];
        end
        
        function qd = Dynamics(obj, t, q, u)
            %%% falling cat dynamics
            M_2 = obj.m2;
            I_xx = obj.I1(3,3);%I_xx in body frame is I_zz in box frame
            I_yy = obj.I1(2,2);
            I_zz = obj.I1(1,1);%I_zz in body frame is I_xx in box frame
            l_2 = obj.l2;
            %% ------------------------------------------
            theta_1 = q(1);
            theta_2 = q(2);
            theta_3 = q(3);
            theta_4 = q(4);
            theta_6 = q(5);
            d_7 = 0;
            theta_1d = q(6);
            theta_2d = q(7);
            theta_3d = q(8);
            theta_4d = q(9);
            theta_6d = q(10);
            d_7d = 0;
            %% ------------------------------------------pre-calculate sin cos functions
            S1 = sin(theta_1);
            S2 = sin(theta_2);
            S3 = sin(theta_3);
            S4 = sin(theta_4);
            S6 = sin(theta_6);
            C1 = cos(theta_1);
            C2 = cos(theta_2);
            C3 = cos(theta_3);
            C4 = cos(theta_4);
            C6 = cos(theta_6);
            %% equation of motion from recursive Newton Euler method
            G = zeros(5,1);
            R = zeros(5); 
            M(1,1) = S2*(I_zz*S2 + M_2*S4*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4) + M_2*C4*S6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)) + C2*(S3*(I_yy*C2*S3 - M_2*C6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)) + C3*(I_xx*C2*C3 - M_2*C4*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4) + M_2*S4*S6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)));
            M(1,2) = C2*(C3*(I_xx*S3 + M_2*C4^2*S3*(d_7 + l_2)^2 + M_2*S4*S6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)) - S3*(I_yy*C3 + M_2*C6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6))) + M_2*C4*C6*S2*(d_7 + l_2)^2*(C3*S6 - C6*S3*S4);
            M(1,3) = S2*(I_zz + M_2*S4^2*(d_7 + l_2)^2 + M_2*C4^2*S6^2*(d_7 + l_2)^2) - M_2*C2*C4*C6*(d_7 + l_2)^2*(S3*S6 + C3*C6*S4);
            M(1,4) = M_2*C6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6);
            M(1,5) = -M_2*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4);
            M(2,1) = S3*(I_xx*C2*C3 - M_2*C4*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4) + M_2*S4*S6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)) - C3*(I_yy*C2*S3 - M_2*C6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6));
            M(2,2) = S3*(I_xx*S3 + M_2*C4^2*S3*(d_7 + l_2)^2 + M_2*S4*S6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)) + C3*(I_yy*C3 + M_2*C6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6));
            M(2,3) = M_2*C4*C6*(d_7 + l_2)^2*(C3*S6 - C6*S3*S4);
            M(2,4) = M_2*C6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6);
            M(2,5) = M_2*C4*S3*(d_7 + l_2)^2;
            M(3,1) = I_zz*S2 + M_2*S4*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4) + M_2*C4*S6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6);
            M(3,2) = M_2*C4*C6*(d_7 + l_2)^2*(C3*S6 - C6*S3*S4);
            M(3,3) = I_zz + M_2*S4^2*(d_7 + l_2)^2 + M_2*C4^2*S6^2*(d_7 + l_2)^2;
            M(3,4) = M_2*C4*C6*S6*(d_7 + l_2)^2;
            M(3,5) = -M_2*S4*(d_7 + l_2)^2;
            M(4,1) = M_2*C6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6);
            M(4,2) = M_2*C6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6);
            M(4,3) = M_2*C4*C6*S6*(d_7 + l_2)^2;
            M(4,4) = M_2*C6^2*(d_7 + l_2)^2;
            M(4,5) = 0;
            M(5,1) = -M_2*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4);
            M(5,2) = M_2*C4*S3*(d_7 + l_2)^2;
            M(5,3) = -M_2*S4*(d_7 + l_2)^2;
            M(5,4) = 0;
            M(5,5) = M_2*(d_7 + l_2)^2;

            C(1,1) = theta_3d*((C2*(C3*(I_yy*C2*S3 - M_2*C6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)) + S3*(I_yy*C2*C3 + M_2*C2*C6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)) - C3*(I_xx*C2*S3 + M_2*C2*C4^2*S3*(d_7 + l_2)^2 + M_2*C2*S4*S6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)) - S3*(I_xx*C2*C3 - M_2*C4*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4) + M_2*S4*S6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6))))/2 - (M_2*C2*C4*C6*S2*(d_7 + l_2)^2*(C3*S6 - C6*S3*S4))/2) + d_7d*((S2*(M_2*S4*(S2*S4 - C2*C3*C4)*(2*d_7 + 2*l_2) + M_2*C4*S6*(2*d_7 + 2*l_2)*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)))/2 - (C2*(C3*(M_2*C4*(S2*S4 - C2*C3*C4)*(2*d_7 + 2*l_2) - M_2*S4*S6*(2*d_7 + 2*l_2)*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)) + M_2*C6*S3*(2*d_7 + 2*l_2)*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)))/2) + theta_2d*((S2*(I_zz*C2 + M_2*S4*(d_7 + l_2)^2*(C2*S4 + C3*C4*S2) + M_2*C4*S6*(d_7 + l_2)^2*(C6*S2*S3 + C2*C4*S6 - C3*S2*S4*S6)))/2 - (S2*(S3*(I_yy*C2*S3 - M_2*C6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)) + C3*(I_xx*C2*C3 - M_2*C4*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4) + M_2*S4*S6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6))))/2 - (C2*(S3*(I_yy*S2*S3 + M_2*C6*(d_7 + l_2)^2*(C6*S2*S3 + C2*C4*S6 - C3*S2*S4*S6)) + C3*(I_xx*C3*S2 + M_2*C4*(d_7 + l_2)^2*(C2*S4 + C3*C4*S2) - M_2*S4*S6*(d_7 + l_2)^2*(C6*S2*S3 + C2*C4*S6 - C3*S2*S4*S6))))/2 + (C2*(I_zz*S2 + M_2*S4*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4) + M_2*C4*S6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)))/2) + M_2*theta_6d*(d_7 + l_2)^2*(C4^2*C6*S2^2*S6 - C2^2*C6*S3^2*S6 - C2*C4*C6^2*S2*S3 + C2*C4*S2*S3*S6^2 - C2^2*C3*C6^2*S3*S4 + C2^2*C3*S3*S4*S6^2 + C2^2*C3^2*C6*S4^2*S6 + 2*C2*C3*C4*C6*S2*S4*S6) + M_2*theta_4d*C6*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4)*(C2*S3*S6 + C4*C6*S2 + C2*C3*C6*S4);
            C(1,2) = theta_3d*((C2*(I_zz + M_2*S4^2*(d_7 + l_2)^2 + M_2*C4^2*S6^2*(d_7 + l_2)^2))/2 - (C2*(S3*(I_xx*S3 + M_2*C4^2*S3*(d_7 + l_2)^2 + M_2*S4*S6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)) + C3*(I_yy*C3 + M_2*C6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)) - S3*(I_yy*S3 + M_2*C6*(d_7 + l_2)^2*(C6*S3 - C3*S4*S6)) - C3*(I_xx*C3 + M_2*C3*C4^2*(d_7 + l_2)^2 - M_2*S4*S6*(d_7 + l_2)^2*(C6*S3 - C3*S4*S6))))/2) + d_7d*((C2*(C3*(M_2*C4^2*S3*(2*d_7 + 2*l_2) + M_2*S4*S6*(C3*C6 + S3*S4*S6)*(2*d_7 + 2*l_2)) - M_2*C6*S3*(C3*C6 + S3*S4*S6)*(2*d_7 + 2*l_2)))/2 + (M_2*C4*C6*S2*(C3*S6 - C6*S3*S4)*(2*d_7 + 2*l_2))/2) + theta_1d*((S2*(I_zz*C2 + M_2*S4*(d_7 + l_2)^2*(C2*S4 + C3*C4*S2) + M_2*C4*S6*(d_7 + l_2)^2*(C6*S2*S3 + C2*C4*S6 - C3*S2*S4*S6)))/2 - (S2*(S3*(I_yy*C2*S3 - M_2*C6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)) + C3*(I_xx*C2*C3 - M_2*C4*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4) + M_2*S4*S6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6))))/2 - (C2*(S3*(I_yy*S2*S3 + M_2*C6*(d_7 + l_2)^2*(C6*S2*S3 + C2*C4*S6 - C3*S2*S4*S6)) + C3*(I_xx*C3*S2 + M_2*C4*(d_7 + l_2)^2*(C2*S4 + C3*C4*S2) - M_2*S4*S6*(d_7 + l_2)^2*(C6*S2*S3 + C2*C4*S6 - C3*S2*S4*S6))))/2 + (C2*(I_zz*S2 + M_2*S4*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4) + M_2*C4*S6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)))/2) - theta_2d*(S2*(C3*(I_xx*S3 + M_2*C4^2*S3*(d_7 + l_2)^2 + M_2*S4*S6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)) - S3*(I_yy*C3 + M_2*C6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6))) - M_2*C2*C4*C6*(d_7 + l_2)^2*(C3*S6 - C6*S3*S4)) - (M_2*theta_6d*(d_7 + l_2)^2*(2*C2*C3^2*S4 + 2*C2*C6^2*S4 + 2*C3*C4*S2 - 2*C3*C4*C6^2*S2 - 4*C2*C3^2*C6^2*S4 - 4*C2*C3*C6*S3*S6 - 2*C4*C6*S2*S3*S4*S6 + 2*C2*C3*C4^2*C6*S3*S6))/2 - M_2*theta_4d*C6*(d_7 + l_2)^2*(C4^2*C6*S2*S3 - C2*C3^2*C4*S6 - C6*S2*S3 + C3*S2*S4*S6 + C2*C3*C4*C6*S3*S4);
            C(1,3) = theta_1d*((C2*(C3*(I_yy*C2*S3 - M_2*C6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)) + S3*(I_yy*C2*C3 + M_2*C2*C6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)) - C3*(I_xx*C2*S3 + M_2*C2*C4^2*S3*(d_7 + l_2)^2 + M_2*C2*S4*S6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)) - S3*(I_xx*C2*C3 - M_2*C4*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4) + M_2*S4*S6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6))))/2 - (M_2*C2*C4*C6*S2*(d_7 + l_2)^2*(C3*S6 - C6*S3*S4))/2) + d_7d*((S2*(M_2*S4^2*(2*d_7 + 2*l_2) + M_2*C4^2*S6^2*(2*d_7 + 2*l_2)))/2 - (M_2*C2*C4*C6*(S3*S6 + C3*C6*S4)*(2*d_7 + 2*l_2))/2) + theta_2d*((C2*(I_zz + M_2*S4^2*(d_7 + l_2)^2 + M_2*C4^2*S6^2*(d_7 + l_2)^2))/2 - (C2*(S3*(I_xx*S3 + M_2*C4^2*S3*(d_7 + l_2)^2 + M_2*S4*S6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)) + C3*(I_yy*C3 + M_2*C6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)) - S3*(I_yy*S3 + M_2*C6*(d_7 + l_2)^2*(C6*S3 - C3*S4*S6)) - C3*(I_xx*C3 + M_2*C3*C4^2*(d_7 + l_2)^2 - M_2*S4*S6*(d_7 + l_2)^2*(C6*S3 - C3*S4*S6))))/2) + M_2*theta_4d*C4*C6^2*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4) + M_2*theta_6d*C4*C6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6) - M_2*theta_3d*C2*C4*C6*(d_7 + l_2)^2*(C3*S6 - C6*S3*S4);
            C(1,4) = (M_2*theta_6d*(d_7 + l_2)^2*(2*C4*C6^2*S2 - 2*C4*S2 - 2*C2*C3*S4 + 2*C2*C3*C6^2*S4 + 2*C2*C6*S3*S6))/2 + (M_2*d_7d*C6*(2*d_7 + 2*l_2)*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6))/2 - M_2*theta_2d*C6*(d_7 + l_2)^2*(C4^2*C6*S2*S3 - C2*C3^2*C4*S6 - C6*S2*S3 + C3*S2*S4*S6 + C2*C3*C4*C6*S3*S4) + M_2*theta_1d*C6*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4)*(C2*S3*S6 + C4*C6*S2 + C2*C3*C6*S4) + M_2*theta_3d*C4*C6^2*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4) - M_2*theta_4d*C6*S6*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4);
            C(1,5) = (M_2*theta_4d*(d_7 + l_2)^2*(2*C4*C6^2*S2 - 2*C4*S2 - 2*C2*C3*S4 + 2*C2*C3*C6^2*S4 + 2*C2*C6*S3*S6))/2 - (M_2*d_7d*(S2*S4 - C2*C3*C4)*(2*d_7 + 2*l_2))/2 - (M_2*theta_2d*(d_7 + l_2)^2*(2*C2*C3^2*S4 + 2*C2*C6^2*S4 + 2*C3*C4*S2 - 2*C3*C4*C6^2*S2 - 4*C2*C3^2*C6^2*S4 - 4*C2*C3*C6*S3*S6 - 2*C4*C6*S2*S3*S4*S6 + 2*C2*C3*C4^2*C6*S3*S6))/2 + M_2*theta_1d*(d_7 + l_2)^2*(C4^2*C6*S2^2*S6 - C2^2*C6*S3^2*S6 - C2*C4*C6^2*S2*S3 + C2*C4*S2*S3*S6^2 - C2^2*C3*C6^2*S3*S4 + C2^2*C3*S3*S4*S6^2 + C2^2*C3^2*C6*S4^2*S6 + 2*C2*C3*C4*C6*S2*S4*S6) + M_2*theta_3d*C4*C6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6);
            C(2,1) = (M_2*theta_6d*(d_7 + l_2)^2*(2*C2*S4 - 2*C2*C3^2*S4 - 2*C2*C6^2*S4 + 2*C3*C4*C6^2*S2 + 4*C2*C3^2*C6^2*S4 + 4*C2*C3*C6*S3*S6 + 2*C4*C6*S2*S3*S4*S6 - 2*C2*C3*C4^2*C6*S3*S6))/2 - d_7d*((S3*(M_2*C4*(S2*S4 - C2*C3*C4)*(2*d_7 + 2*l_2) - M_2*S4*S6*(2*d_7 + 2*l_2)*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)))/2 - (M_2*C3*C6*(2*d_7 + 2*l_2)*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6))/2) - theta_1d*((S2*(I_zz*C2 + M_2*S4*(d_7 + l_2)^2*(C2*S4 + C3*C4*S2) + M_2*C4*S6*(d_7 + l_2)^2*(C6*S2*S3 + C2*C4*S6 - C3*S2*S4*S6)))/2 - (S2*(S3*(I_yy*C2*S3 - M_2*C6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)) + C3*(I_xx*C2*C3 - M_2*C4*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4) + M_2*S4*S6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6))))/2 - (C2*(S3*(I_yy*S2*S3 + M_2*C6*(d_7 + l_2)^2*(C6*S2*S3 + C2*C4*S6 - C3*S2*S4*S6)) + C3*(I_xx*C3*S2 + M_2*C4*(d_7 + l_2)^2*(C2*S4 + C3*C4*S2) - M_2*S4*S6*(d_7 + l_2)^2*(C6*S2*S3 + C2*C4*S6 - C3*S2*S4*S6))))/2 + (C2*(I_zz*S2 + M_2*S4*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4) + M_2*C4*S6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)))/2) - theta_3d*((C3*(I_yy*C2*C3 + M_2*C2*C6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)))/2 - (S3*(I_yy*C2*S3 - M_2*C6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)))/2 + (I_zz*C2)/2 + (S3*(I_xx*C2*S3 + M_2*C2*C4^2*S3*(d_7 + l_2)^2 + M_2*C2*S4*S6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)))/2 - (C3*(I_xx*C2*C3 - M_2*C4*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4) + M_2*S4*S6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)))/2 + (M_2*S4*(d_7 + l_2)^2*(C2*S4 + C3*C4*S2))/2 + (M_2*C4*S6*(d_7 + l_2)^2*(C6*S2*S3 + C2*C4*S6 - C3*S2*S4*S6))/2) - M_2*theta_4d*C4*C6*(d_7 + l_2)^2*(C2*S6 - C2*C3^2*S6 + C4*C6*S2*S3 + C2*C3*C6*S3*S4);
            C(2,2) = d_7d*((S3*(M_2*C4^2*S3*(2*d_7 + 2*l_2) + M_2*S4*S6*(C3*C6 + S3*S4*S6)*(2*d_7 + 2*l_2)))/2 + (M_2*C3*C6*(C3*C6 + S3*S4*S6)*(2*d_7 + 2*l_2))/2) + theta_3d*((S3*(I_xx*C3 + M_2*C3*C4^2*(d_7 + l_2)^2 - M_2*S4*S6*(d_7 + l_2)^2*(C6*S3 - C3*S4*S6)))/2 + (C3*(I_xx*S3 + M_2*C4^2*S3*(d_7 + l_2)^2 + M_2*S4*S6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)))/2 - (S3*(I_yy*C3 + M_2*C6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)))/2 - (C3*(I_yy*S3 + M_2*C6*(d_7 + l_2)^2*(C6*S3 - C3*S4*S6)))/2) - M_2*theta_6d*(d_7 + l_2)^2*(C3^2*C6*S6 - C3*C6^2*S3*S4 + C3*S3*S4*S6^2 - C6*S3^2*S4^2*S6) + M_2*theta_4d*C4*C6*S3*(d_7 + l_2)^2*(C3*S6 - C6*S3*S4);
            C(2,3) = theta_2d*((S3*(I_xx*C3 + M_2*C3*C4^2*(d_7 + l_2)^2 - M_2*S4*S6*(d_7 + l_2)^2*(C6*S3 - C3*S4*S6)))/2 + (C3*(I_xx*S3 + M_2*C4^2*S3*(d_7 + l_2)^2 + M_2*S4*S6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)))/2 - (S3*(I_yy*C3 + M_2*C6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)))/2 - (C3*(I_yy*S3 + M_2*C6*(d_7 + l_2)^2*(C6*S3 - C3*S4*S6)))/2) - theta_1d*((C2*(I_zz + M_2*S4^2*(d_7 + l_2)^2 + M_2*C4^2*S6^2*(d_7 + l_2)^2))/2 - (S3*(I_yy*C2*S3 - M_2*C6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)))/2 + (C3*(I_yy*C2*C3 + M_2*C2*C6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)))/2 + (S3*(I_xx*C2*S3 + M_2*C2*C4^2*S3*(d_7 + l_2)^2 + M_2*C2*S4*S6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)))/2 - (C3*(I_xx*C2*C3 - M_2*C4*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4) + M_2*S4*S6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)))/2 + (M_2*C4*C6*S2*(d_7 + l_2)^2*(S3*S6 + C3*C6*S4))/2) + (M_2*d_7d*C4*C6*(C3*S6 - C6*S3*S4)*(2*d_7 + 2*l_2))/2 + (M_2*theta_4d*C6^2*S3*(2*S4^2 - 2)*(d_7 + l_2)^2)/2 - M_2*theta_3d*C4*C6*(d_7 + l_2)^2*(S3*S6 + C3*C6*S4) + M_2*theta_6d*C4*C6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6);
            C(2,4) = (M_2*d_7d*C6*(C3*C6 + S3*S4*S6)*(2*d_7 + 2*l_2))/2 - M_2*C4*C6*(d_7 + l_2)^2*(theta_1d*C2*S6 + theta_2d*C6*S4 - theta_4d*S3*S6 + theta_3d*C4*C6*S3 - theta_2d*C3*S3*S6 - theta_1d*C2*C3^2*S6 - theta_2d*C3^2*C6*S4 + theta_1d*C4*C6*S2*S3 + theta_1d*C2*C3*C6*S3*S4) - (M_2*theta_6d*(d_7 + l_2)^2*(2*S3*S4 - 2*C6^2*S3*S4 + 2*C3*C6*S6))/2;
            C(2,5) = (M_2*(d_7 + l_2)^2*(theta_2d*sin(2*theta_6) + 2*theta_1d*C2*S4 - 2*theta_4d*S3*S4 - 2*theta_4d*C3*C6*S6 - 2*theta_2d*C3*S3*S4 + 2*theta_3d*C3*C4*C6^2 - 2*theta_1d*C2*C3^2*S4 - 2*theta_1d*C2*C6^2*S4 - 4*theta_2d*C3^2*C6*S6 - 2*theta_2d*C4^2*C6*S6 + 2*theta_4d*C6^2*S3*S4 + 4*theta_1d*C2*C3^2*C6^2*S4 + 2*theta_2d*C3^2*C4^2*C6*S6 + 2*theta_1d*C3*C4*C6^2*S2 + 4*theta_2d*C3*C6^2*S3*S4 + 4*theta_1d*C2*C3*C6*S3*S6 + 2*theta_3d*C4*C6*S3*S4*S6 + 2*theta_1d*C4*C6*S2*S3*S4*S6 - 2*theta_1d*C2*C3*C4^2*C6*S3*S6))/2 + (M_2*d_7d*C4*S3*(2*d_7 + 2*l_2))/2;
            C(3,1) = d_7d*((M_2*S4*(S2*S4 - C2*C3*C4)*(2*d_7 + 2*l_2))/2 + (M_2*C4*S6*(2*d_7 + 2*l_2)*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6))/2) - theta_1d*((C2*(C3*(I_yy*C2*S3 - M_2*C6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)) + S3*(I_yy*C2*C3 + M_2*C2*C6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)) - C3*(I_xx*C2*S3 + M_2*C2*C4^2*S3*(d_7 + l_2)^2 + M_2*C2*S4*S6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)) - S3*(I_xx*C2*C3 - M_2*C4*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4) + M_2*S4*S6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6))))/2 - (M_2*C2*C4*C6*S2*(d_7 + l_2)^2*(C3*S6 - C6*S3*S4))/2) + theta_2d*((C3*(I_yy*C2*C3 + M_2*C2*C6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)))/2 - (S3*(I_yy*C2*S3 - M_2*C6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)))/2 + (I_zz*C2)/2 + (S3*(I_xx*C2*S3 + M_2*C2*C4^2*S3*(d_7 + l_2)^2 + M_2*C2*S4*S6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)))/2 - (C3*(I_xx*C2*C3 - M_2*C4*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4) + M_2*S4*S6*(d_7 + l_2)^2*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6)))/2 + (M_2*S4*(d_7 + l_2)^2*(C2*S4 + C3*C4*S2))/2 + (M_2*C4*S6*(d_7 + l_2)^2*(C6*S2*S3 + C2*C4*S6 - C3*S2*S4*S6))/2) + M_2*theta_4d*C6*(d_7 + l_2)^2*(C2*C3*C6 - C2*C3*C4^2*C6 + C4*C6*S2*S4 + C2*S3*S4*S6) + M_2*theta_6d*C4*(d_7 + l_2)^2*(C2*S3 - C2*C6^2*S3 + C4*C6*S2*S6 + C2*C3*C6*S4*S6);
            C(3,2) = theta_1d*((C2*(S3*(I_xx*S3 + M_2*C4^2*S3*(d_7 + l_2)^2 + M_2*S4*S6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)) + C3*(I_yy*C3 + M_2*C6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)) - S3*(I_yy*S3 + M_2*C6*(d_7 + l_2)^2*(C6*S3 - C3*S4*S6)) - C3*(I_xx*C3 + M_2*C3*C4^2*(d_7 + l_2)^2 - M_2*S4*S6*(d_7 + l_2)^2*(C6*S3 - C3*S4*S6))))/2 + (I_zz*C2)/2 + (M_2*S4*(d_7 + l_2)^2*(C2*S4 + C3*C4*S2))/2 + (M_2*C4*S6*(d_7 + l_2)^2*(C6*S2*S3 + C2*C4*S6 - C3*S2*S4*S6))/2 + (M_2*C4*C6*S2*(d_7 + l_2)^2*(S3*S6 + C3*C6*S4))/2) - theta_2d*((S3*(I_xx*C3 + M_2*C3*C4^2*(d_7 + l_2)^2 - M_2*S4*S6*(d_7 + l_2)^2*(C6*S3 - C3*S4*S6)))/2 + (C3*(I_xx*S3 + M_2*C4^2*S3*(d_7 + l_2)^2 + M_2*S4*S6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)))/2 - (S3*(I_yy*C3 + M_2*C6*(d_7 + l_2)^2*(C3*C6 + S3*S4*S6)))/2 - (C3*(I_yy*S3 + M_2*C6*(d_7 + l_2)^2*(C6*S3 - C3*S4*S6)))/2) - M_2*theta_4d*C6*(d_7 + l_2)^2*(C3*S4*S6 - C6*S3 + C4^2*C6*S3) + M_2*theta_6d*C4*(d_7 + l_2)^2*(C3*C6^2 - C3 + C6*S3*S4*S6) + (M_2*d_7d*C4*C6*(C3*S6 - C6*S3*S4)*(2*d_7 + 2*l_2))/2;
            C(3,3) = d_7d*((M_2*S4^2*(2*d_7 + 2*l_2))/2 + (M_2*C4^2*S6^2*(2*d_7 + 2*l_2))/2) + M_2*C4*C6*(d_7 + l_2)^2*(theta_4d*C6*S4 + theta_6d*C4*S6);
            C(3,4) = M_2*C6*(d_7 + l_2)^2*(theta_2d*C6*S3 - theta_4d*S4*S6 + theta_1d*C2*C3*C6 + theta_3d*C4*C6*S4 - theta_2d*C3*S4*S6 - theta_2d*C4^2*C6*S3 + theta_1d*C4*C6*S2*S4 + theta_1d*C2*S3*S4*S6 - theta_1d*C2*C3*C4^2*C6) + M_2*theta_6d*C4*(d_7 + l_2)^2*(C6^2 - 1) + (M_2*d_7d*C4*C6*S6*(2*d_7 + 2*l_2))/2;
            C(3,5) = (M_2*C4*(d_7 + l_2)^2*(2*theta_4d*C6^2 - 2*theta_4d - 2*theta_2d*C3 + 2*theta_1d*C2*S3 + 2*theta_2d*C3*C6^2 + 2*theta_3d*C4*C6*S6 - 2*theta_1d*C2*C6^2*S3 + 2*theta_1d*C4*C6*S2*S6 + 2*theta_2d*C6*S3*S4*S6 + 2*theta_1d*C2*C3*C6*S4*S6))/2 - (M_2*d_7d*S4*(2*d_7 + 2*l_2))/2;
            C(4,1) = (M_2*d_7d*C6*(2*d_7 + 2*l_2)*(C4*S2*S6 - C2*C6*S3 + C2*C3*S4*S6))/2 + M_2*theta_6d*C6*(d_7 + l_2)^2*(C2*S3*S6 + C4*C6*S2 + C2*C3*C6*S4) - M_2*theta_3d*C6*(d_7 + l_2)^2*(C2*C3*C6 - C2*C3*C4^2*C6 + C4*C6*S2*S4 + C2*S3*S4*S6) - M_2*theta_1d*C6*(d_7 + l_2)^2*(S2*S4 - C2*C3*C4)*(C2*S3*S6 + C4*C6*S2 + C2*C3*C6*S4) + M_2*theta_2d*C4*C6*(d_7 + l_2)^2*(C2*S6 - C2*C3^2*S6 + C4*C6*S2*S3 + C2*C3*C6*S3*S4);
            C(4,2) = M_2*C4*C6*(d_7 + l_2)^2*(theta_1d*C2*S6 + theta_2d*C6*S4 - theta_2d*C3*S3*S6 - theta_1d*C2*C3^2*S6 - theta_2d*C3^2*C6*S4 + theta_1d*C4*C6*S2*S3 + theta_1d*C2*C3*C6*S3*S4) - M_2*theta_6d*C6*(d_7 + l_2)^2*(C3*S6 - C6*S3*S4) + (M_2*d_7d*C6*(C3*C6 + S3*S4*S6)*(2*d_7 + 2*l_2))/2 + M_2*theta_3d*C6*(d_7 + l_2)^2*(C3*S4*S6 - C6*S3 + C4^2*C6*S3);
            C(4,3) = (M_2*d_7d*C4*C6*S6*(2*d_7 + 2*l_2))/2 - M_2*C6*(d_7 + l_2)^2*(theta_2d*C6*S3 - theta_6d*C4*C6 + theta_1d*C2*C3*C6 + theta_3d*C4*C6*S4 - theta_2d*C3*S4*S6 - theta_2d*C4^2*C6*S3 + theta_1d*C4*C6*S2*S4 + theta_1d*C2*S3*S4*S6 - theta_1d*C2*C3*C4^2*C6);
            C(4,4) = (M_2*d_7d*C6^2*(2*d_7 + 2*l_2))/2 - (M_2*theta_6d*sin(2*theta_6)*(d_7 + l_2)^2)/2;
            C(4,5) = M_2*C6*(d_7 + l_2)^2*(theta_3d*C4*C6 - theta_4d*S6 - theta_2d*C3*S6 + theta_1d*C4*C6*S2 + theta_1d*C2*S3*S6 + theta_2d*C6*S3*S4 + theta_1d*C2*C3*C6*S4);
            C(5,1) = - (M_2*d_7d*(S2*S4 - C2*C3*C4)*(2*d_7 + 2*l_2))/2 - M_2*theta_1d*(d_7 + l_2)^2*(C4^2*C6*S2^2*S6 - C2^2*C6*S3^2*S6 - C2*C4*C6^2*S2*S3 + C2*C4*S2*S3*S6^2 - C2^2*C3*C6^2*S3*S4 + C2^2*C3*S3*S4*S6^2 + C2^2*C3^2*C6*S4^2*S6 + 2*C2*C3*C4*C6*S2*S4*S6) - (M_2*theta_2d*(d_7 + l_2)^2*(2*C2*S4 - 2*C2*C3^2*S4 - 2*C2*C6^2*S4 + 2*C3*C4*C6^2*S2 + 4*C2*C3^2*C6^2*S4 + 4*C2*C3*C6*S3*S6 + 2*C4*C6*S2*S3*S4*S6 - 2*C2*C3*C4^2*C6*S3*S6))/2 - M_2*theta_4d*C6*(d_7 + l_2)^2*(C2*S3*S6 + C4*C6*S2 + C2*C3*C6*S4) - M_2*theta_3d*C4*(d_7 + l_2)^2*(C2*S3 - C2*C6^2*S3 + C4*C6*S2*S6 + C2*C3*C6*S4*S6);
            C(5,2) = (M_2*d_7d*C4*S3*(2*d_7 + 2*l_2))/2 - (M_2*(d_7 + l_2)^2*(theta_2d*sin(2*theta_6) - 2*theta_3d*C3*C4 + 2*theta_1d*C2*S4 - 2*theta_4d*C3*C6*S6 - 2*theta_2d*C3*S3*S4 + 2*theta_3d*C3*C4*C6^2 - 2*theta_1d*C2*C3^2*S4 - 2*theta_1d*C2*C6^2*S4 - 4*theta_2d*C3^2*C6*S6 - 2*theta_2d*C4^2*C6*S6 + 2*theta_4d*C6^2*S3*S4 + 4*theta_1d*C2*C3^2*C6^2*S4 + 2*theta_2d*C3^2*C4^2*C6*S6 + 2*theta_1d*C3*C4*C6^2*S2 + 4*theta_2d*C3*C6^2*S3*S4 + 4*theta_1d*C2*C3*C6*S3*S6 + 2*theta_3d*C4*C6*S3*S4*S6 + 2*theta_1d*C4*C6*S2*S3*S4*S6 - 2*theta_1d*C2*C3*C4^2*C6*S3*S6))/2;
            C(5,3) = - M_2*C4*(d_7 + l_2)^2*(theta_4d*C6^2 - theta_2d*C3 + theta_1d*C2*S3 + theta_2d*C3*C6^2 + theta_3d*C4*C6*S6 - theta_1d*C2*C6^2*S3 + theta_1d*C4*C6*S2*S6 + theta_2d*C6*S3*S4*S6 + theta_1d*C2*C3*C6*S4*S6) - (M_2*d_7d*S4*(2*d_7 + 2*l_2))/2;
            C(5,4) = -M_2*C6*(d_7 + l_2)^2*(theta_3d*C4*C6 - theta_4d*S6 - theta_2d*C3*S6 + theta_1d*C4*C6*S2 + theta_1d*C2*S3*S6 + theta_2d*C6*S3*S4 + theta_1d*C2*C3*C6*S4);
            C(5,5) = (M_2*d_7d*(2*d_7 + 2*l_2))/2;
            %% ------------------------------------------
            dy = q(6:10,:);
            tau = [0;0;0;u(1);u(2)];
            qd = [dy; M\(tau-C*dy-R*dy-G)];
        end
        
        function q_next = rk(obj, q, u, dt)
           k1 = obj.Dynamics(0, q,             u);
           k2 = obj.Dynamics(0, q + 0.5*dt*k1, u);
           k3 = obj.Dynamics(0, q + 0.5*dt*k2, u);
           k4 = obj.Dynamics(0, q +     dt*k3, u);
           q_next = q + dt/6*(k1+2*k2+2*k3+k4);
        end
        function [fx,fu] = getLinSys(obj,qbar,ubar,dt)
            dh = 1e-4;
            Nx = numel(qbar);
            Nu = numel(ubar);
            Jx = zeros(Nx,Nx);
            Ju = zeros(Nx,Nu);
            Hx = eye(Nx) * dh;
            Hu = eye(Nu) * dh;
            for i=1:Nx
                Jx(:,i) = (obj.rk(qbar + Hx(:,i), ubar, dt) - obj.rk(qbar - Hx(:,i), ubar, dt)) / (2*dh);
            end

            for i=1:Nu
                Ju(:,i) = (obj.rk(qbar, ubar+ Hu(:,i), dt) - obj.rk(qbar , ubar - Hu(:,i), dt)) / (2*dh);
            end
            fx = Jx;
            fu = Ju;
        end
        
    end
    methods (Static)
        [fx, fu] = Finite_Diff(func, qbar, ubar, dt, interval)
    end
end

