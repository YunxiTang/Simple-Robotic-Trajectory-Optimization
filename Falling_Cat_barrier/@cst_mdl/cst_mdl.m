classdef cst_mdl < handle
    %CST_MODEL cost model
    
    properties
        Q,
        R,
        Qf,
        Rf,
    end
    
    methods
        function obj = cst_mdl(Q, R, Qf, Rf)
            %CST_MODEL Construct an instance of this class
            fprintf("[INFO]: Initializing Quadratic Cost Function. \n");
            obj.Q = Q;             % path state cost
            obj.R = R;             % path input cost (integral{1/2*(u'*R*u*dt)})
            obj.Qf = Qf;           % terminal cost (1/2*(x-x_star)'*Qf*(x-x_star))
            obj.Rf = Rf;
        end  
        function [] = update_Qf(obj,H)
            obj.Qf = H;
        end
    end
    
    methods (Static)
       %%% After code generation to class, Add the function declaration
       %%% here !!!
       l = l_cost(in1,in2,in3,in4);
       lf = lf_cost(in1);
       [l,lx,lu,lxx,lux,lxu,luu] = l_info(in1,in2,in3,in4);
       [lf,lfx,lfxx] = lf_info(in1);
       [Qx,Qu,Qxx,Quu,Qux,Qxu] = Q_info(rbt,cst,x,u,Vx,Vxx,params,iter);
       [Qx,Qu,Qxx,Quu,Qux,Qxu] = Qcs_info(rbt,cst,x,u,Vx,Vxx,params,cons)
    end
end

