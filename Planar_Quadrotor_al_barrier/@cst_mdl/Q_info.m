function [Qx,Qu,Qxx,Quu,Qux,Qxu,Quu_hat,Qux_hat] = Q_info(rbt,cst,x,u,Vx,Vxx,params,iter)
%Q_INFO
V_reg = 0.00;
if isfield(params,'shooting_phase') 
   V_reg = 0.00 + 1.0 / exp(iter);
end
nx = params.nx;
Vxx_hat = Vxx + V_reg * eye(nx) * (params.Reg_Type == 2);
[fx,fu] = rbt.getLinSys(x,u, params.dt);
[~,lx,lu,lxx,lux,lxu,luu] = cst.l_info(x,u);
Qx = lx + fx.' * Vx;
Qu = lu + fu.' * Vx;
Qxx = lxx + fx.' * Vxx * fx;
Quu = luu + fu.' * Vxx * fu;
Qux = lux + fu.' * Vxx * fx;
Qxu = lxu + fx.' * Vxx * fu;

Quu_hat = luu + fu.' * Vxx_hat * fu;
Qux_hat = lux + fu.' * Vxx_hat * fx;
end
