function [Qx,Qu,Qxx,Quu,Qux,Qxu] = Q_info(rbt,cst,x,u,Vx,Vxx,params,iter)
%Q_INFO
V_reg = 0.00;
if isfield(params,'shooting_phase') && params.shooting_phase > 1
   V_reg = 0.00 + 0.3 / exp(iter);
end
nx = params.nx;
Vxx = Vxx + V_reg * eye(nx) * (params.Reg_Type == 2);
[fx,fu] = rbt.getLinSys(x,u,params.dt);
[~,lx,lu,lxx,lux,lxu,luu] = cst.l_info(x,u);
Qx = lx + fx.' * Vx;
Qu = lu + fu.' * Vx;
Qxx = lxx + fx.' * Vxx * fx;
Quu = luu + fu.' * Vxx * fu;
Qux = lux + fu.' * Vxx * fx;
Qxu = lxu + fx.' * Vxx * fu;
end
