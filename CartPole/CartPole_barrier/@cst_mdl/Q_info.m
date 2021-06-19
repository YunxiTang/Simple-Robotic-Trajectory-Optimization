function [Qx,Qu,Qxx,Quu,Qux,Qxu] = Q_info(rbt,cst,x,u,Vx,Vxx,params,iter)
%Q_INFO
V_reg = 0.00;
if isfield(params,'shooting_phase') 
   V_reg = 0.00 + 1.0 / exp(iter);
end
nx = params.nx;
Vxx = Vxx + V_reg * eye(nx) * (params.Reg_Type == 2);
if params.full_ddp == 1
    [fx,fu,fxx,fuu,fux,fxu] = rbt.getLinSys(x,u, params.dt);
else
    [fx,fu] = rbt.getLinSys(x,u, params.dt);
    fxx = zeros(params.nx, params.nx, params.nx);
    fuu = zeros(params.nu, params.nu, params.nu);
    fxu = zeros(params.nx, params.nx, params.nu);
    fux = zeros(params.nu, params.nu, params.nx);
end
[~,lx,lu,lxx,lux,lxu,luu] = cst.l_info(x,u);
Qx = lx + fx.' * Vx;
Qu = lu + fu.' * Vx;
Qxx = lxx + fx.' * Vxx * fx + Vx * fxx;
Quu = luu + fu.' * Vxx * fu + Vx * fuu;
Qux = lux + fu.' * Vxx * fx + Vx * fux;
Qxu = lxu + fx.' * Vxx * fu + Vx * fxu;
end
