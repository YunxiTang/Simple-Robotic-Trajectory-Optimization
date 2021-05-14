function [Qx,Qu,Qxx,Quu,Qux,Qxu,Quu_hat,Qux_hat] = Qcs_info(rbt,cst,x,u,xref,uref,Vx,Vxx,params,path_constraint)
%Qcs_INFO : For multiple shooting DDP/SLQ
if nargin > 7
    nx = params.nx;
    Vxx_hat = Vxx + 0.01 * eye(nx) * (params.Reg_Type == 2);
end
[fx,fu] = rbt.getLinSys(x,u,params.dt);
[~, Penal_x, Penal_u, Penal_xu, Penal_ux, Penal_xx, Penal_uu] = path_constraint.penalty_info(x, u);
[~,lx,lu,lxx,lux,lxu,luu] = cst.l_info(x,u,xref,uref);
Qx = lx + Penal_x + fx.' * (Vx);
Qu = lu + Penal_u + fu.' * (Vx);
Qxx = lxx + Penal_xx + fx.' * Vxx * fx;
Quu = luu + Penal_uu + fu.' * Vxx * fu;
Qux = lux + Penal_ux + fu.' * Vxx * fx;
Qxu = lxu + Penal_xu + fx.' * Vxx * fu;
if nargin > 7
    Quu_hat = luu + Penal_uu + fu.' * Vxx_hat * fu;
    Qux_hat = lux + Penal_ux + fu.' * Vxx_hat * fx;
else
    Quu_hat = Quu;
    Qux_hat = Qux;
end
end
