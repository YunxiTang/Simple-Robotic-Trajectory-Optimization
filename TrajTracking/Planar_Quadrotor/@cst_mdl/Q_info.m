function [Qx,Qu,Qxx,Quu,Qux,Qxu,Quu_hat,Qux_hat] = Q_info(rbt,cst,x,u,xref,uref,Vx,Vxx,params)
%Q_INFO
if nargin > 6
    nx = params.nx;
    Vxx_hat = Vxx + 0.01 * eye(nx) * (params.Reg_Type == 2);
end
[fx,fu] = rbt.getLinSys(x,u);
[~,lx,lu,lxx,lux,lxu,luu] = cst.l_info(x,u,xref,uref);
Qx = lx + fx.' * Vx;
Qu = lu + fu.' * Vx;
Qxx = lxx + fx.' * Vxx * fx;
Quu = luu + fu.' * Vxx * fu;
Qux = lux + fu.' * Vxx * fx;
Qxu = lxu + fx.' * Vxx * fu;
if nargin > 6
    Quu_hat = luu + fu.' * Vxx_hat * fu;
    Qux_hat = lux + fu.' * Vxx_hat * fx;
else
    Quu_hat = Quu;
    Qux_hat = Qux;
end
end
