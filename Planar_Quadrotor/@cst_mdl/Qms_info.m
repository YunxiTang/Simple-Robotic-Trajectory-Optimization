function [Qx,Qu,Qxx,Quu,Qux,Qxu,Quu_hat,Qux_hat] = Qms_info(rbt,cst,x,u,Vx,Vxx,params,dft)
%Qms_INFO : For multiple shooting DDP/SLQ
nx = params.nx;
Vxx_hat = Vxx + 0.1 * eye(nx) * (params.Reg_Type == 2);
[fx,fu] = rbt.getLinSys(x,u);
[~,lx,lu,lxx,lux,lxu,luu] = cst.l_info(x,u);
Qx = lx + fx.' * (Vx + Vxx * dft);
Qu = lu + fu.' * (Vx + Vxx * dft);
Qxx = lxx + fx.' * Vxx * fx;
Quu = luu + fu.' * Vxx * fu;
Qux = lux + fu.' * Vxx * fx;
Qxu = lxu + fx.' * Vxx * fu;

Quu_hat = luu + fu.' * Vxx_hat * fu;
Qux_hat = lux + fu.' * Vxx_hat * fx;
end
