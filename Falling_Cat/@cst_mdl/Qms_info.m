function [Qx,Qu,Qxx,Quu,Qux,Qxu] = Qms_info(rbt,cst,x,u,Vx,Vxx,params,dft)
%Qms_INFO : For multiple shooting DDP/SLQ
if nargin == 7
    nx = params.nx;
    Vxx = Vxx + 0.1 * eye(nx) * (params.Reg_Type == 2);
    dft = zeros(params.nx, 1);
end
[fx,fu] = rbt.getLinSys(x,u);
[~,lx,lu,lxx,lux,lxu,luu] = cst.l_info(x,u);
Qx = lx + fx.' * (Vx + Vxx * dft);
Qu = lu + fu.' * (Vx + Vxx * dft);
Qxx = lxx + fx.' * Vxx * fx;
Quu = luu + fu.' * Vxx * fu;
Qux = lux + fu.' * Vxx * fx;
Qxu = lxu + fx.' * Vxx * fu;
end
