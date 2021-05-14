function [Qx,Qu,Qxx,Quu,Qux,Qxu,Quu_hat,Qux_hat] = Qcms_info(rbt,cst,constraint,lambda,Imu,x,u,Vx,Vxx,params,iter)
%Qms_INFO : For constrained multiple shooting DDP/SLQ

reg = 1.0 / exp(iter);
nx = params.nx;
Vxx_hat = Vxx + reg * eye(nx) * (params.Reg_Type == 2);

[fx,fu] = rbt.getLinSys(x,u);
[~,lx,lu,lxx,lux,lxu,luu] = cst.l_info(x,u);
c = constraint.c_ineq(x, u);
[cx,cu] = constraint.algrad(x,u);

Qx = lx + fx.' * (Vx) + cx' * (lambda + Imu*c)*params.dt;
Qu = lu + fu.' * (Vx) + cu' * (lambda + Imu*c)*params.dt;
Qxx = lxx + fx.' * Vxx * fx + cx'*Imu*cx*params.dt;
Quu = luu + fu.' * Vxx * fu + cu'*Imu*cu*params.dt;
Qux = lux + fu.' * Vxx * fx + cu'*Imu*cx*params.dt;
Qxu = lxu + fx.' * Vxx * fu + cx'*Imu*cu*params.dt;

Quu_hat = luu + fu.' * Vxx_hat * fu + cu'*Imu*cu*params.dt;
Qux_hat = lux + fu.' * Vxx_hat * fx + cu'*Imu*cx*params.dt;
end
