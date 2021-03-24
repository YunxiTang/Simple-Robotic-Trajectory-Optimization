function [Qx,Qu,Qxx,Quu,Qux,Qxu] = Qcms_info(rbt,cst,constraint,lambda,Imu,x,u,Vx,Vxx,params)
%Qms_INFO : For constrained multiple shooting DDP/SLQ
nx = params.nx;
Vxx = Vxx + 0.1 * eye(nx) * (params.Reg_Type == 2);

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
end
