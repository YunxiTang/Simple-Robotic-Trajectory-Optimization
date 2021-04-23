function [Qx,Qu,Qxx,Quu,Qux,Qxu] = Qcss_info(rbt,cst,constraint,lambda,Imu,x,u,Vx,Vxx,params)
%Qcss_INFO
    if nargin > 9
        nx = params.nx;
        Vxx = Vxx + 0.1 * eye(nx) * (params.Reg_Type == 2);
    end
    [fx,fu] = rbt.getLinSys(x,u);
    [~,lx,lu,lxx,lux,lxu,luu] = cst.l_info(x,u);
    c = constraint.c_ineq(x, u);
    [cx,cu] = constraint.algrad(x,u);
    Qx = lx + fx.' * Vx + cx' * (lambda + Imu*c);
    Qu = lu + fu.' * Vx + cu' * (lambda + Imu*c);
    Qxx = lxx + fx.' * Vxx * fx + cx'*Imu*cx;
    Quu = luu + fu.' * Vxx * fu + cu'*Imu*cu;
    Qux = lux + fu.' * Vxx * fx + cu'*Imu*cx;
    Qxu = lxu + fx.' * Vxx * fu + cx'*Imu*cu;
end
