function [Qx,Qu,Qxx,Quu,Qux,Qxu] = Q_info(rbt,cst,x,u,Vx,Vxx,params)
%Q_INFO
    if nargin > 6
        nx = params.nx;
        Vxx = Vxx + 0.1 * eye(nx) * (params.Reg_Type == 2);
    end
    [fx,fu] = rbt.getLinSys(x,u);
    [~,lx,lu,lxx,lux,lxu,luu] = cst.l_info(x,u);
    Qx = lx + fx.' * Vx;
    Qu = lu + fu.' * Vx;
    Qxx = lxx + fx.' * Vxx * fx;
    Quu = luu + fu.' * Vxx * fu;
    Qux = lux + fu.' * Vxx * fx;
    Qxu = lxu + fx.' * Vxx * fu;
end
