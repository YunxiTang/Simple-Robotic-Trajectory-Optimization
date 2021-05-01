function [J] = compute_cost(tsol,xsol,usol,params)
%COMPUTE_COST 
    dt = tsol(2)-tsol(1);
    t0 = tsol(1);
    x0 = xsol(:,1);
    tf = tsol(end);
    xf = xsol(:,end);
    J_path = sum(pathObjFunc(tsol(1:end-1),xsol(:,1:end-1),usol(:,1:end-1),params) * dt);
    J_bnd = bndObjFunc(t0,x0,tf,xf,params);
    J = J_path + J_bnd;
end

