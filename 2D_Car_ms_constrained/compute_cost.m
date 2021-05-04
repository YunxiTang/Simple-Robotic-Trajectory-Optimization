function [J] = compute_cost(xsol,usol,cst_mdl,params)
%COMPUTE_COST
Nt = size(xsol,2);
J_path = 0.0;

for i=1:Nt-1
    xi = xsol(:,i);
    ui = usol(:,i);
    J_path = J_path + cst_mdl.l_cost(xi,ui);
end
J = J_path + cst_mdl.lf_cost(xsol(:,end));
end

