function [J] = compute_cost(xsol,usol,cst_mdl,params)
%COMPUTE_COST
Nt = size(xsol,2);
J_path = 0.0;
xref = params.xf;
uref = zeros(params.nu,1);
for i=1:Nt-1
    xi = xsol(:,i);
    ui = usol(:,i);
    J_path = J_path + cst_mdl.l_cost(xi,ui,xref,uref);
end
J = J_path + cst_mdl.lf_cost(xsol(:,end));
end

