%%
% This function takes the configuration of the 3-link model and plots the 
% 3-link model. 
% q = [q1, q2 ,q3] the generalized coordinates. Try different angles to see
% if your formulas for x1, z1, etc. makes sense. Example: q = [-pi/6, pi/6,
% pi/8]
% r0 is the position of the stance foot in the global frame. 
%%
function visualize(q, r0, step_number)

    % default r0 = [0; 0]
    if nargin == 1
        r0 = [0; 0];
        step_number = -1;
    end
    
    if nargin == 2
        step_number = -1;
    end

    x0 = r0(1);
    z0 = r0(2);
    
    [~, ~, ~, l1, l2, l3, ~] = set_parameters;
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    
    x1 = (l1*sin(q1))/2 + x0;
    z1 = (l1*cos(q1))/2 + z0;
    x2 = l1*sin(q1) - (l2*sin(q2))/2 + x0;
    z2 = l1*cos(q1) - (l2*cos(q2))/2 + z0;
    x3 = l1*sin(q1) + (l3*sin(q3))/2 + x0;
    z3 = l1*cos(q1) + (l3*cos(q3))/2 + z0;
    
    x_h = l1*sin(q1) + x0;
    z_h = l1*cos(q1) + z0;
    
    x_t = l1*sin(q1) + l3*sin(q3) + x0;
    z_t = l1*cos(q1) + l3*cos(q3) + z0;
    
    x_swf = l1*sin(q1) - l2*sin(q2) + x0;
    z_swf = l1*cos(q1) - l2*cos(q2) + z0;
    %% External perturbation 
    [ext_step_number, value] = ext_perturbation_parameters();
    if step_number == ext_step_number
        F = value / 100;
    else
        F = 0;
    end
    %% 
    % Here plot a schematic of the configuration of three link biped at the
    % generalized coordinate q = [q1, q2, q3]:
    lw = 4;
    % links
    plot([x0, x_h], [z0, z_h], 'linewidth', lw); 
    hold on
    plot([x_h, x_t], [z_h, z_t], 'linewidth', lw); 
    plot([x_h, x_swf], [z_h, z_swf], 'linewidth', lw);
    % plot a line for "ground"
    plot([-1 + x_h, 1 + x_h], [0, 0], 'color', 'black','linewidth', 1);
    axis 'square'
    xlim([-1 + x_h, 1 + x_h]);
    ylim([-0.8, 1.2]);
    % visualize external perturbation
    delete(findall(gcf,'type','annotation'));
    if step_number == ext_step_number
        h = annotation('arrow', 'X', [1 / 2 - 0.1 * F, 1 / 2 + 0.02], 'Y', [(z_h + 0.8) / 2, (z_h + 0.8) / 2]);
        h.Color = 'red';
    end
    plot([x_h - F / 200, x_h], [z_h, z_h], 'linewidth', 2 * lw, 'color', 'r')
    % point masses
    mz = 40;
    plot(x1, z1, '.', 'markersize', mz); 
    hold on
    plot(x2, z2, '.', 'markersize', mz); 
    plot(x3, z3, '.', 'markersize', mz);
end