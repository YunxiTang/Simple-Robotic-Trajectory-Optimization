%%
% This function animates the solution of the equations of motion of the
% three link biped. 
% sln is the solution computed by solve_eqns.m
%%
function rbt_animate(sln, rbt)

figure(963);
skip = 20;
tic();
num_steps = length(sln.T);
r0 = [0; 0];
for j = 1:num_steps
    Y = sln.Y{j};
    [~, N] = size(Y);
    for i = 1:skip:N
        q = Y(1:3, i);
        pause(0.002);
        visualize(q, r0, j, rbt);
        hold on
    end
    [x0, ~, ~, ~] = kin_swf(q, rbt);
    r0 = r0 + [x0; 0];
end
t_anim = toc();
real_time_factor = sln.T{end}(end) / t_anim;
fprintf('Real time factor:');
disp(real_time_factor);
end