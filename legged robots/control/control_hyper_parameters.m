% You can set any hyper parameters of the control function here; you may or
% may not want to use the step_number as the input of the function. 
function [params] = control_hyper_parameters(step_number)
clearvars params
% try
%     load('SAVES/params.mat', 'params');
% catch
%     fprintf('\ncould not load params, using default.\n');
    params.sw_target = 0.15166;
    params.kp_t = 500; 
    params.kd_t = 20;
    params.kp_s = 900; 
    params.kd_s = 30;
    params.alfa = 16.29; 
    params.t_target = 0.14;
    params.f = 6; %%Hz

% end
s = 20;
if step_number > s; step_number = s; end
params.t_target = params.t_target + 0.06 / s * step_number;
% params.sw_target = params.sw_target - 0.01 / s * step_number; 



% if step_number > s; step_number = s; end
% params.t_target = params.t_target / s * step_number;

end


%     params.sw_target = 0.14166;
%     params.kp_t = 700; 
%     params.kd_t = 20;
%     params.kp_s = 700; 
%     params.kd_s = 30;
%     params.alfa = 16.29; 
%     params.t_target = 0.11;
%     params.f = 6; %%Hz
