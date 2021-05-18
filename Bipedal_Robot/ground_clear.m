function [z_swf] = ground_clear(x, rbt)
%GROUND_CLEAR 
q = x(1:3);
[~, z_swf, ~, ~] = kin_swf(q, rbt);
end

