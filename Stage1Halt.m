function [position,isterminal,direction] = Stage1Halt(t,x)
%Stage 1 Halt event function. 
%   Stop Stage 1 simulation when length of rod reaches max length
global r l
position = r*(x(7)-x(9))-l; % Terminate at length = length_max
isterminal = 1; % Terminate simulation
direction = 0; 
end

