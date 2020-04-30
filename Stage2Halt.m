function [position,isterminal,direction] = Stage2Halt(t,x)
%Stage 2 Halt event function. 
%   Stop Stage 2 simulation when rope reaches 90 degrees (pi/2)

position = x(7)-pi/2; % Terminate at gamma = pi/2
isterminal = 1; % Terminate simulation
direction = 0; 
end


