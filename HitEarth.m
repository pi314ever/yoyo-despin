function [position,isterminal,direction] = HitEarth(t,x)
%Free Fall Halt event function. 
%   Stop Free Fall simulation when object hits surface of Earth or is too
%   far away from Earth
global Re
dist = sqrt(x(1).^2+x(3).^2+x(5).^2);
% Terminate when position hits Earth or is too far away
position = (Re-dist).*(Re*3-dist); 
isterminal = 1; % Terminate simulation
direction = 0; 
end

