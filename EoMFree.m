function [dx] = EoMFree(t,x)
%Equations of Motion for free fall (translational motion only)
%   Single mass free-fall equations of motion
%   INDEX   x    y    z   
%   f(t)    1    3    5   
%   f'(t)   2    4    6    

global G me
dx = ones(6,1); % Initialize output

for ii = 1:3
    dx(2*ii-1) = x(2*ii);
end

R = sqrt(x(1).^2+x(3).^2+x(5).^2);

dx(2) = -G*me*x(1)./R.^3;
dx(4) = -G*me*x(3)./R.^3;
dx(6) = -G*me*x(5)./R.^3;
end

