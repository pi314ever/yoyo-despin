function [dx] = EoM1(t,x)
%Equation of Motion for Yo-yo despin maneuver stage 1. 
%   Stage 1: Yo-yo string length variable, always tangent to rocket body
%   INDEX   x    y    z    beta phi
%   f(t)    1    3    5    7    9    
%   f'(t)   2    4    6    8    10    
%   x,y,z are in Newtonian frame
%   beta and phi are with reference to the spin-independent body frame (B)

global mr mw G me
dx = ones(10,1); % Initialize output

for ii = 1:5
    dx(2*ii-1) = x(2*ii);
end

R = sqrt(x(1).^2+x(3).^2+x(5).^2);

dx(2) = -G*me*x(1)./R.^3;
dx(4) = -G*me*x(3)./R.^3;
dx(6) = -G*me*x(5)./R.^3;
dx(8) = (x(8).^2-2*x(8).*(x(8)-x(10)))./(x(7)-x(9)); % Beta dd
dx(10) = -mw*x(8).^2.*(x(7)-x(9))./(0.5*mr+mw); % Phi dd
end

