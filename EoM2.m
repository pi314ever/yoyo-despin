function [dx] = EoM2(t,x)
%Equation of Motion for yo-yo despin Stage 2
%   Stage 2: String length constant, rotating about single pivot point.
%   INDEX   x    y    z    gamma phi
%   f(t)    1    3    5    7    9    
%   f'(t)   2    4    6    8    10    
%   x,y,z are in Newtonian frame
%   gamma and phi are with respect to the spin-independent body frame (B)

global l r mr mw G me
dx = ones(10,1); % Initialize output

for ii = 1:5
    dx(2*ii-1) = x(2*ii);
end

R = sqrt(x(1).^2+x(3).^2+x(5).^2);
I = 1/2*mr*r^2;
A = cos(x(8))/l.*x(10).^2*r;
B = -(sin(x(7))*r/l+1);

dx(2) = -G*me*x(1)./R.^3;
dx(4) = -G*me*x(3)./R.^3;
dx(6) = -G*me*x(5)./R.^3;
dx(10) = (-mw*l^2*A-l*sin(x(7)).*A-l*cos(x(7)).*(x(8)+2*x(10)).*x(8)*r)./...
    (I+mw*l^2*(B+1)+mw*r^2+l*sin(x(7)).*(B+2)*r);
dx(8) = A+B.*dx(10);

end

