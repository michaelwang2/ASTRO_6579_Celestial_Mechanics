function zdot = refOrbit(t,z,p)
% z = [r; rd; h; th; I; Omg]

% unpack
kJ2 = p.kJ2;
mu = p.mu;
r = z(1);
rd = z(2);
h = z(3);
th = z(4);
I = z(5);
Omg = z(6);

% dynamics
zdot = [rd;...
    (-mu/(r^2) + (h^2)/(r^3) - kJ2*(1 - 3*(sin(I)^2)*(sin(th)^2))/(r^4));...
    -kJ2*(sin(I)^2)*sin(2*th)/(r^3);...
    (h/(r^2) + 2*kJ2*(cos(I)^2)*(sin(th)^2)/(h*r^3));...
    -kJ2*sin(2*I)*sin(2*th)/(2*h*r^3);...
    -2*kJ2*cos(I)*(sin(th)^2)/(h*r^3)];
end