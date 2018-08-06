function zdot = relDyns(t,z,p)
% z = [r; rd; h; th; I; Omg; x1; y1; z1; x1d; y1d; z1d; ...]

% unpack
kJ2 = p.kJ2;
mu = p.mu;
r = z(1);
rd = z(2);
h = z(3);
th = z(4);
I = z(5);
Omg = z(6);

% LVLH frame dynamics
zdot = [rd;...
    (-mu/(r^2) + (h^2)/(r^3) - kJ2*(1 - 3*(sin(I)^2)*(sin(th)^2))/(r^4));...
    -kJ2*(sin(I)^2)*sin(2*th)/(r^3);...
    (h/(r^2) + 2*kJ2*(cos(I)^2)*(sin(th)^2)/(h*r^3));...
    -kJ2*sin(2*I)*sin(2*th)/(2*h*r^3);...
    -2*kJ2*cos(I)*(sin(th)^2)/(h*r^3)];

% relative dynamics
for j = 7:6:length(z)
    % relative positions and velocities
    xj = z(j); yj = z(j+1); zj = z(j+2);
    xdj = z(j+3); ydj = z(j+4); zdj = z(j+5);
    
    % variables
    wz = h/r^2;
    wx = -kJ2*sin(2*I)*sin(th)/(h*r^3);
    rj = norm([(r+xj); yj; zj]);
    rjZ = (r+xj)*sin(I)*sin(th) + yj*sin(I)*cos(th) + zj*cos(I);
    xi = 2*kJ2*sin(I)*sin(th)/r^4;
    xij = 2*kJ2*rjZ/rj^5;
    eta2 = mu/r^3 + kJ2/r^5 - 5*kJ2*(sin(I)^2)*(sin(th)^2)/r^5;
    eta2j = mu/rj^3 + kJ2/rj^5 - 5*kJ2*(rjZ^2)/rj^7;
    az = -2*h*rd/r^3 - kJ2*(sin(I)^2)*sin(2*th)/r^5;
    ax = -kJ2*sin(2*I)*cos(th)/r^5 + 3*rd*kJ2*sin(2*I)*sin(th)/(h*r^4) ...
        - 8*(kJ2^2)*(sin(I)^3)*cos(I)*(sin(th)^2)*cos(th)/((r^6)*h^2);
    
    zdot(j:(j+2)) = [xdj; ydj; zdj];
    zdot((j+3):(j+5)) = [2*ydj*wz - xj*(eta2j - wz^2) + yj*az - zj*wx*wz - (xij - xi)*sin(I)*sin(th) - r*(eta2j - eta2);...
        -2*xdj*wz + 2*zdj*wx - xj*az - yj*(eta2j- wz^2 - wx^2) + zj*ax - (xij - xi)*sin(I)*cos(th);...
        -2*ydj*wx - xj*wx*wz - yj*ax - zj*(eta2j - wx^2) - (xij - xi)*cos(I)];
end
end