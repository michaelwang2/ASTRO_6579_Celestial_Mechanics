function zdot = relDynsSmallr(t,z,p)
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
    eta2 = mu/r^3 + kJ2/r^5 - 5*kJ2*(sin(I)^2)*(sin(th)^2)/r^5;
    az = -2*h*rd/r^3 - kJ2*(sin(I)^2)*sin(2*th)/r^5;
    ax = -kJ2*sin(2*I)*cos(th)/r^5 + 3*rd*kJ2*sin(2*I)*sin(th)/(h*r^4) ...
        - 8*(kJ2^2)*(sin(I)^3)*cos(I)*(sin(th)^2)*cos(th)/((r^6)*h^2);
    
    % relative velocities
    zdot(j:(j+2)) = [xdj; ydj; zdj];
    
    % relative accelerations
    A1 = [0, 2*wz, 0;...
        -2*wz, 0, 2*wx;...
        0, -2*wx, 0];
    A2 = [(2*eta2 + wz^2 + 2*kJ2*(1 - (sin(I)^2)*(sin(th)^2))/r^5), (az + 4*kJ2*(sin(I)^2)*sin(2*th)/r^5), -5*wx*wz;...
        (4*kJ2*(sin(I)^2)*sin(2*th)/r^5 - az), (2*kJ2*(sin(I)^2)*(cos(th)^2)/r^5 + eta2 - wz^2 - wx^2), (ax - kJ2*sin(2*I)*cos(th)/r^5);...
        (-5*wx*wz), (kJ2*sin(2*I)*cos(th)/r^5 + ax), (eta2 - wx^2 + 2*kJ2*(cos(I)^2)/r^5)];
    zdot((j+3):(j+5)) = A1*[xdj; ydj; zdj] + A2*[xj; yj; zj];
end
end