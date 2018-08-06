function [energy, state, center] = LVLH_to_ECI(zarray, numSat, p)
r = zarray(:, 1);
rd = zarray(:, 2);
h = zarray(:, 3);
th = zarray(:, 4);
I = zarray(:, 5);
Omg = zarray(:, 6);

% compute energy, state, and center
energy = zeros(length(r), numSat);
state = zeros(length(r), 6*numSat);
center = zeros(length(r), 6);
for i = 1:length(r)
    R_Omg = [cos(Omg(i)), -sin(Omg(i)), 0;...
        sin(Omg(i)), cos(Omg(i)), 0;...
        0, 0, 1];
    R_I = [1, 0, 0;...
        0, cos(I(i)), -sin(I(i));...
        0, sin(I(i)), cos(I(i))];
    R_th = [cos(th(i)), -sin(th(i)), 0;...
        sin(th(i)), cos(th(i)), 0;...
        0, 0, 1];
    I_R_LVLH = R_Omg*R_I*R_th;
    
    % angular velocities
    wz = h(i)/r(i)^2;
    wx = -p.kJ2*sin(2*I(i))*sin(th(i))/(h(i)*r(i)^3);
    
    % center of LVLH frame
    center(i, 1:3) = (I_R_LVLH*[r(i); 0; 0])';
    vI = I_R_LVLH*[rd(i); h(i)/r(i); 0];
    center(i, 4:6) = vI';
    
    for j = 1:numSat
        % relative positions and velocities
        xj = zarray(i, 6*j+1); yj = zarray(i, 6*j+2); zj = zarray(i, 6*j+3);
        xdj = zarray(i, 6*j+4); ydj = zarray(i, 6*j+5); zdj = zarray(i, 6*j+6);
        rj = norm([(r(i) + xj), yj, zj]);
        rjZ = (r(i) + xj)*sin(I(i))*sin(th(i)) + yj*sin(I(i))*cos(th(i)) + zj*cos(I(i));
        
        % position
        state(i, (6*j-5):(6*j-3)) = (I_R_LVLH*[(xj + r(i)); yj; zj])';

        % velocity
        vjI = I_R_LVLH*[(xdj + rd(i) - yj*wz);...
            (ydj + (r(i) + xj)*wz - zj*wx);...
            (zdj + yj*wx)];
        state(i, (6*j-2):(6*j)) = vjI';
        
        % total energy
        energy(i, j) = 0.5*(vjI')*vjI + (-p.mu/rj - p.kJ2/(3*rj^3) + p.kJ2*(rjZ^2)/(rj^5));
    end
end
end