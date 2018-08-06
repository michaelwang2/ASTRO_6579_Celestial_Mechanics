function err = errNorm_lsqnonlin(z0)
% z0 = [x0; T]
% constants
p.mu = 398600.4418; % km^3/s^2
p.J2 = 1.08262668e-3;
p.Re = 6.3781e3; % km
p.kJ2 = (3*(p.J2)*(p.mu)*(p.Re)^2)/2; % km^5/s^2

% Full Nonlinear Relative Dynamics
opts.RelTol = 5e-12;
opts.AbsTol = 5e-12;
[t, zarray] = ode45(@relDyns, [0, z0(end)], z0(1:(end-1)), opts, p);
t = t./3600;
r = zarray(:, 1);
rd = zarray(:, 2);
h = zarray(:, 3);
th = zarray(:, 4);
I = zarray(:, 5);
Omg = zarray(:, 6);
[energy, state, center] = LVLH_to_ECI(zarray, 1, p);

err = (zarray(1,:) - zarray(end, :));
end