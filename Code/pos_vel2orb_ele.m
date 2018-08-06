function [a, e, Omega, I, w, nu] = pos_vel2orb_ele(r, vI, mu)
% this function converts position and velocity to orbital elements
a = 0;
e = 0;
Omega = 0;
I = 0;
w = 0;
nu = 0;

% inertial basis vectors
e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1];

% specific angular momentum of P
hP = cross(r, vI);

% eccentricity vector
e_vec = cross(vI, hP)./mu - r./norm(r);
e = norm(e_vec);
e_vec = e_vec./e;

% semi-major axis
a = (norm(hP)^2)/(mu*(1 - e^2));

% true anomaly
nu = atan2(sqrt(e^2 - ((norm(hP)^2)/(mu*norm(r)) - 1)^2), (norm(hP)^2)/(mu*norm(r)) - 1);

% line of nodes
n = cross(e3, hP);
if sum(n) < 1e-9
    w = atan2(e_vec(2), e_vec(1));
   return; 
end
n = n./norm(n);

% Longitude of ascending node
Omega = atan2(norm(n - dot(n, e1).*e1), dot(n, e1));

% Inclination
I = atan2(norm(hP - dot(hP, e3).*e3), dot(hP, e3));

% Argument of periapsis
w = atan2(norm(e_vec - dot(e_vec, n).*n), dot(e_vec, n));
end