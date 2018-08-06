%% Relative Dynamics
clear; close; clc;

syms th I Omg thd Id Omgd  real

% calculate angular velocity vector
R_Omg = [cos(Omg), -sin(Omg), 0;...
    sin(Omg), cos(Omg), 0;...
    0, 0, 1];
R_I = [1, 0, 0;...
    0, cos(I), -sin(I);...
    0, sin(I), cos(I)];
R_th = [cos(th), -sin(th), 0;...
    sin(th), cos(th), 0;...
    0, 0, 1];
I_R_LVLH = simplify(R_Omg*R_I*R_th);
w = R_th.'*R_I.'*[0;0;Omgd] + R_th.'*[Id;0;0] + [0;0;thd];

% find gradient of U
syms mu r X Y Z kJ2   real
DU_ECI = [(mu*X/r^3 + kJ2*X/r^5 - 5*kJ2*(Z^2)*X/r^7);...
    (mu*Y/r^3 + kJ2*Y/r^5 - 5*kJ2*(Z^2)*Y/r^7);...
    (mu*Z/r^3 + kJ2*Z/r^5 + kJ2*(2*Z*r^2 - 5*Z^3)/r^7)];% in ECI frame
DU_LVLH = simplify(expand(I_R_LVLH'*DU_ECI));
DU_LVLH = subs(DU_LVLH, Z^2 / r^2, sin(th)^2 * sin(I)^2);

%% Simulate Nonlinear Relative Dynamic
clear; close; clc;

% constants
p.mu = 398600.4418; % km^3/s^2
p.J2 = 1.08262668e-3;
p.Re = 6.3781e3; % km
p.kJ2 = (3*(p.J2)*(p.mu)*(p.Re)^2)/2; % km^5/s^2

% initial conditions
load('ImageSat_Inertial_Position_Velocity3.mat');
x = data(:, 1);
y = data(:, 2);
z = data(:, 3);
vx = data(:, 4);
vy = data(:, 5);
vz = data(:, 6);
start = 1;
rvec = [x(start); y(start); z(start)]; % km
vvec = [vx(start); vy(start); vz(start)]; % km/s
% [a0, e0, Omega0, I0, w0, nu0] = pos_vel2orb_ele(rvec, vvec, p.mu);
[a0 ,e0 ,E0 ,I0 , w0 , Omega0 ,P , tau ,A , B] = vec2orbElem (rvec, vvec, p.mu);
nu0 = 2*atan2(sqrt((1+e0)/(1-e0))*tan(E0/2), 1);
% z = [r; rd; h; th; I; Omg]
r0 = norm(rvec);
rd0 = vvec'*(rvec./norm(rvec));
h0 = norm(cross(rvec, vvec));
th0 = w0 + nu0;
z0 = [r0; rd0; h0; th0; I0; Omega0];

% Convert from cartesian to hybrid
state = zeros(length(time(start:end)), 7);
for i = start:length(time)
   rvec = [x(i); y(i); z(i)]; % km
   vvec = [vx(i); vy(i); vz(i)]; % km/s
   r = norm(rvec);
   rd = vvec'*(rvec./norm(rvec));
   h = norm(cross(rvec, vvec));
   % [a, e, Omega, I, w, nu] = pos_vel2orb_ele(rvec, vvec, p.mu);
   [a ,e ,E ,I , w , Omega ,P , tau ,A , B] = vec2orbElem (rvec, vvec, p.mu);
   nu = 2*atan2(sqrt((1+e)/(1-e))*tan(E/2), 1);
   th = w + nu;
   state(i-start+1, :) = [time(i).*60, r, rd, h, th, I, Omega];
end

% integration
opts.RelTol = 1e-13;
opts.AbsTol = 1e-13;
[t, zarray] = ode45(@refOrbit, time(start:end).*60, z0, opts, p);
r = zarray(:, 1);
rd = zarray(:, 2);
h = zarray(:, 3);
th = zarray(:, 4);
I = zarray(:, 5);
Omg = zarray(:, 6);

% plot
figure;
subplot(3,2,1); hold on;
plot(t, r-state(:, 2));
ylabel('$$r\ [km]$$', 'interpreter', 'latex');
grid on; box on;

subplot(3,2,2); hold on;
plot(t, rd-state(:, 3));
ylabel('$$\dot{r}\ [km/s]$$', 'interpreter', 'latex');
grid on; box on;

subplot(3,2,3); hold on;
plot(t, h-state(:, 4));
ylabel('$$h\ [km^2/s]$$', 'interpreter', 'latex');
grid on; box on;

subplot(3,2,4); hold on;
plot(t, rad2deg(wrapToPi(th))-rad2deg(state(:, 5)));
ylabel('$$\theta\ [^\circ]$$', 'interpreter', 'latex');
grid on; box on;

subplot(3,2,5); hold on;
plot(t, rad2deg(wrapToPi(I))-rad2deg(state(:, 6)));
ylabel('$$I\ [^\circ]$$', 'interpreter', 'latex');
xlabel('$$Time\ [s]$$', 'interpreter', 'latex');
grid on; box on;

subplot(3,2,6); hold on;
plot(t, rad2deg(wrapToPi(Omg))-rad2deg(state(:, 7)));
ylabel('$$\Omega\ [^\circ]$$', 'interpreter', 'latex');
xlabel('$$Time\ [s]$$', 'interpreter', 'latex');
grid on; box on;

%% Convert Hybrid elements to ECI
clear; close; clc;

% constants
p.mu = 398600.4418; % km^3/s^2
p.J2 = 1.08262668e-3;
p.Re = 6.3781e3; % km
p.kJ2 = (3*(p.J2)*(p.mu)*(p.Re)^2)/2; % km^5/s^2

% initial conditions
load('ImageSat_Inertial_Position_Velocity.mat');
x = data(:, 1);
y = data(:, 2);
z = data(:, 3);
vx = data(:, 4);
vy = data(:, 5);
vz = data(:, 6);
rvec = [x(1); y(1); z(1)]; % km
vvec = [vx(1); vy(1); vz(1)]; % km/s
% [a0, e0, Omega0, I0, w0, nu0] = pos_vel2orb_ele(rvec, vvec, p.mu);
[a0 ,e0 ,E0 ,I0 , w0 , Omega0 ,P , tau ,A , B] = vec2orbElem (rvec, vvec, p.mu);
nu0 = 2*atan2(sqrt((1+e0)/(1-e0))*tan(E0/2), 1);
% z = [r; rd; h; th; I; Omg]
r0 = norm(rvec);
rd0 = vvec'*(rvec./norm(rvec));
h0 = norm(cross(rvec, vvec));
th0 = w0 + nu0;
z0 = [r0; rd0; h0; th0; I0; Omega0];

% integration
opts.RelTol = 5e-14;
opts.AbsTol = 5e-14;
[t, zarray] = ode45(@refOrbit, time.*60, z0, opts, p);
t = t./60;
r = zarray(:, 1);
rd = zarray(:, 2);
h = zarray(:, 3);
th = zarray(:, 4);
I = zarray(:, 5);
Omg = zarray(:, 6);

% convert from LVLH to ECI
state = zeros(length(r), 6);
energy = zeros(length(r), 1);
energySTK = zeros(length(r), 1);
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
    
    % position
    state(i, 1:3) = (I_R_LVLH*[r(i); 0; 0])';
    
    % velocity
    vI = I_R_LVLH*[rd(i); h(i)/r(i); 0];
    state(i, 4:6) = (vI)';
    
    % energy
    energy(i) = 0.5*vI'*vI;
    energySTK(i) = 0.5*norm([vx(i), vy(i), vz(i)])^2;
end

% plot
spacing = -0.18;
figure;
subplot(3,2,1);
semilogy(t, abs((x-state(:, 1))));
yy = ylabel('$$X_{err}$$', 'interpreter', 'latex', 'rotation', 0);
set(yy, 'Units', 'Normalized', 'Position', [spacing, 0.5, 0]);
grid on; box on;

subplot(3,2,3);
semilogy(t, abs((y-state(:, 2))));
yy = ylabel('$$Y_{err}$$', 'interpreter', 'latex', 'rotation', 0);
set(yy, 'Units', 'Normalized', 'Position', [spacing, 0.5, 0]);
grid on; box on;

subplot(3,2,5);
semilogy(t, abs((z-state(:, 3))));
yy = ylabel('$$Z_{err}$$', 'interpreter', 'latex', 'rotation', 0);
set(yy, 'Units', 'Normalized', 'Position', [spacing, 0.5, 0]);
xlabel('$$Time\ [min]$$', 'interpreter', 'latex');
grid on; box on;

subplot(3,2,2);
semilogy(t, abs((vx-state(:, 4))));
yy = ylabel('$$Vx_{err}$$', 'interpreter', 'latex', 'rotation', 0);
set(yy, 'Units', 'Normalized', 'Position', [spacing, 0.5, 0]);
grid on; box on;

subplot(3,2,4);
semilogy(t, abs((vy-state(:, 5))));
yy = ylabel('$$Vy_{err}$$', 'interpreter', 'latex', 'rotation', 0);
set(yy, 'Units', 'Normalized', 'Position', [spacing, 0.5, 0]);
grid on; box on;

subplot(3,2,6);
semilogy(t, abs((vz-state(:, 6))));
yy = ylabel('$$Vz_{err}$$', 'interpreter', 'latex', 'rotation', 0);
set(yy, 'Units', 'Normalized', 'Position', [spacing, 0.5, 0]);
xlabel('$$Time\ [min]$$', 'interpreter', 'latex');
grid on; box on;

figure; hold on;
plot(t, energy-energySTK);
title('Specific Energy');
xlabel('$$Time\ [min]$$', 'interpreter', 'latex');
yy = ylabel('$$[MJ/kg]$$', 'interpreter', 'latex', 'rotation', 0);
set(yy, 'Units', 'Normalized', 'Position', [-0.11, 0.5, 0]);
grid on; box on;

% figure; hold on;
% plot3(x, y, z);
% plot3(state(:,1), state(:, 2), state(:,3));
% grid on; box on; axis equal;

%% Simulate Relative Dynamics
clear; close; clc;

% constants
p.mu = 398600.4418; % km^3/s^2
p.J2 = 1.08262668e-3;
p.Re = 6.3781e3; % km
p.kJ2 = (3*(p.J2)*(p.mu)*(p.Re)^2)/2; % km^5/s^2

% initial conditions
load('ImageSat_Inertial_Position_Velocity.mat');
x = data(:, 1);
y = data(:, 2);
z = data(:, 3);
vx = data(:, 4);
vy = data(:, 5);
vz = data(:, 6);
rvec = [x(1); y(1); z(1)]; % km
vvec = [vx(1); vy(1); vz(1)]; % km/s
[a0, e0, Omega0, I0, w0, nu0] = pos_vel2orb_ele(rvec, vvec, p.mu);
% [a0 ,e0 ,E0 ,I0 , w0 , Omega0 ,P , tau ,A , B] = vec2orbElem (rvec, vvec, p.mu);
% nu0 = 2*atan2(sqrt((1+e0)/(1-e0))*tan(E0/2), 1);
% z = [r; rd; h; th; I; Omg; x1; y1; z1; x1d; y1d; z1d; ...]
r0 = norm(rvec);
rd0 = vvec'*(rvec./norm(rvec));
h0 = norm(cross(rvec, vvec));
th0 = w0 + nu0;
% 100; 100; 100; .1; .1; -.1
z0 = [r0; rd0; h0; th0; I0; Omega0;...
    -1; 0; 0; .01; 0; 0];

% integration
opts.RelTol = 5e-14;
opts.AbsTol = 5e-14;
[t, zarray] = ode45(@relDyns, time.*60, z0, opts, p);
t = t./60;
r = zarray(:, 1);
rd = zarray(:, 2);
h = zarray(:, 3);
th = zarray(:, 4);
I = zarray(:, 5);
Omg = zarray(:, 6);

% convert LVLH to ECI
numSat = (length(z0) - 6)/6;
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
        
        % position
        state(i, (6*j-5):(6*j-3)) = (I_R_LVLH*[(xj + r(i)); yj; zj])';

        % velocity
        vjI = I_R_LVLH*[(xdj + rd(i) - yj*wz);...
            (ydj + (r(i) + xj)*wz - zj*wx);...
            (zdj + yj*wx)];
        state(i, (6*j-2):(6*j)) = vjI';
    end
end

% plot in ECI frame
figure; hold on;
plot3(center(:, 1), center(:, 2), center(:, 3));
plot3(state(:,1), state(:,2), state(:,3));
grid on; box on; axis equal;

% plot in LVLH frame
figure;
plot3(zarray(:, 7), zarray(:, 8), zarray(:, 9));
grid on; box on; axis equal;
title('LVLH Frame');

% Time series
figure;
subplot(3,1,1);
plot(t, zarray(:, 7));
ylabel('xj [km]');

subplot(3,1,2);
plot(t, zarray(:, 8));
ylabel('yj [km]');

subplot(3,1,3);
plot(t, zarray(:, 9));
ylabel('zj [km]');
xlabel('$$Time\ [min]$$', 'interpreter', 'latex');


%% Test1 to test relative dynamics
clear; close; clc;

% constants
p.mu = 398600.4418; % km^3/s^2
p.J2 = 1.08262668e-3;
p.Re = 6.3781e3; % km
p.kJ2 = (3*(p.J2)*(p.mu)*(p.Re)^2)/2; % km^5/s^2

% Test 1 ephemeris
load('Test1.mat');
x = data(:, 1);
y = data(:, 2);
z = data(:, 3);
vx = data(:, 4);
vy = data(:, 5);
vz = data(:, 6);

% Full Nonlinear Relative Dynamics
opts.RelTol = 5e-14;
opts.AbsTol = 5e-14;
[t, zarray] = ode45(@relDyns, time.*60, z0, opts, p);
t = t./3600;
r = zarray(:, 1);
rd = zarray(:, 2);
h = zarray(:, 3);
th = zarray(:, 4);
I = zarray(:, 5);
Omg = zarray(:, 6);

% Linear J2 Relative Dynamics
[tr, zarrayr] = ode45(@relDynsSmallr, time.*60, z0, opts, p);
tr = tr./3600;
rr = zarrayr(:, 1);
rdr = zarrayr(:, 2);
hr = zarrayr(:, 3);
thr = zarrayr(:, 4);
Ir = zarrayr(:, 5);
Omgr = zarrayr(:, 6);

% convert LVLH to ECI (Nonlinear)
numSat = (length(z0) - 6)/6;
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

% convert LVLH to ECI (Linear)
energyr = zeros(length(rr), numSat);
stater = zeros(length(rr), 6*numSat);
centerr = zeros(length(rr), 6);
for i = 1:length(rr)
    R_Omgr = [cos(Omgr(i)), -sin(Omgr(i)), 0;...
        sin(Omgr(i)), cos(Omgr(i)), 0;...
        0, 0, 1];
    R_Ir = [1, 0, 0;...
        0, cos(Ir(i)), -sin(Ir(i));...
        0, sin(Ir(i)), cos(Ir(i))];
    R_thr = [cos(thr(i)), -sin(thr(i)), 0;...
        sin(thr(i)), cos(thr(i)), 0;...
        0, 0, 1];
    I_R_LVLHr = R_Omgr*R_Ir*R_thr;
    
    % angular velocities
    wzr = hr(i)/rr(i)^2;
    wxr = -p.kJ2*sin(2*Ir(i))*sin(thr(i))/(hr(i)*rr(i)^3);
    
    % center of LVLH frame
    centerr(i, 1:3) = (I_R_LVLHr*[rr(i); 0; 0])';
    vIr = I_R_LVLHr*[rdr(i); hr(i)/rr(i); 0];
    centerr(i, 4:6) = vIr';
    
    for j = 1:numSat
        % relative positions and velocities
        xjr = zarrayr(i, 6*j+1); yjr = zarrayr(i, 6*j+2); zjr = zarrayr(i, 6*j+3);
        xdjr = zarrayr(i, 6*j+4); ydjr = zarrayr(i, 6*j+5); zdjr = zarrayr(i, 6*j+6);
        rjr = norm([(rr(i) + xjr), yjr, zjr]);
        rjZr = (rr(i) + xjr)*sin(Ir(i))*sin(thr(i)) + yjr*sin(Ir(i))*cos(thr(i)) + zjr*cos(Ir(i));
        
        % position
        stater(i, (6*j-5):(6*j-3)) = (I_R_LVLHr*[(xjr + rr(i)); yjr; zjr])';

        % velocity
        vjIr = I_R_LVLHr*[(xdjr + rdr(i) - yjr*wzr);...
            (ydjr + (rr(i) + xjr)*wzr - zjr*wxr);...
            (zdjr + yjr*wxr)];
        stater(i, (6*j-2):(6*j)) = vjIr';
        
        % total energy
        energyr(i, j) = 0.5*(vjIr')*vjIr + (-p.mu/rjr - p.kJ2/(3*rjr^3) + p.kJ2*(rjZr^2)/(rjr^5));
    end
end

% plot in ECI frame (Nonlinear)
spacing = -0.14;
figure;
subplot(3,2,1);
semilogy(t, abs((x-state(:, 1))./x)); 
yy = ylabel('$$X_{err}\ [\%]$$', 'interpreter', 'latex', 'rotation', 0);
set(yy, 'Units', 'Normalized', 'Position', [spacing, 0.5, 0]);
grid on; box on;

subplot(3,2,3); 
semilogy(t, abs((y-state(:, 2))./y)); 
yy = ylabel('$$Y_{err}\ [\%]$$', 'interpreter', 'latex', 'rotation', 0);
set(yy, 'Units', 'Normalized', 'Position', [spacing, 0.5, 0]);
grid on; box on;

subplot(3,2,5);
semilogy(t, abs((z-state(:, 3))./z));
yy = ylabel('$$Z_{err}\ [\%]$$', 'interpreter', 'latex', 'rotation', 0);
set(yy, 'Units', 'Normalized', 'Position', [spacing, 0.5, 0]);
xlabel('$$Time\ [hrs]$$', 'interpreter', 'latex');
grid on; box on;

subplot(3,2,2);
semilogy(t, abs((vx-state(:, 4))./vx));
yy = ylabel('$$Vx_{err}\ [\%]$$', 'interpreter', 'latex', 'rotation', 0);
set(yy, 'Units', 'Normalized', 'Position', [spacing, 0.5, 0]);
grid on; box on;

subplot(3,2,4);
semilogy(t, abs((vy-state(:, 5))./vy));
yy = ylabel('$$Vy_{err}\ [\%]$$', 'interpreter', 'latex', 'rotation', 0);
set(yy, 'Units', 'Normalized', 'Position', [spacing, 0.5, 0]);
grid on; box on;

subplot(3,2,6);
semilogy(t, abs((vz-state(:, 6))./vz));
yy = ylabel('$$Vz_{err}\ [\%]$$', 'interpreter', 'latex', 'rotation', 0);
set(yy, 'Units', 'Normalized', 'Position', [spacing, 0.5, 0]);
xlabel('$$Time\ [hrs]$$', 'interpreter', 'latex');
grid on; box on;

% plot energy
figure; hold on;
plot(t, energy)
% plot(tr, energyr)
grid on; box on;
title('$$Total\ Specific\ Energy\ [MJ/kg]$$', 'interpreter', 'latex');
xlabel('$$Time\ [hrs]$$', 'interpreter', 'latex');

% plot orbits
[X,Y,Z] = sphere(100);
figure; hold on;
plot3(x, y, z);
plot3(state(:, 1), state(:, 2), state(:, 3));
%s = surf(X.*p.Re,Y.*p.Re,Z.*p.Re);
%s.FaceColor = 'b';
grid on; box on; axis equal;
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
legend('STK', 'Nonlinear J2 Model');

%% Compare Nonlinear and Linear Model (close approach vs. high orbit)
clear; close; clc;
str = 'closed'; % choose 'h' for high orbit or 'c' for close approach 

% constants
p.mu = 398600.4418; % km^3/s^2
p.J2 = 1.08262668e-3;
p.Re = 6.3781e3; % km
p.kJ2 = (3*(p.J2)*(p.mu)*(p.Re)^2)/2; % km^5/s^2

if strcmp(str, 'c')
    % close approach
    load('z0.mat');
    z0(7:end) = [0.001; 0; 0; 0; 0; 0];
    % 20500
    tspan = linspace(0, 20500, 2880);
elseif strcmp(str, 'h')
    % high orbit
    load('z02.mat');
    z0 = [z0; 0.001; 0; 0; 0; 0; 0];
    % 150000
    tspan = linspace(0, 15000, 28800);
elseif strcmp(str, 'closed')
    % high orbit closed orbit
    load('OptSol14.mat');
    z0 = sol(1:12);
    tspan = linspace(0, 85000, 28800);
end


% Full Nonlinear Relative Dynamics
opts.RelTol = 5e-14;
opts.AbsTol = 5e-14;
[t, zarray] = ode45(@relDyns, tspan, z0, opts, p);
t = t./3600;
r = zarray(:, 1);
rd = zarray(:, 2);
h = zarray(:, 3);
th = zarray(:, 4);
I = zarray(:, 5);
Omg = zarray(:, 6);

% Linear J2 Relative Dynamics
[tr, zarrayr] = ode45(@relDynsSmallr, tspan, z0, opts, p);
tr = tr./3600;
rr = zarrayr(:, 1);
rdr = zarrayr(:, 2);
hr = zarrayr(:, 3);
thr = zarrayr(:, 4);
Ir = zarrayr(:, 5);
Omgr = zarrayr(:, 6);

% convert LVLH to ECI (Nonlinear)
numSat = (length(z0) - 6)/6;
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

% convert LVLH to ECI (Linear)
energyr = zeros(length(rr), numSat);
stater = zeros(length(rr), 6*numSat);
centerr = zeros(length(rr), 6);
for i = 1:length(rr)
    R_Omgr = [cos(Omgr(i)), -sin(Omgr(i)), 0;...
        sin(Omgr(i)), cos(Omgr(i)), 0;...
        0, 0, 1];
    R_Ir = [1, 0, 0;...
        0, cos(Ir(i)), -sin(Ir(i));...
        0, sin(Ir(i)), cos(Ir(i))];
    R_thr = [cos(thr(i)), -sin(thr(i)), 0;...
        sin(thr(i)), cos(thr(i)), 0;...
        0, 0, 1];
    I_R_LVLHr = R_Omgr*R_Ir*R_thr;
    
    % angular velocities
    wzr = hr(i)/rr(i)^2;
    wxr = -p.kJ2*sin(2*Ir(i))*sin(thr(i))/(hr(i)*rr(i)^3);
    
    % center of LVLH frame
    centerr(i, 1:3) = (I_R_LVLHr*[rr(i); 0; 0])';
    vIr = I_R_LVLHr*[rdr(i); hr(i)/rr(i); 0];
    centerr(i, 4:6) = vIr';
    
    for j = 1:numSat
        % relative positions and velocities
        xjr = zarrayr(i, 6*j+1); yjr = zarrayr(i, 6*j+2); zjr = zarrayr(i, 6*j+3);
        xdjr = zarrayr(i, 6*j+4); ydjr = zarrayr(i, 6*j+5); zdjr = zarrayr(i, 6*j+6);
        rjr = norm([(rr(i) + xjr), yjr, zjr]);
        rjZr = (rr(i) + xjr)*sin(Ir(i))*sin(thr(i)) + yjr*sin(Ir(i))*cos(thr(i)) + zjr*cos(Ir(i));
        
        % position
        stater(i, (6*j-5):(6*j-3)) = (I_R_LVLHr*[(xjr + rr(i)); yjr; zjr])';

        % velocity
        vjIr = I_R_LVLHr*[(xdjr + rdr(i) - yjr*wzr);...
            (ydjr + (rr(i) + xjr)*wzr - zjr*wxr);...
            (zdjr + yjr*wxr)];
        stater(i, (6*j-2):(6*j)) = vjIr';
        
        % total energy
        energyr(i, j) = 0.5*(vjIr')*vjIr + (-p.mu/rjr - p.kJ2/(3*rjr^3) + p.kJ2*(rjZr^2)/(rjr^5));
    end
end

% plot comparisons between Nonlinear and Linear
spacing = -0.12;
figure; 
subplot(4,4,[1,2]);
semilogy(t, abs(state(:, 1) - stater(:, 1)));
grid on; box on;
yy = ylabel('$$|X_{err}|\ [km]$$', 'interpreter', 'latex', 'rotation', 0);
set(yy, 'Units', 'Normalized', 'Position', [spacing, 0.5, 0]);
title('$$Nonlinear\ vs.\ Linear\ J2\ Relative\ Dynamics$$', 'interpreter', 'latex');

subplot(4,4,[5,6]);
semilogy(t, abs(state(:, 2) - stater(:, 2)));
grid on; box on;
yy = ylabel('$$|Y_{err}|\ [km]$$', 'interpreter', 'latex', 'rotation', 0);
set(yy, 'Units', 'Normalized', 'Position', [spacing, 0.5, 0]);

subplot(4,4,[9,10]);
semilogy(t, abs(state(:, 3) - stater(:, 3)));
grid on; box on;
yy = ylabel('$$|Z_{err}|\ [km]$$', 'interpreter', 'latex', 'rotation', 0);
set(yy, 'Units', 'Normalized', 'Position', [spacing, 0.5, 0]);

subplot(4,4,[13,14]);
semilogy(t, vecnorm(zarrayr(:, 7:9)'));
grid on; box on;
xlabel('$$Time\ [hrs]$$', 'interpreter', 'latex');
yy = ylabel('$$\|\mathbf{\rho}\|\ [km]$$', 'interpreter', 'latex', 'rotation', 0);
set(yy, 'Units', 'Normalized', 'Position', [spacing, 0.5, 0]);

% plot Linear
subplot(4,4,[3,4,7,8,11,12,15,16]); hold on;
plot3(stater(:,1), stater(:, 2), stater(:,3));
plot3(state(:,1), state(:, 2), state(:,3));
grid on; box on; axis equal;
title('Nonlinear vs. Linear J2 Relative Dynamics');
xlabel('$$X\ [km]$$','interpreter','latex');
ylabel('$$Y\ [km]$$','interpreter','latex');
zlabel('$$Z\ [km]$$','interpreter','latex');
legend('Linear', 'Nonlinear');

% plot energies
figure; hold on;
plot(t, energyr);
plot(t, energy);
grid on; box on; 
legend('Linear', 'Nonlinear');
xlabel('$$Time\ [hrs]$$', 'interpreter', 'latex');
title('Total Specific Energy [MJ/kg]');

%% Comparison with CLohessy Wiltshire
clear; close; clc;

% constants
p.mu = 398600.4418; % km^3/s^2
p.J2 = 1.08262668e-3;
p.Re = 6.3781e3; % km
p.kJ2 = (3*(p.J2)*(p.mu)*(p.Re)^2)/2; % km^5/s^2

% % high circular orbit
% load('z03.mat');
% zrel0 = [1000; 0; 0; 0; 0; 0];
% z0 = [z0; zrel0];
% tspan = linspace(0, 100000, 2880);
% numSat = (length(z0) - 6)/6;
% rvec = [25000.1; 1.28e-18; 5.969e-19]; % km
% vvec = [-2.25613e-22; 3.61888; 1.68751]; % km/s

load('OptSol14.mat');
z0 = sol(1:12);
zrel0 = z0(7:end);
numSat = (length(z0) - 6)/6;
tspan = linspace(0, 100000, 2880);

% Full Nonlinear Relative Dynamics
opts.RelTol = 5e-14;
opts.AbsTol = 5e-14;
[t, zarray] = ode45(@relDyns, tspan, z0, opts, p);
t = t./3600;
r = zarray(:, 1);
rd = zarray(:, 2);
h = zarray(:, 3);
th = zarray(:, 4);
I = zarray(:, 5);
Omg = zarray(:, 6);
[energy, state, center] = LVLH_to_ECI(zarray, numSat, p);

% Linear J2 Relative Dynamics
[tr, zarrayr] = ode45(@relDynsSmallr, tspan, z0, opts, p);
tr = tr./3600;
rr = zarrayr(:, 1);
rdr = zarrayr(:, 2);
hr = zarrayr(:, 3);
thr = zarrayr(:, 4);
Ir = zarrayr(:, 5);
Omgr = zarrayr(:, 6);
[energyr, stater, centerr] = LVLH_to_ECI(zarrayr, numSat, p);

% Clohessy Wiltshire
Ic = deg2rad(25.0001);
c.mu = p.mu;
c.a = 20000; % km semi-major axis
c.kJ2 = p.kJ2;
c.I = Ic;
[tc, zarrayc] = ode45(@ClohessyWiltshire, tspan, zrel0, opts, c);
tc = tc./3600;
xc = zarrayc(:, 1);
yc = zarrayc(:, 2);
zc = zarrayc(:, 3);
xdc = zarrayc(:, 4);
ydc = zarrayc(:, 5);
zdc = zarrayc(:, 6);
nc = sqrt(c.mu/c.a^3);
Omgc = 0;
R_Omg = [cos(Omgc), -sin(Omgc), 0;...
    sin(Omgc), cos(Omgc), 0;...
    0, 0, 1];
R_I = [1, 0, 0;...
    0, cos(Ic), -sin(Ic);...
    0, sin(Ic), cos(Ic)];
statec = zeros(length(tc), 6);
energyc = zeros(length(tc), 1);
for i = 1:length(tc)
    thc = nc*tc(i);
    R_th = [cos(thc), -sin(thc), 0;...
        sin(thc), cos(thc), 0;...
        0, 0, 1];
    I_R_LVLH = R_Omg*R_I*R_th;
    statec(i, 1:3) = (I_R_LVLH*[(c.a+xc(i)); yc(i); zc(i)])';
    statec(i, 4:6) = (I_R_LVLH*[(xdc(i) - yc(i)*nc); ((c.a + xc(i))*nc + xdc(i)); zdc(i)])';
end

% plot
figure; 
subplot(3,2,1); hold on;
plot(t, state(:, 1), 'k-');
plot(tr, stater(:, 1), 'b-');
plot(tc, statec(:,1), 'g-');
grid on; box on; 
subplot(3,2,3); hold on;
plot(t, state(:, 2), 'k-');
plot(tr, stater(:, 2), 'b-');
plot(tc, statec(:,2), 'g-');
grid on; box on; 
subplot(3,2,5); hold on;
plot(t, state(:, 3), 'k-');
plot(tr, stater(:, 3), 'b-');
plot(tc, statec(:,3), 'g-');
grid on; box on; 

subplot(3,2,2); hold on;
plot(t, state(:, 4), 'k-');
plot(tr, stater(:, 4), 'b-');
plot(tc, statec(:,4), 'g-');
grid on; box on; 
subplot(3,2,4); hold on;
plot(t, state(:, 5), 'k-');
plot(tr, stater(:, 5), 'b-');
plot(tc, statec(:,5), 'g-');
grid on; box on; 
subplot(3,2,6); hold on;
plot(t, state(:, 6), 'k-');
plot(tr, stater(:, 6), 'b-');
plot(tc, statec(:,6), 'g-');
grid on; box on; 

% plot 3D
figure; hold on;
plot3(state(:, 1), state(:, 2), state(:, 3));
plot3(stater(:, 1), stater(:, 2), stater(:, 3));
plot3(statec(:, 1), statec(:, 2), statec(:, 3));
grid on; box on; axis equal;


%% Frozen Orbit (OptSol1-2) for Leader Satellite
% (OptSol7) - +/- y0 satellite will stay there since orbit precess by same
% mount
% (OptSol9) double z
% (OptSol11) kissing orbits
% (OptSol12) X rings
% (OptSol14) 3 rings
% (OptSol17) heart shape tall 
% (OptSol19) Test1
clear; close all; clc;

% constants
p.mu = 398600.4418; % km^3/s^2
p.J2 = 1.08262668e-3;
p.Re = 6.3781e3; % km
p.kJ2 = (3*(p.J2)*(p.mu)*(p.Re)^2)/2; % km^5/s^2
p.kJ2 = 0;

% I = 63.4
% load('z04.mat');
% zrel0 = [0; 0.1; 0.1; -.01; 0; 0.01];
% z0 = [z0; zrel0];
% tspan = linspace(0, 1000000, 2880);
% numSat = (length(z0) - 6)/6;
% rvec = [20000; 0; 0]; % km
% vvec = [0; 1.996498038566532; 3.992996077133058]; % km/s
load('OptSol14.mat');
z0 = sol(1:(end-1));
tspan = linspace(0, 400000, 1e5);
numSat = length(z0)/6 - 1;

% Full Nonlinear Relative Dynamics
opts.RelTol = 5e-14;
opts.AbsTol = 5e-14;
[time, zarray] = ode45(@relDyns, tspan, z0, opts, p);
t = time./3600;
r = zarray(:, 1);
rd = zarray(:, 2);
h = zarray(:, 3);
th = zarray(:, 4);
I = zarray(:, 5);
Omg = zarray(:, 6);
[energy, state, center] = LVLH_to_ECI(zarray, numSat, p);

% animation
animation = false;
if animation
    fig = figure; hold on;
    cc = plot3(0, 0, 0, 'k*');
    col = 'kbg'; k = 1;
    for i = 1:6:(size(state, 2))
        plot3(zarray(:, i+6), zarray(:, i+7), zarray(:, i+8), [col(k), '-']);
        k = k + 1;
    end
    legend('Leader','Follower1', 'Follower2', 'Follower3');
    title('LVLH Frame');
    xlabel('X [km]');
    ylabel('Y [km]');
    zlabel('Z [km]');
    tic;
    current = toc;
    while current < time(end)
        % interpolate
        current = toc*1000;
        sat1x = interp1(time, zarray(:, 7), current);
        sat1y = interp1(time, zarray(:, 8), current);
        sat1z = interp1(time, zarray(:, 9), current);
        sat2x = interp1(time, zarray(:, 13), current);
        sat2y = interp1(time, zarray(:, 14), current);
        sat2z = interp1(time, zarray(:, 15), current);
        sat3x = interp1(time, zarray(:, 19), current);
        sat3y = interp1(time, zarray(:, 20), current);
        sat3z = interp1(time, zarray(:, 21), current);
        s1 = plot3(sat1x, sat1y, sat1z, 'k.','MarkerSize', 30);
        s2 = plot3(sat2x, sat2y, sat2z, 'b.','MarkerSize', 30);
        s3 = plot3(sat3x, sat3y, sat3z, 'g.','MarkerSize', 30);
        axis equal; box on; grid on;
        view(3);
        drawnow;
        if ~ishandle(fig)
            close all;
            break;
        end
        delete(s1); delete(s2); delete(s3);
    end
else
    % plot
    figure; 
    subplot(1,2,1); hold on;
    plot3(center(:, 1), center(:, 2), center(:, 3));
    for i = 1:6:(size(state, 2))
        plot3(state(:, i), state(:, i+1), state(:, i+2));
    end
    grid on; box on; axis equal;
    title('ECI frame');
    xlabel('X [km]');
    ylabel('Y [km]');
    zlabel('Z [km]');
    legend('Leader', 'Follower1', 'Follower2', 'Follower3');

    % LVLH
    subplot(1,2,2); hold on;
    cc = plot3(0, 0, 0, 'k*');
    for i = 1:6:(size(state, 2))
        plot3(zarray(:, i+6), zarray(:, i+7), zarray(:, i+8));
    end
    grid on; box on; axis equal;
    title('LVLH');
    xlabel('X [km]');
    ylabel('Y [km]');
    zlabel('Z [km]');
    legend('Center', 'Follower1', 'Follower2', 'Follower3');
    view(3);
end

%% Initial conditions for natural motions (high orbit)
clear; close all; clc;

% load('Test1.mat');

mu = 398600.4418; % km^3/s^2
% z0(1) = 7e3;
% z0(7) = 0;
% z0(11) = -0.05;
% z0(12) = 0;
% a = z0(1); % from Test 1
% x0 = [z0; T];
a = 20000; % km
T = 2*pi*sqrt(a^3 / mu);
x0 = [a;0;89286.1066235951;0;1.10714871779409;0;...
    0; 10; 0; -0.1; 0; 0;...
    0; -10; 0; 0.1; 0; 0;...
    T];
% x0 = [a;0;89286.10662359512;0;0.43633231300;0;...
%     0; 0; 0.5; 0; -0.1; -0.1;...
%     T];
options = optimoptions('fminunc', 'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e6,...
    'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, 'Display', 'iter');
sol = fminunc(@errNorm, x0, options);
% options = optimoptions('lsqnonlin', 'Display', 'iter', 'FunctionTolerance', 1e-12,...
%     'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e6, 'OptimalityTolerance', 1e-12, 'Algorithm', 'levenberg-marquardt');
% [sol, resnorm] = lsqnonlin(@errNorm_lsqnonlin, x0, [], [], options);

save('OptSol20.mat', 'sol');

