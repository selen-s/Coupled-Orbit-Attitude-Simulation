% coupled orbit attitude MAIN - 3D
%
% propagate coupled orbit attitude simulation in ode45.
% written by Selen Serdar, last updated 4/5/2026

clc;clear;close all;
tic

N = 20; % number of links

% INITIAL CONDITIONS ======================================================
omega_tot = 2 * pi / (22.5 * 60);  % How fast is the tether spinning? [rad/s]
% make sure the R matrices and omegas work
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
omegas = zeros(3*N, 1) * omega_tot; % angular velocity of each link [x, y, z] [rad/s]
x = zeros(3*N,1); % position of each link relative to world (x,y,z) [m]
phi = zeros(3*N,1); % external forces [neglect for now]
tau = zeros(3*N,1); % torques [neglect for now]

% to calculate the moment of inertia:
r = 0.044;  % What is the radius of the tether cable everywhere? [m]
L = 300000;  % How long is the tether in total? [m]
rho = 1440;  % What is the density of the tether material [kg/m^3]
mc = 50000;  % How massive is the counterweight? [kg]
mt = 10000;  % How massive is the catch mechanism? [kg]
masses = [mc ones(1,N-2) mt]; % vector of masses [kg]
link_length = L / N; % length of each link [m]

J = zeros(3,3,N);
% form each link's inertia tensor
for k = 1:N
    m = masses(k); 
    J(:,:,k) = [1/12 *m * (3*r^2 + link_length^2) 0 0;
        0 1/12 * m *(3*r^2 + link_length^2) 0;
        0 0 1/2*m*r^2]; 
end

% set up d vectors in link frame
dplus = [link_length/2;0;0]; % how far is the current link's CoM from the next link's joint? [m]
dminus = [link_length/2;0;0]; % how far is current link's CoM from previous link's joint? [m]
 
% quaternions 
qlist = zeros(N,4);
for k = 1:N
    qlist(k, :) = [1 0 0 0];
end
qvec = reshape(qlist, [4*N,1]);

% assemble the initial state vector (column)
sv0 = [omegas; x; qvec]; % concatenate phi and tau if needed

% timespan [s]
tspan = 0:.01:10;

% PROPAGATE ===============================================================
options = odeset('RelTol',1e-12,'AbsTol',1e-12); % set ode45 options

[tout, yout] = ode45(@(t, sv)coupled_orbit_3d(t, sv,masses, J, N, dplus, dminus), tspan, sv0, options);

omegas_final = yout(:, 1:3); % angular velocity output
x_final = yout(:, 4:6); % pos output

% sort x, y, z, components for plotting
% angular velocity [rad/s]
wx = omegas_final(:, 1);
wy = omegas_final(:, 2);
wz = omegas_final(:, 3);

% velocity [m/s]
xx = x_final(:, 1);
xy = x_final(:, 2);
xz = x_final(:, 3);

% plot the results 
figure(1)
sgtitle("COA Simulation Results")

subplot(3,2, 1)
plot(tspan, wx)
title("Angular velocity (x) [rad/s]")
xlabel("time (s)")
ylabel("angular velocity (rad/s)")

subplot(3,2,2)
plot(tspan, xx)
title("Velocity (x) [m/s]")
xlabel("time (s)")
ylabel("velocity (m/s)")

subplot(3,2,3)
plot(tspan, wy)
title("Angular velocity (y) [rad/s]")
xlabel("time (s)")
ylabel("angular velocity (rad/s)")


subplot(3,2,4)
plot(tspan, xy)
title("Velocity (y) [m/s]")
xlabel("time (s)")
ylabel("velocity (m/s)")

subplot(3,2,5)
plot(tspan, wz)
title("Angular velocity (z) [rad/s]")
xlabel("time (s)")
ylabel("angular velocity (rad/s)")


subplot(3,2,6)
plot(tspan, xz)
title("Velocity (z) [m/s]")
xlabel("time (s)")
ylabel("velocity (m/s)")

toc