% coupled orbit attitude MAIN - 3D
%
% propagate coupled orbit attitude simulation in ode45.
% written by Selen Serdar, last updated 4/3/2026

clc;clear;close all;

N = 20; % number of links

% INITIAL CONDITIONS ======================================================

omega_tot = 2 * pi / (22.5 * 60);  % How fast is the tether spinning? [rad/s]
% make sure the R matrices and omegas work
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
omegas = zeros(3*N,1); % angular velocity of each link [x, y, z] [rad/s]
R = zeros(3*N,1); % rotation matrix of each link (link to world)
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

% create a cell array for J
J = {};
% form each link's inertia tensor
for k = 1:N
    m = masses(k); 
    J(k) = [1/12 *m * (3*r^2 + link_length^2) 0 0;
        0 1/12 * m *(3*r^2 + link_length^2) 0;
        0 0 1/2*m*R^2]; 
end

% set up d vectors 
dplus = ones(1,N-1); % how far is the current link's CoM from the next link's joint? [m]
dminus = ones(1,N-1); % how far is current link's CoM from previous link's joint? [m]

% assemble the initial state vector (column)
sv0 = [omegas x R phi tau]';

% timespan [s]
tspan = 0:.01:10;

% PROPAGATE ===============================================================
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

[tout, yout] = ode45(@(t, sv, masses, J, N, dplus, dminus)coupled_orbit_1d, tspan, sv0, options);

omegas_final = yout(:, 1:3);
x_final = yout(:, 4:6);

% plot the results 
% figure(1)
% plot(tspan, omegas_final)
% hold on
% plot(tspan, x_final)
% legend("angular velocity", "position")