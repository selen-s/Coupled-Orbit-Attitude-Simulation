% coupled orbit attitude MAIN
%
% propagate coupled orbit attitude simulation in ode45.
% written by Selen Serdar, last updated 3/29/2026
clc;clear;close all;
N = 20; % number of links
% initial conditions
omegas = zeros(1,N); % angular velocity of each link 
R = zeros(1,N); % rotation matrix of each link (link to world)
x = zeros(1,N); % position of each link
phi = zeros(1,N); % external forces
tau = zeros(1,N); % torques 

% to calculate the moment of inertia:
r = 0.044;  % What is the radius of the tether cable everywhere? [m]
L = 300000;  % How long is the tether in total? [m]
rho = 1440;  % What is the density of the tether material [kg/m^3]

mc = 50000;  % How massive is the counterweight? [kg]
mt = 10000;  % How massive is the catch mechanism? [kg]

masses = [mc ones(1,N-2) mt]; % vector of masses [kg]

J = 

%omega = 2 * pi / (22.5 * 60);  % How fast is the tether spinning? [rad/s]



% set up d vectors
dplus = ones(1,N);
dminus = ones(1,N);

% assemble the initial state vector
sv0 = [omegas x R phi tau]';

% timespan 
tspan = 0:.01:10;

% propagate
options = odeset('RelTol',1e-12,'AbsTol',1e-12);

[tout, yout] = ode45(@(t, sv, masses, J, N, dplus, dminus)coupled_orbit_1d, tspan, sv0, options);


figure(1)
plot(tspan, yout)