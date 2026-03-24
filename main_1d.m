% coupled orbit attitude MAIN
clc;clear;close all;
N = 100;
omegas = zeros(1,N);
R = zeros(N);
x = zeros(1,N);
masses = ones(1,N); 
phi = zeros(1,N);
tau = zeros(1,N);

% moi
J = ones(1,N);

% set up d vectors
dplus = ones(1,N);
dminus = ones(1,N);

% assemble state vector
sv0 = [omegas x R phi tau]';

% [tout, yout] = ode45(      );