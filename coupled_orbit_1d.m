%% coupled_orbit_1d 
% compute the coupled orbit attitude problem. 
% sv_derivative: inputs are current time(t) in s, state vector (sv), link
%masses (masses) in kg, moment of inertia tensor (J) , and number of links
%(N)
% outputs 
% written by Selen Serdar, last updated 3/29/2026
function [sv_derivative] = coupled_orbit_1d(t, sv, masses, J, N)

% categorize inputs from state vector.
w = sv(1:N); % angular velocities (rad/s)
x = sv(N+1:2*N); % positions (m)
R = sv(2*N+1:3*N); % rotation matrix 
phi = sv(3*N+1:4*N); % external forces (N)
tau = sv(4*N+1:5*N); % torques (Nm)

% create the matrix of constants:
% constM = constant_matrix(); 
massM = mass_matrix(w, x, R, phi, tau, N); 

% return the derivative of state vector
sv_derivative = inv(massM) * constM; 
end