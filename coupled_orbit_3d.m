%% coupled_orbit_3d
% compute the coupled orbit attitude problem in 3D. 
% sv_derivative: inputs are current time(t) in s, state vector (sv), link
%masses (masses) in kg, moment of inertia tensor (J) , and number of links
%(N)
% outputs 
% written by Selen Serdar, last updated 4/3/2026
function [sv_derivative] = coupled_orbit_3d(t, sv, masses, J, N, dp, dm)

% categorize inputs from state vector. these will go into the coeff matrix
omega = sv(1:3*N); % angular velocities (rad/s)
x = sv(3*N+1:6*N); % positions (m)
R = sv(6*N+1:9*N); % rotation matrix 
phi = sv(9*N+1:12*N); % external forces (N)
tau = sv(12*N+1:15*N); % torques (Nm)

% create the matrix of constants:
% constM = constant_matrix(); 
constM = sparse(eye(15*N));
massM = mass_matrix_3d(dp, dm, N, masses, R, J); 

% return the derivative of state vector
sv_derivative = constM * sv / massM ; 

end
