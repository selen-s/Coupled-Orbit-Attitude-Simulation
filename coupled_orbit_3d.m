%% coupled_orbit_3d
% compute the coupled orbit attitude problem in 3D. 
% sv_derivative: inputs are current time(t) in s, state vector (sv), link
%masses (masses) in kg, moment of inertia tensor (J) , and number of links
%(N)
% outputs 
% written by Selen Serdar, last updated 4/5/2026
function [sv_dot] = coupled_orbit_3d(t, sv, masses, J, N, dp, dm)

% categorize inputs from state vector. these will go into the coeff matrix
x = sv(1:3*N); % positions (m)
x_dot = sv(3*N+1:6*N); % velocities (m/s)
omega = sv(6*N+1:9*N); % omegas (rad/s)
R = sv(6*N+1:15*N); % rotation matrices

R = reshape(R, [3, 3, N]); % put R into useable shapes 

massM = mass_matrix_3d_4(dp, dm, N, masses, R, J); % get the nass natrix
coeffM = coeff_matrix_3d_4(dp, dm , N, R, J, omega); % get coeff matrix

% get the derivative of everything
omega = reshape(omega, [3,1,N]);
for i = 1:N
    % R dot:
    Rdot(:,:,i) = -skew(omega(:,i)) * R(:,:,i);
end

Rdot_lin = Rdot(:); % flatten R dot
wx_derivative = coeffM' / massM; % take derivative of omega and x terms
omega_dot = wx_derivative(1:3*N)'; % omega dot
xddot = wx_derivative(3*N+1:6*N)'; % x double dot

sv_dot = [x_dot; xddot; omega_dot; Rdot_lin]; % state derivative

end