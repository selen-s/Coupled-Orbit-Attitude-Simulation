%% coeff_matrix_3d
%{
takes inputs of: dp (dist to next link CoM), dm (dist to prev link CoM), N
(link number), masses (link masses), R (rot matrix), J (moi tensor)
outputs: C, matrix of coefficients to be used in solving for state vector
author: Selen Serdar, last updated: 4/6/2026
%}

function C = coeff_matrix_3d_4(dp, dm, N, Rlist, Jlist, omegas)

c = zeros(6*N, 1);
omegas = reshape(omegas, [3, N]);
count = 1; 

for i = 1:3:3*N - 3
    
    % get all current/next link values
    R = Rlist(:,:,count);
    nextR = Rlist(:,:,count+1);
    J = Jlist(:,:,count);
    omega = omegas(:, count); % get current link omega vector
    next_omega = omegas(:, count+1); % next

    %   eqn 26 RHS
    c(i:i+2) = c(i:i+2) + R' * cross(omega, cross(omega, dp)) - nextR' * cross(next_omega, cross(next_omega, dm)); 

    if count == N % which link is it on
        count = 1;
    else
        count = count + 1; 

    end
end

for i = 3*N-3 : 3*N
    c(i:i+2) = 0; % eqn 28
end

for i = 3*N+1 : 6*N - 6 % eqn 29 - 31 - assume no phi, tau
    c(i:i+2) = cross(omega, J*omega);
end

C = sparse(c);
end