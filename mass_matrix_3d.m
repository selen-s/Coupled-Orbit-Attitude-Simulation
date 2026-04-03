%% mass_matrix_3d 
%n forms the mass matrix with inputs in 3D
% dp = distance plus, dm = distance minus, N  = number of links, masses =
% masses of links, Rlist = rotation matrices at current timestep
% written by Selen Serdar, last updated 4/3/2026

function M = mass_matrix_3d(dp, dm, N, masses, Rlist, J)

m = zeros(6*N,6*N); % preallocate matrix
tracksum_row = zeros(3,6*N); % tracking summations - does not reset each iteration

%% First half of the mass matrix ==========================================
for idx = 1:N % compute eqn 29 
 
    Jmat = J{idx}; 
    Jn = Jmat(1:3,1:3);
    R = Rlist(idx);
    %% account for the summations of previous iterations
    if idx == 1
        tracksum_row(1:3, idx + N) = [0;0;0];
    else
        tracksum_row(1:3, idx - 1 + N) = [masses(idx - 1); masses(idx-1); masses(idx-1)] ; 
    end

    row = zeros(3,6*N);
    %% put in the stepwise-varied terms (coefficients of ddotX and dotw)
    % for the 0th link:
    if idx == 1 
        row(1:3, idx:idx+2) = Jn(1:3,1:3); % angular accel term 
        row(1:3,idx+N:idx+2+N) = -skew(dp) * R * masses(idx); % linear accel term 
        extra = 1; 
    % compute for the middle links: 
    elseif idx > 1 && idx < N
        row(1:3, idx:idx+2) = Jn;
        row(1:3, idx+N:idx+2+N) = -skew(dp) * R * masses(idx); 
        extra = -skew(dp - dm) * R;
     % for the final link: 
    elseif idx == N
        row(1:3, idx:idx+2) = Jn; % angular accb
        row(1:3, idx+N:idx+2+N) = -skew(dm) * R; % linear acc
        extra = -skew(dm) * R;
    end
    
    row = row + extra * tracksum_row ; % incorporate factor multiplying the sums

    % append the new row to m:
    m = [m ; row];
end
 clear row
%% Second half of the mass matirx =========================================
% use eqn 26 - same for each link
for idx = 1:6*N-1
    R1 = Rlist(idx+1);
    row = zeros(3, 3*N); % re initialize row
    % angular velocity multipliers:
    row(1:3, idx:idx+2) = R' * skew(dp);
    row(1:3, idx+1:idx+3) = -R1' * skew(dm); 
    
    % ddotx multipliers:
    row(1:3, 1) = 1;
    row(1:3, idx) = 1 ;

    if idx == 6*N-6
        row(1:3, idx+1:idx+3) = 1;
    end
    % append the new row to m: 
    m = [m; row];

end

M = sparse(m); % return a sparse mass matrix

end