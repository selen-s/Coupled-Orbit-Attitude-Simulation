%% mass_matrix_3d 
%n forms the mass matrix with inputs in 3D
% dp = distance plus, dm = distance minus, N  = number of links, masses =
% masses of links, Rlist = rotation matrices at current timestep
% written by Selen Serdar, last updated 4/3/2026

function M = mass_matrix_3d(dp, dm, N, masses, Rlist)

m = zeros(2*N); % preallocate matrix
tracksum_row = zeros(2*N); % tracking summations - does not reset each iteration

for idx = 1:N % compute eqn 29 
 
    Jn = J(idx); 

    R = Rlist(idx);
    %% account for the summations of previous iterations
    if idx == 1
        tracksum_row(idx + N) = 0;
    else
        tracksum_row(idx - 1 + N) = masses(idx - 1) ; 
    end
    %% put in the stepwise-varied terms (coefficients of ddotX and ddotw)
    % for the 0th link:
    if idx == 1
        row(idx) = Jn; % angular accel term 
        row(idx + N) = -skew(dp(idx)) * R * masses(idx); % linear accel term 
        
        extra = 1; 
    % compute for the middle links: 
    elseif idx > 1 && idx < N
        row(idx) = Jn;
        row(idx + N) = -skew(dp(idx)) * R * masses(idx); 
        extra = -skew(dp(idx) - dm(idx)) * R;
     % for the final link: 
    elseif idx == N
        row(idx) = Jn; % angular accb
        row(idx + N) = -skew(dm(idx)) * R; % linear acc
        extra = -skew(dm(idx)) * R;
    end
    
    row = row + tracksum_row * extra; % incorporate factor multiplying the sums

    % append the new row to m:
    m = [m ; row];
end

M = sparse(m); % return a sparse mass matrix

end