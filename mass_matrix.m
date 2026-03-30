%% mass_matrix 
%n forms the mass matrix with inputs
% dp = distance plus, dm = distance minus, N  = number of links, J = MOI
% tensor, masses = masses of links
% written by Selen Serdar, last updated 3/29/2026
%
function M = mass_matrix(dp, dm, N, J, masses)

m = zeros(2*N); % preallocate matrix
tracksum_row = zeros(2*N); % tracking summations - does not reset each iteration

for idx = 1:N % compute eqn 29 
 
    Jn = J{idx}; % get MOI for current link

    R = get_R(N); % get rotation mtx for the link 

    %% account for the summations of previous iterations
    tracksum_row(idx) = masses(idx) ;
    row = tracksum_row; 

    %% put in the stepwise-varied terms (coefficients of ddotX and ddotw)
    % for the 0th link:
    if idx == 1
        row(idx) = Jn; % angular accel term 
        row(idx + N) = skew(dp(N)) * R * masses(idx); % linear accel term 

    % compute for the middle links: 
    elseif idx > 1 && idx < N
        row(idx) = ;
        row(idx + N) = ; 

     % for the final link: 
    elseif idx == N
        row(idx) = Jn; % angular acc
        row(idx + N) = -skew(dm(idx)) * R; % linear acc
    end

    % append the new row to m:
    m = [m ; row];
end

M = sparse(m); % return a sparse mass matrix

end