%% mass_matrix_3d 
% forms the mass matrix with inputs in 3D
% dp = distance plus, dm = distance minus, N  = number of links, masses =
% masses of links, qlist = quaternions at current timestep
% written by Selen Serdar, last updated 4/4/2026

function M = mass_matrix_3d(dp, dm, N, masses, qlist, J)

m = zeros(6*N, 6*N); % preallocate mass matrix size
tracksum_row = zeros(3,6*N); % tracking summations in terms - does not reset each iteration

% convert the quats back into DCM for use 
for k = 1:N
    q = qlist(4*(k-1)+1 : 4*k)'; 
    q = q / norm(q); % normalize 
    Rlist(:,:,k) = quat2dcm(q);
end

%% First half of the mass matrix ==========================================
c = 1; % count link number
for idx = 1:3*N % compute eqn 29 
    
    % get the current values for J, R, mass
    Jn = J(:,:,c);
    R = Rlist(:,:,c);
    mass = masses(c); 
    
    %% account for the summations of previous iterations
    if idx == 1 || c == 1
        tracksum_row(1:3, idx + N) = [0;0;0];
    else
        prevm = masses(c-1); % previous link mass [kg]
        tracksum_row(1:3, idx + N) = [prevm; prevm; prevm] ; 
    end

    row = zeros(3,6*N); % initialize row 

    %% put in the stepwise-varied terms (coefficients of ddotX and dotw)
    % for the 0th link:
    if idx == 1 
        row(1:3, idx:idx+2) = Jn(1:3,1:3); % angular vel term 
        row(1:3,idx+N:idx+2+N) = -skew(dp) * R * mass; % linear accel term 
        extra = 1; 
        
    % compute for the middle links: 
    elseif idx > 1 && idx < N
        row(1:3, idx:idx+2) = Jn; % angular vel term 
        row(1:3, idx+N:idx+2+N) = -skew(dp) * R * mass; % linear acc 
        extra = skew(dm - dp) * R;

    % for the final link: 
    elseif idx == N
        row(1:3, idx:idx+2) = Jn; % angular velocity
        row(1:3, idx+N:idx+2+N) = 0; % linear acc
        extra = -skew(dm) * R;
    end
    
    row = row + extra * tracksum_row ; % incorporate factor multiplying the sums
    
    % put the row into m 
    row_ind = 3*(c-1)+1 : 3*c;
    m(row_ind,:) = row;

    % decide which link it is calculating
    if c == 20
        c = 1; % reset c 
    else
        c = c + 1;
    end
end

clear row idx c R tracksum_row

%% Second half of the mass matirx =========================================
% use eqn 26 - same for each link
tracksum_row = zeros(3, 6*N); % reinitialize sum tracking 
c = 1; % counting variable 1 - 20 

for idx = 1:6*N-3
    
    R = Rlist(:,:,c); % get the current R 
    if c < 20
        nextR = Rlist(:,:,c+1); % get the next R 
    else
        nextR = Rlist(:,:,1);
    end
    row = zeros(3, 6*N); % reinitialize row


    % angular velocity multipliers:
    if idx > 1
        tracksum_row(1:3, idx:idx+2) = R' * skew(dp);
        tracksum_row(1:3, idx+1:idx+3) = nextR'* skew(dm) ;
    end

    % combine row with the reappearing terms
    row(1:3,idx:idx+2) = row(1:3,idx:idx+2) + tracksum_row(1:3, idx:idx+2);
    row(1:3,idx+1:idx+3) = row(1:3, idx+1:idx+3) + tracksum_row(1:3, idx+1:idx+3);

    % put the row into m 
    row_ind = 3*N + 3*(c-1)+1 : 3*N + 3*c;
    m(row_ind,:) = row;

    % reset c 
    if c == 20
        c = 1; 
    else
        c = c + 1;
    end
end

% put in the 1s for linear acceleration
m(3*N:6*N,3*N+1) = -1;
for k = 3*N+1 : 6*N
    m(k, k) = 1;     
end

M = sparse(m); % return a sparse mass matrix

end