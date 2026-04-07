% mass matrix N=4

%% mass_matrix_3d 
% forms the mass matrix with inputs in 3D
% dp = distance plus, dm = distance minus, N  = number of links, masses =
% masses of links, R = rotation matrices at current timestep
% written by Selen Serdar, last updated 4/6/2026

function M = mass_matrix_3d_4(dp, dm, N, masses, R, J)

m = zeros(6*N, 6*N); % preallocate mass matrix size


for i = 1:N
    m(3*i-2: 3*i, 3*i-2:3*i) = J(:,:,i);

    A = -masses(i) * skew(dp) * R(:,:,i);
    m(3*i-2:3*i , 3*i-2+3*N:3*i+3*N) = A;

    for beta = 1:1:i-1
        B = -skew(dp - dm) * R(:,:,i) * masses(beta);

        m(3*i-2:3*i , 3*beta-2+3*N:3*beta+3*N) = B;
    end

    m(3*N+1:3*N+3, 3*i-2+3*N:3*i+3*N) = masses(i)*eye(3);
end

for i = 2:N
    m(3*i-2+3*N:3*i+3*N, 3*N+1:3*N+3) = -eye(3);
    m(3*i-2+3*N:3*i+3*N, 3*i-2+3*N:3*i+3*N) = eye(3);
    for beta = 1:i-1
        A = R(:,:,i-1)' * skew(dp);
        B = -R(:,:,i)' * skew(dm);

        m(3*i-2+3*N:3*i+3*N, 3*beta-2:3*beta) = m(3*i-2+3*N:3*i+3*N, 3*beta-2:3*beta) + A;
        m(3*i-2+3*N:3*i+3*N, 3*beta+1:3*beta+3) =m(3*i-2+3*N:3*i+3*N, 3*beta+1:3*beta+3) + B;
    end
end

M = sparse(m);

end