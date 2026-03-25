function [] = coupled_orbit_1d(t, sv, masses, J, N)

omega = sv(1:N);
x = sv(N+1:2*N);
R = sv(2*N+1:3*N);
phi = sv(3*N+1:4*N);
tau = sv(4*N+1:5*N);

numelems = 4*N;
% massM = mass_matrix( ); 



end