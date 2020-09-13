function wave_delay_time_series(L,t,a,m)
%% WAVE_DELAY_TIME_SERIES plots the time series for a wave leaving the
%% interior region in the projected cooling algorithm.
    % L is the lattice size.
    % t is the number time steps after t = 0 on the graph.
    % a is the distance between adjacent lattice points. Typically a = 1.
    % m is the mass of the particle.

%% Construct the Hamiltonian
H = hamiltonian_construction(L,a,m);

%% Define the starting delta.
v = zeros(L,1);
v(1,1) = 1;

%% Perform time evolution and plot each step.
for jj = 0:t
    plot((abs(expm(-i*H*jj)*v)).^2 + 0.2*jj);
    hold on
end