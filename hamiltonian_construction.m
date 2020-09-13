function H = hamiltonian_construction(L,a,m)
%% HAMILTONIAN_CONSTRUCTION creates the lattice Hamiltonian for projected
%% cooling. We set Planck's constant equal to 1.
    % L is the lattice size.
    % a is the distance between adjacent lattice points. Typically a = 1.
    % m is the mass of the particle.

%% Set values of H.
H = zeros(L);
for ii = 1:L
    H(ii,ii) = -2;
end
for ii = 1:L-1
    H(ii+1,ii) = 1;
    H(ii,ii+1) = 1;
end
H(L,1) = 1;
H(1,L) = 1;

%% Multiply H by corresponding constant.
H = -1/(2*m*a^2) * H;
end