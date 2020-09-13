function v_evolved = exponentiate(v_in,H_in,dt_in,subdivisions)
%% EXPONENTIATE performs time evolution on an inpute state using a Taylor
%% approximation of the input Hamiltonian.
    % v_in is the input state.
    % H_in is the input Hamiltonian.
    % dt is the time step for which time evolution occurs. dt_in is defined
    % as dt multiplied by a factor of -i when using EXPONENTIATE.
    % subdivisions is the number of times the Taylor approximation is
    % used in the time step.

%% Split the time step based on the number of subdivisions.
ddt = dt_in/subdivisions;

%% Define v_evolved to be equal to the time evolution of v_in for the first
%% time subdivision.
v_evolved = v_in + ddt*H_in*v_in + ddt^2/2*H_in*(H_in*v_in) ...
    + ddt^3/6*H_in*(H_in*(H_in*v_in));

%% Continue evolving v_evolved for each time subdivision in dt.
for nt = 1:subdivisions-1
    v_evolved = v_evolved + ddt*H_in*v_evolved + ddt^2/2*H_in ...
        *(H_in*v_evolved) + ddt^3/6*H_in*(H_in*(H_in*v_evolved));
end