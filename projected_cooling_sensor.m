function [first_initial_overlap,overlaps] = projected_cooling_sensor(L,D, ...
    maxeig,nrepeat,H_obj)
%% PROJECTED_COOLING_SENSOR performs the projected cooling sensor algorithm
%% multiple times on a DxD Hamiltonian over a lattice of size L. The outputs
%% are (1) the initial overlap with the randomized state and (2) the overlaps
%% after each iteration of the algorithm.
    % L is the reservoir size.
    % D is the dimension of the Hamiltonian.
    % maxeig is the maximum energy eigenvalue.
    % nrepeat is the number of iterations the algorithm is performed.
    % H_obj is the Hamiltonian of interest. Set H_obj = 0 for arbitrary
    % DxD Hamiltonians.

%% Set the run time proportional to L, and break each time step dt = 0.5 into
%% 10 subdivisions.
dt = 0.5;
Lt = floor(L/dt*0.8);
subdivisions = 10;

%% Define the reservoir Hamiltonian H_res.
r = [0:L-1];
H_res = - sparse(mod(r+1,L)+1,r+1,0.5) - sparse(mod(r-1,L)+1,r+1,0.5);

%% Create arbitrary Hamiltonian if H_obj is not previously defined.
if H_obj == 0
    H_obj = rand(D,D) + i*rand(D,D);
    H_obj = 0.5*(H_obj + H_obj');
end

%% Shift and scale H_obj so that the maximum energy of H_obj_new is maxeig.
list = sort(eig(H_obj));
H_object_new = H_obj-(list(1)+list(2))/2*eye(D);
scale = max(eig(H_object_new));
H_object_new = H_object_new/scale*maxeig;

%% Define system Hamiltonian H as the tensor product between H_res and H_obj_new.
H = sparse(D*L,D*L);
for ii = 0:D-1 
    H(ii*L+[0:L-1]+1,ii*L+[0:L-1]+1) = H_res;
end
for jj = 0:D-1
    for ii = 0:D-1
        H(ii*L+1,jj*L+1) = H(ii*L+1,jj*L+1) + H_object_new(ii+1,jj+1)/2.0;
        H(ii*L+2,jj*L+1) = H(ii*L+2,jj*L+1) + H_object_new(ii+1,jj+1)/2.0;
        H(ii*L+1,jj*L+2) = H(ii*L+1,jj*L+2) + H_object_new(ii+1,jj+1)/2.0;
        H(ii*L+2,jj*L+2) = H(ii*L+2,jj*L+2) + H_object_new(ii+1,jj+1)/2.0;
    end
end

%% Calculate the exact ground state to check for overlap later.
[vv,dd] = eig(H_obj);
[~,ord] = sort(diag(dd));
index = find(abs(vv(:,ord(1))) == max(abs(vv(:,ord(1)))));
vv_exact = vv(:,ord(1))/vv(index,ord(1));

%% Define the random values of the initial state.
vobj_init = zeros(D,1);
for ii = 0:D-1
    vobj_init(ii+1) = (rand-0.5) + i*(rand-0.5);
end
v_init = zeros(D*L,1);

%% Perform the projected cooling sensor algorithm nrepeat times.
for ntrial = 1:nrepeat

    %Define initial state v_init using random values from vobj_init.
    for ii = 0:D-1
        v_init(ii*L+1) = vobj_init(ii+1,1);
        v_init(ii*L+2) = vobj_init(ii+1,1);       
    end
    
    %Perform time evolution for first time step dt.
    psi1(:,1) = exponentiate(v_init,H,-i*dt,subdivisions);    

    %Perform time evolution for remaining time steps.
    for nt = 1:Lt
        psi1(:,nt+1) = exponentiate(psi1(:,nt),H,-i*dt,subdivisions);
    end

    %Graph the convergence of the state(s) vs. time.
    figure(ntrial)
    hold on
    for ii = 0:D-1
        plot(dt*[1:Lt],log(abs(psi1(ii*L+1,1:Lt)+psi1(ii*L+2,1:Lt)).^2))
    end
    xlabel('{\it t}','FontSize',14)
    ylabel('ln({\it|v_{PCS} |^2})','FontSize',14)
    hold off
 
    %Find the reconstructed ground state and display with trial number.
    v_PC = zeros(D,1);
    for ii = 0:D-1
        v_PC(ii+1,1) = psi1(ii*L+1,Lt) + psi1(ii*L+2,Lt);
    end
    index = find(abs(v_PC) == max(abs(v_PC)));
    v_PC = v_PC/v_PC(index);
    disp(ntrial)
    disp([v_PC vv_exact])
    
    %Calculate initial and final overlap values and display for each trial.
    initial_overlap(ntrial,1) = (abs(vv_exact'*vobj_init))^2 ...
        /((vv_exact'*vv_exact)*(vobj_init'*vobj_init));
    final_overlap(ntrial,1) = (abs(v_PC'*vv_exact))^2 ...
        /((vv_exact'*vv_exact)*(v_PC'*v_PC));
    disp(initial_overlap(ntrial,1))
    disp(final_overlap(ntrial,1))
    
    %Set initial state as reconstructed state for next iteration.
    vobj_init = v_PC;
end

%% Define initial overlap for algorithm and final overlaps for each ntrial.
first_initial_overlap = initial_overlap(1);
overlaps = final_overlap;

%% Graph log error of overlap vs. iteration.
figure(nrepeat+1)
plot(log(1-final_overlap))
xlabel('Iteration','FontSize',14)
ylabel('ln(1-{\it O})','FontSize',14)
end