clear all;
close all;
clc;

% Escrit per: Joel Campo, Albert Chacón
% Vehicles Aeroespacials. MUEA.
% Task 1: Matrix structural analysis of a optical mount

load("fe_model.mat");

%% TASK 1_1
dimension = 2; %Sobre quina dimensió s'aplica la gravetat

DoF = 6;
nodes_fix = [10735; 13699; 16620; 19625; 22511; 4747];

% PREALLOCATING
fixnodes = zeros(size(nodes_fix,1)*DoF,3);
posicio = zeros(size(nodes_fix,1),1);
in_D = zeros(size(nodes_fix,1)*DoF,1);

% VALORS GENERALS
g = 9.81e3;

% FIXNODES
for i = 1:size(nodes_fix,1)
    for j = 1:DoF
        fixnodes(j+(DoF*(i-1)),1) = nodes_fix(i);
        fixnodes(j+(DoF*(i-1)),2) = j;
    end
    clear j
end
clear i

% INDEX LOCATION
for i = 1:size(nodes_fix,1)
    [posicio(i)] = Pos_Find(nodes_fix(i),DoF);
end
clear i

for i = 1:size(nodes_fix,1)
    for j = 1:DoF
    in_D(j+(DoF*(i-1)),1) = posicio(i) + (j-1); % Dirichlet Index
    end
    clear j
end
clear i

n_tot = 1:size(K,1); % Dummy vector d'1 a Nod*DoF
in_N = setdiff(n_tot,in_D); % Neumann Index

% Calcul de les K segregades
K_DD = K(in_D,in_D);
K_NN = K(in_N,in_N);
K_DN = K(in_D,in_N);
K_ND = K(in_N,in_D);

% Calcul u_D
u_D = fixnodes(:,3);

% Calcul F_N
F = zeros(size(n_tot,2),1);

for i=1:(size(F,1)/6)
    F(dimension+6*(i-1)) = M(dimension+6*(i-1),dimension+6*(i-1))*g;
end

F_N = F(in_N);

% CALCULATIONS

u_N = K_NN\(F_N - K_ND * u_D);
F_D = K_DD * u_D + K_DN * u_N;

F_D = reshape(F_D,[6,6]);

% MASS COMPROVATION

true_mass=0;

for i=1:(size(F,1)/6)
   true_mass = true_mass + M(dimension+6*(i-1),dimension+6*(i-1))*g;
end

calc_mass = sum(F_D(2,:));

error = true_mass+calc_mass;

%% TASK 1_3 EIGENMODES

M_DD = M(in_D,in_D);
M_NN = M(in_N,in_N);
M_DN = M(in_D,in_N);
M_ND = M(in_N,in_D);

N_mod = 5;

[V_r,D_r] = eigs(K_NN,M_NN,N_mod,'smallestabs');
[V_u,D_u] = eigs(K,M,(N_mod+DoF),'smallestabs');

eig_val_r = diag(D_r); % Eigenvalues restricted
eig_val_u = diag(D_u); % Eigenvalues unconstrained

[eig_val_r,ordre_r] = sort(eig_val_r); % Endreçar de menor a major
V_r = V_r(:,ordre_r); % Endreçar amb el mateix ordre

[eig_val_u,ordre_u] = sort(eig_val_u); % Endreçar de menor a major
V_u = V_u(:,ordre_u); % Endreçar amb el mateix ordre

eig_mod_r = V_r(:,1:N_mod); % Eigenmodes collected. En el nostre cas V=eig_mod 
                        % pero preparem per cas general que haguem extret 
                        % més i volguem uns en concret

eig_mod_r2 = zeros(size(K,1),5);
eig_mod_r2(in_N,:) = eig_mod_r;

eig_mod_u = V_u(:,(DoF+1):(N_mod+DoF));


freq_r = sqrt(eig_val_r)/(2*pi); 
freq_u = sqrt(eig_val_u)/(2*pi);

% Preparing to Plot RESTRICTED
u_restr = zeros((size(K,1)/DoF),DoF,N_mod);

k = 1;
for n = 1:N_mod
    for i = 1:size(K,1)/DoF
        for j = 1:DoF
                u_restr(i,j,n) = eig_mod_r2(j+(i-1)*6,n);
        end
    end
end

clear i j n;

fillhdf('template.h5','restricted_1.h5',u_restr(:,:,1));
fillhdf('template.h5','restricted_2.h5',u_restr(:,:,2));
fillhdf('template.h5','restricted_3.h5',u_restr(:,:,3));
fillhdf('template.h5','restricted_4.h5',u_restr(:,:,4));
fillhdf('template.h5','restricted_5.h5',u_restr(:,:,5));

% Preparing to Plot UNCONSTRAINED
u_unconst = zeros((size(K,1)/DoF),DoF,N_mod);

for n = 1:N_mod
    for i = 1:(size(K,1)/DoF)
        for j = 1:DoF
            u_unconst(i,j,n) = eig_mod_u(j+(i-1)*6,n);
        end
    end
end
 
fillhdf('template.h5','unconstrained_1.h5',u_unconst(:,:,1));
fillhdf('template.h5','unconstrained_2.h5',u_unconst(:,:,2));
fillhdf('template.h5','unconstrained_3.h5',u_unconst(:,:,3));
fillhdf('template.h5','unconstrained_4.h5',u_unconst(:,:,4));
fillhdf('template.h5','unconstrained_5.h5',u_unconst(:,:,5));

clear i j n;