clear all;
close all;
clc;

% Escrit per: Joel Campo, Albert Chacón
% Vehicles Aeroespacials. MUEA.
% Task 1: Matrix structural analysis of a optical mount

load("fe_model.mat");

%% TASK 1
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
    F(dimension+6*(i-1)) = M(i,i)*g;
end

F_N = F(in_N);

% CALCULATIONS

u_N = K_NN\(F_N - K_ND * u_D);
F_D = K_DD * u_D + K_DN * u_N;
