clear all;
close all;
clc;

% Escrit per: Joel Campo, Albert Chacón
% Vehicles Aeroespacials. MUEA.
% Task 1: Matrix structural analysis of a optical mount

load("fe_model.mat");

%% TASK 1
dimension = 2; %Dimensio sobre la qual s'aplica el desplaçament del node a causa de la shim
support2disp = 1:6; %Suport sobre el qual s'aplica el desplaçament
refNode = 1305;
refNodeDof = 6*(refNode-1)+1;

DoF = 6;
nodes_fix = [10735; 13699; 16620; 19625; 22511; 4747];

% PREALLOCATING
fixnodes = zeros(size(nodes_fix,1)*DoF,3);
posicio = zeros(size(nodes_fix,1),1);
in_D = zeros(size(nodes_fix,1)*DoF,1);

% VALORS GENERALS
g = 0;

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

u_RefNode = zeros(size(support2disp));

for i=1:length(support2disp)

    fixnodes(:,3) = 0;

    % Imposició del desplaçament causat per una shim a una pota
    fixnodes(6*(support2disp(i)-1) + dimension,3) = 1;

    % Calcul u_D
    u_D = fixnodes(:,3);

    % Calcul F_N
    F = zeros(size(n_tot,2),1);

    for j=1:(size(F,1)/6)
        F(dimension+6*(j-1)) = M(dimension+6*(j-1),dimension+6*(j-1))*g;
    end

    F_N = F(in_N);

    % CALCULATIONS

    u_N = K_NN\(F_N - K_ND * u_D);
    F_D = K_DD * u_D + K_DN * u_N;

    % Reference Node Displacement

    posRef = find(in_N==refNodeDof);
    u_RefNode(i,:) = u_N(posRef:(posRef+5));
end
clear i

% System resolution

TargetDisp = [0; 0; 0; 0.0005; 0; -0.0002];

u_RefNode = transpose(u_RefNode);

Shims = u_RefNode\TargetDisp;

% Result Comprovation

u_D2 = zeros(length(nodes_fix)*DoF,1);
F2 = zeros(size(n_tot,2),1);

for i=1:length(support2disp)
    u_D2(6*(support2disp(i)-1) + dimension) = Shims(i);
end
clear i

for j=1:(size(F2,1)/6)
    F2(dimension+6*(j-1)) = M(dimension+6*(j-1),dimension+6*(j-1))*g;
end
clear j

F2_N = F2(in_N);

u_N2 = K_NN\(F2_N - K_ND * u_D2);
F_D2 = K_DD * u_D2 + K_DN * u_N;

u_RefNode2 = u_N2(posRef:(posRef+5));