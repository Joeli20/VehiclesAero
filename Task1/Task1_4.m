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
nodes_fix = [4747 10735 13699 16620 19625 22511];

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

%% TASK 4

Nod_ref = 1305;

F_ext = zeros(size(K_NN,1),1);
F_ext((Nod_ref-1)*6 + 1,1) = g;

M_DD = M(in_D,in_D);
M_NN = M(in_N,in_N);
M_DN = M(in_D,in_N);
M_ND = M(in_N,in_D);

dump_rat = 0.02; % 2%
freq = 0:2000;
omega = 2*pi*freq;

N_mod = 5; % Review

[V,D] = eigs(K_NN,M_NN,N_mod,'smallestabs');

eig_val = diag(D); % Eigenvalues restricted

[eig_val,ordre] = sort(eig_val); % Endreçar de menor a major
V = V(:,ordre); % Endreçar amb el mateix ordre

eig_mod = V(:,1:N_mod); % Eigenmodes collected.

freq_eig = sqrt(eig_val); %OJO

M_mod = eig_mod' * M_NN * eig_mod; % Transp * M * eig

K_DD = (1 + dump_rat*sqrt(-1)) * K_NN;
K_mod = eig_mod' * K_DD * eig_mod; % Transp * K * eig

f_mod = eig_mod' * F_ext;

for i = length(omega)
        Q_mod = -(omega(i)^2).* M_mod + K_mod;
        X_i = Q_mod \ f_mod;
        X_barret = eig_mod * X_i;
        X_barret_total = zeros(size(M,1), 1);
        X_barret_total(in_N) = X_barret;
        X_total(i) = X_barret_total((Nod_ref-1)*6+1);
end

X_mod = abs(X_total);
X_ang = angle(X_total);

figure
hold on
plot(freq, X_mod);
xlabel('Freqüència (Hz)')
ylabel('Amplitud (mm)')
title('Mòduls')
legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5')
grid on
hold off

figure
hold on
plot(freq, -1*X_ang);
xlabel('Freqüència (Hz)')
ylabel('Fase (rad)')
title('Angles')
legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5')
grid on
hold off