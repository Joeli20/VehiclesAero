clear all;
close all;
clc;

% Escrit per: Joel Campo, Albert Chacón
% Vehicles Aeroespacials. MUEA.
% Task 1: Matrix structural analysis of a optical mount

load("fe_model.mat");

DoF = 6;
nodes_fix = [10735; 13699; 16620; 19625; 22511; 4747];

% PREALLOCATING
fixnodes = zeros(size(nodes_fix,1)*DoF,3);

for i = 1:size(nodes_fix,1)
    for j = 1:DoF
        fixnodes(j+(DoF*(i-1)),1) = nodes_fix(i);
        fixnodes(j+(DoF*(i-1)),2) = j;
    end
    clear j
end
clear i
