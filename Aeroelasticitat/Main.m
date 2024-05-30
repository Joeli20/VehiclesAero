clear all;
close all;
clc;

% Escrit per: Joel Campo, Albert Chacón
% Vehicles Aeroespacials. MUEA.
% Aeroelasticity Project

%% Hypotesis
% - Weight loads are negligible compared to aerodynamic loads.
% - All spring deformations are small (small angles hypothesis holds).
% - Panels are flat and thin.
% - Incompressible and inviscid flow (consider an air density of 𝜌𝜌= 1.25 kg/m3).
% - Induced velocities in the spanwise direction are negligible (strip theory).

%% INPUT DATA
% Geometric
c_1 = 0.65; % m
c_2 = 0.50; % m
h_0 = 1.50; % m
h_1 = 1.25; % m
h_2 = 1.75; % m
x_1 = 0.22; % m
x_2 = 0.20; % m
s   = 0.17; % m

% Material
mu  = 2.0; % kg/m^2
K_1 = 4.0; % kN/m
K_2 = 3.0; % kN/m
K_3 = 4.0; % kN/m
K_4 = 2.0; % kN/m
K_5 = 2.5; % kN/m

% Others
rho = 1.25; %kg/m^3

%% PART A - Aerodynamic modelling

%% PART B - Structural dynamics modelling

%% PART C - Aeroelastic analysis
