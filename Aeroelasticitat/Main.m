clear all;
close all;
clc;

% Written by: Joel Campo, Albert Chacón
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
k_1 = 4000; % N/m
k_2 = 3000; % N/m
k_3 = 4000; % N/m
k_4 = 2000; % N/m
k_5 = 2500; % N/m

% Others
rho = 1.25; %kg/m^3

%% PART B - Structural dynamics modelling
% Previous calculations
m_1 = mu * h_1 * c_1;
m_2 = mu * h_2 * c_1;
m_3 = mu * h_2 * c_2;

Is_1 = (1/12) * m_1 * c_1^2 + m_1 * (0.5 * c_1 - x_1)^2;
Is_2 = (1/12) * m_2 * c_1^2 + m_2 * (0.5 * c_1 - x_1)^2;
Is_3 = (1/12) * m_3 * c_2^2 + m_3 * (0.5 * c_2 - x_2)^2;

% Stiffness Matrix computation
K = zeros(5,5);

K(1,1) = k_3 * x_1^2;
K(1,2) = -1 * (k_3 * x_1^2);

K(2,1) = -1 * (k_3 * x_1^2);
K(2,2) = k_3 * x_1^2 + k_4 * (c_1 - x_1)^2;
K(2,4) = -1 * ((k_4/h_0) * (c_1 - x_1)^2);
K(2,5) = (k_4/h_0) * (c_1 - x_1)^2;

K(3,3) = k_5 * x_2^2;
K(3,4) = -1 * ((k_5/h_0) * x_2^2);
K(3,5) = (k_5/h_0) * x_2^2;

K(4,2) = -1 * ((k_4/h_0) * (c_1 - x_1)^2);
K(4,3) = -1 * ((k_5/h_0) * x_2^2);
K(4,4) = k_1 + ((k_4/(h_0^2)) * (c_1 - x_1)^2) + ((k_5/(h_0^2)) * x_2^2);
K(4,5) = -1 * ((k_4/(h_0^2)) * (c_1 - x_1)^2) - ((k_5/(h_0^2)) * x_2^2);

K(5,2) = (k_4/h_0) * (c_1 - x_1)^2;
K(5,3) = (k_5/h_0) * x_2^2;
K(5,4) = -1 * ((k_4/(h_0^2)) * (c_1 - x_1)^2) - ((k_5/(h_0^2)) * x_2^2);
K(5,5) = k_2 + ((k_4/(h_0^2)) * (c_1 - x_1)^2) + ((k_5/(h_0^2)) * x_2^2);

% Mass Matrix computation
M = zeros(5,5);

M(1,1) = Is_1;
M(1,4) = m_1 * (x_1 - c_1/2);

M(2,2) = Is_2;
M(2,4) = m_2 * (x_1 - c_1/2);

M(3,3) = Is_3;
M(3,5) = m_3 * (x_2 - c_2/2);

M(4,1) = m_1 * (x_1 - c_1/2);
M(4,2) = m_2 * (x_1 - c_1/2);
M(4,4) = m_1 + m_2;

M(5,3) = m_3 * (x_2 - c_2/2);
M(5,5) = m_3;

% EIGENVALUES COMPUTATION
n_eigval = 5;
[mode,value] = eigs(K,M,n_eigval);

freq = diag(sqrt(value));

% Mode shapes (required results)
modeshape = zeros(5);

for i = 1:n_eigval
    modeshape(i,1) = (mode(4,i) + mode(5,i))/2;
    modeshape(i,2) = (mode(4,i) - mode(5,i))/h_0;
    modeshape(i,3) = (mode(1,i) - mode(2,i));
    modeshape(i,4) = (mode(2,i) - modeshape(i,2));
    modeshape(i,5) = (mode(3,i) - modeshape(i,2));
end

MODE_1 = mode(:,1);
MODE_2 = mode(:,2);
MODE_3 = mode(:,3);
MODE_4 = mode(:,4);
MODE_5 = mode(:,5);

% Plotting (if needed)
n_mode = MODE_1;
plot_structure(n_mode(4),n_mode(5),n_mode(1),n_mode(2),n_mode(3));

%% PART C - Aeroelastic analysis
% Force Vector computation (matrix)
S = zeros(5,3);

S(1,1) = h_1 * (x_1 - c_1/4);

S(2,2) = h_2 * (x_1 - c_1/4);

S(3,3) = h_2 * (x_2 - c_2/4);

S(4,1) = h_1;
S(4,2) = h_2;

S(5,3) = h_2;

% Aerodynamics Matrix computation

C = zeros(3,3);

C(1,1) = pi * c_1;

C(2,2) = 1;
C(2,3) = 1;

C(3,2) = 1;
C(3,3) = 1;

% C11 = 1*pi*c1;
% C12 = 0;
% C13 = 0;
% C21 = 0;
% C22 = -1/((c2/(4*pi*(0.75*c1-0.25*c2-x1-h0+x2)*(0.75*c2-0.25*c1+x1+h0-x2)))-1/(pi*c1));
% C23 = ((1*c2)/(2*(0.75*c1-0.25*c2-x1-h0+x2)))/((c2/(4*pi*(0.75*c1-0.25*c2-x1-h0+x2)*(0.75*c2-0.25*c1+x1+h0-x2)))-1/(pi*c1));
% C31 = 0;
% C32 = ((1*c1)/(2*(0.75*c2-0.25*c1+x1+h0-x2)))/((c1/(4*pi*(0.75*c1-0.25*c2-x1-h0+x2)*(0.75*c2-0.25*c1+x1+h0-x2)))-1/(pi*c2));
% C33 = -1/((c1/(4*pi*(0.75*c1-0.25*c2-x1-h0+x2)*(0.75*c2-0.25*c1+x1+h0-x2)))-1/(pi*c2));

C = rho * C;

C_calcul = zeros(3,5);

C_calcul(1,1) = 1;
C_calcul(2,2) = 1;
C_calcul(3,3) = 1;

% Divergence COUPLED SYSTEM
U_inf = 1;

l = rho * U_inf^2 * C * C_calcul;
A = U_inf^2 * S * C * C_calcul;

% Divergence Speed Static Cond (a)
[mode_A,value_A] = eigs(K,A,n_eigval);

U_Div_eigen = diag(sqrt(value_A));
U_Div = real(U_Div_eigen(3));

% Potential control reversal conditions (b)

