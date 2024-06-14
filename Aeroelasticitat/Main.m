clc;
clear all;
close all;


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
% n_mode = MODE_1;
% plot_structure(n_mode(4),n_mode(5),n_mode(1),n_mode(2),n_mode(3));

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

C(2,2) = (-1)/((c_2/(4 * pi * ((3/4) * c_1 - (1/4) * c_2 - x_1 - h_0 + x_2) * ((3/4) * c_2 - (1/4) * c_1 + x_1 + h_0 - x_2))) - 1/(pi * c_1));
C(2,3) = ((1 * c_2)/(2 * ((3/4) * c_1 - (1/4) * c_2 - x_1 - h_0 + x_2)))/((c_2/(4 * pi * ((3/4) * c_1 - (1/4) * c_2 - x_1 - h_0 + x_2) * ((3/4) * c_2 - (1/4) * c_1 + x_1 + h_0 - x_2))) - 1/(pi * c_1));

C(3,2) = ((1 * c_1)/(2 * ((3/4) * c_2 - (1/4) * c_1 + x_1 + h_0 - x_2)))/((c_1/(4 * pi * ((3/4) * c_1 - (1/4) * c_2 - x_1 - h_0 + x_2) * ((3/4) * c_2 - (1/4) * c_1 + x_1 + h_0 - x_2))) - 1/(pi * c_2));
C(3,3) = (-1)/((c_1/(4 * pi * ((3/4) * c_1 - (1/4) * c_2 - x_1 - h_0 + x_2) * ((3/4) * c_2 - (1/4) * c_1 + x_1 + h_0 - x_2))) - 1/(pi * c_2));

C = rho * C;

C_coupling = zeros(3,5);

C_coupling(1,1) = 1;
C_coupling(2,2) = 1;
C_coupling(3,3) = 1;

% Divergence COUPLED SYSTEM
U_inf = 1;

l = rho * U_inf^2 * C * C_coupling;
A = S * C * C_coupling;

% Divergence Speed Static Cond (a)
[mode_A, value_A] = eigs(K,A,n_eigval);

U_Div_eigen = sqrt(diag(value_A));
U_Div = real(U_Div_eigen(3));

% Potential control reversal conditions (b)
% Small calculations
C_delta = 1.3149;
delta_l_1 = 0.784375;
delta_l_f = 0.530875;
xp_delta_1 = (c_1 - 3/4 * s)/(4 - x_1);
%xp_delta = x_1 - ((c_1 - s)/(4 * delta_l_1 + (c_1 - 3/4 * s) * delta_l_f));
xp_delta = x_1 - ((c_1 - s)/4 * delta_l_1 + (c_1 - 3 * s/4) * delta_l_f);

% Coupled System resolution
S_delta = zeros(5,1);

S_delta(1) = h_1 * xp_delta;
S_delta(4) = h_1;

f_delta = C_delta * S_delta;

D_ini = zeros(1,3);
D_ini(1) = C(1,1)/rho;
D_ini(2) = (C(2,2) + C(2,3))/rho;
D_ini(3) = (C(3,2) + C(3,3))/rho;
D = D_ini * C_coupling;

% Delta Lift Ratio Computation
U_discret = 0:0.01:U_Div;
Delta_Lift_Ratio = zeros(length(U_discret),1);

for i = 1:length(U_discret)
    Delta_Lift_Ratio(i) = 1 + (U_discret(i)^2/C_delta) * (D/(K - U_discret(i)^2 * A) * f_delta);
end

% For better plotting (calculated after the first plotting)
U_r = [12 12]; % Matching first U_R = 0 --> Aproximately
DeltaL_U_r = [-200 200]; % Matching final asymptote
U_d = [U_Div U_Div];

% Required plot
figure 
hold on
plot(U_d, DeltaL_U_r,'--')
plot(U_r, DeltaL_U_r,'--')
plot(U_discret, Delta_Lift_Ratio)
xlabel('U_{\infty}')
ylabel('\DeltaL''/\DeltaL')
ylim([-100 100])
legend('U_D','U_R')
hold off
grid on

% Flutter (c)











B11 = (pi*c1^2/4)*rho;
B22 = (pi*c1^2/4)*rho;
B33 = (pi*c2^2/4)*rho;
b1 = x1-c1*3/4;
b2 = x1-c1*3/4;
b3 = x2-c2*3/4;

A0 = A;
B0 = [b1 0 0 1 0;
    0 b2 0 1 0;
    0 0 b3 0 1];

B1 = [B11 0 0 0 0;
    0 B22 0 0 0;
    0 0 B33 0 0];


S1 = [h1*(x1-3*c1/4) S12 S13;
    S21 h2*(x1-3*c1/4) S23;
    S31 S32 h2*(x2-3*c2/4);
    S41 S42 S43;
    S51 S52 S53];


U_range = 0:0.01:100;  
A1 = S1*B1 - S*C*B0;
A = zeros(10,10,length(U_range));
B = zeros(10,10,length(U_range));
eigsReal = zeros(10,length(U_range));
maxRealEig = zeros(length(U_range),1);
eigsImaginary = zeros(length(U_range),1);

YPlotReal = zeros(10, length(U_range));
YPlotIm = zeros(10, length(U_range));

for i = 1:length(U_range)
    A(1:5,1:5,i) = K-U_range(i)^2*A0;
    A(6:10,1:5,i) = zeros(5);
    A(1:5,6:10,i) = zeros(5);
    A(6:10,6:10,i) = eye(5);
    
    B(1:5,1:5,i) = U_range(i)*A1;
    B(6:10,1:5,i) = eye(5);
    B(1:5,6:10,i) = -M;
    B(6:10,6:10,i) = zeros(5);
    
    [MODES3, EIGENVAL3] = eig(A(:,:,i),B(:,:,i));
    eigsReal(:,i) = real(diag(EIGENVAL3));
    [maxRealEig(i),pos] = max(eigsReal(:,i)); 
    diagEigen = diag(EIGENVAL3);
%     eigsImaginary(i) = imag(diagEigen(pos));

    YPlotReal(:,i) = (eigsReal(:,i)*c2)/(2*U_range(i));
    YPlotIm(:,i) = imag(diag(EIGENVAL3))./(2*pi);
end

flutter = 6.23;

figure 
hold on
xline(flutter,'--')
plot(U_range,YPlotReal,'.')
grid on
legend('Flutter')
xlabel('U_{\infty}')
ylabel('p^Rc/2U_{\infty} ')
hold off

figure 
hold on
plot(U_range,YPlotIm,'.');
xlabel('U_{\infty}')
ylabel('p^I/2\pi ')
grid on
hold off
