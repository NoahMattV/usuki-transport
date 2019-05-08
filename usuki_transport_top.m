% Top module for usuki-transport repo
% Noah Van Der Weide 4/22/2019
% UW Madison - ECE 845

% 2D transport, identifying slice-by-slice along x

% Problem 1: Quantum Point Contact
% GaAs Nanowire
clc;
close all;
clear; 
global Nx;
global Ny;

QW_width = 500; % y
QW_length = 1000; % x
QPC_width = 50; % along y
QPC_length = 50; % along x

QW_w = QW_width*1E-9;
QW_l = QW_length*1E-9;
QPC_w = QPC_width*1E-9;
QPC_l = QPC_length*1E-9;

V = qpcprofile(QW_width,QW_length,QPC_width,QPC_length);

%hold on;
surf(V);
axis tight
camlight
lighting phong
shading interp
xlabel('x (nm)');
ylabel('y (nm)');
zlabel('t (J)');
%hold off;

e = 1.602E-19;

W = 500E-9; % 500 nm wide
L = 1000E-9; % 1000 nm long
QPC_L = 50E-9; % QPC contact is 50 nm long

%Ef = 50E-3*e; % Fermi energy [J] (part b)
a = 2E-9; % spacing
h = 6.626E-34; % [J-s]
hbar = h/(2*pi); % Reduced planck's constant [J-s]
m0 = 9.11E-31; % [kg]
m = 0.067*m0; % [kg]

t = (hbar^2)/(2*m*a^2); % Hopping Energy [J]

phi_l = zeros(Ny,1); % column of phi_l,1 to phi_l,Ny

% Define H^(0) matrix (Ny-by-Ny)
t_vector = zeros(Ny-1,1);
t_vector(:) = t;
H_0 = diag(4*t+ V(:,1));
H_0p1 = diag(t_vector,1);
H_0n1 = diag(t_vector,-1);
H_0 = H_0 + H_0p1;
H_0 = H_0 + H_0n1;

C_1 = zeros(Ny,1);
C_2 = zeros(Ny,1);

T_11 = zeros(Ny,Nx); % 0
T_12 = zeros(Ny,Nx); % I hat
T_21 = zeros(Ny,Nx); % -I hat
T_22 = zeros(Ny,Nx); % H hat ^-1 l,l+1 (H hat l - E*I hat)

P_2 = zeros(Ny,1);
P_1 = zeros(Ny,1);

%C_1(1,1) = -I_hat;
%C_2(1,1) = 0;

for i = 1:Ny
  P_2(i,1) = (T_21(i,1)*C_2(i,1) + T_22(i,1))^-1; %invert matrix?
  P_1(i,1) = -1*P_2(i,1)*T_21(i,1)*C_1(i,1);
end

% Mode-to-Slice (0 -> 1)

% Slice-to-Slice (1 -> Nx)


% Slice-to-Mode (Nx+1, Nx+2, collect info at Nx+2)


% -----------------
% Problem 1
% -----------------

% a) Plot the conductance of the QPC as a fucntion of the Fermi Energy
% from 0 <= Ef <= 50 meV for the QPC being 50 nm wide.
Ef = zeros(100,1);


% b) For the Fermi energy of 50 meV, vary the width of the QPC between
% 10 and 100 nm, and plot the conductance as a function of the QPC width.


% c) For a value of the Fermi energy such that 3 modes are propagating through
% the QPC of width 50 nm, plot the carrier density throughout the structure.


% d) For the Fermi energy of 50 meV (what is the 2D sheet density ns
% corresponding to that Fermi energy?), vary the position of the bottom of the
% potential well in the QPC as bottom gate E = 60meV − eV ,
% where 0 ≤ Vgate < 60 mV . Plot the conductance as a function of the
% applied gate bias. QPC width is 50 nm.


% e) Plot the carrier density throughout the structure for a gate bias
% corresponding to 1 propagating mode, and then for a bias corresponding
% to 3 propagating modes.


% -----------------
% Problem 2
% -----------------
