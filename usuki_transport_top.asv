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
%surf(V);
%axis tight
%camlight
%lighting phong
%shading interp
%xlabel('x (nm)');
%ylabel('y (nm)');
%zlabel('t (J)');
%hold off;

e = 1.602E-19;
Ef_steps = 25;
Ef = linspace(0,50E-3*e,Ef_steps);
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

I = eye(Ny);
Zero = zeros(Ny);

% Define H^(0) matrix (Ny-by-Ny)
t_vector = zeros(Ny,1);
t_vector(:) = t;
neg_t_vector = zeros(Ny-1,1);
neg_t_vector(:) = -t;
H_0 = diag(4*t + V(:,1)); % V(:,1) ? 
H_0p1 = diag(neg_t_vector,1);
H_0n1 = diag(neg_t_vector,-1);
H_0 = H_0 + H_0p1;
H_0 = H_0 + H_0n1;

H = zeros(Nx, Ny, Ny);
%H = zeros(Nx);
for i = 1:Nx
   H(i,:,:) = diag(4*t + V(:,i)) + H_0p1 + H_0n1;
   %H(i,:,:) = H(i,:,:) + H_0p1;
   %H(i,:,:) = H(i,:,:) + H_0n1;
end


H_0_0p1 = diag(t_vector);
H_0_0n1 = diag(t_vector);


T_11 = zeros(Ny); % 0
T_12 = eye(Ny); % I
T_21 = -inv(H_0_0p1)*H_0_0n1; % No B = -I
%Temp = Ef.*T_12;
T_22 = (inv(H_0_0p1))*(H_0 - Ef(25).*T_12);

%T = zeros(Ef_steps, Nx, Ny, Ny);
T22 = zeros(Nx, Ny, Ny);
T = zeros(Nx, 2*Ny, 2*Ny);
for i = 1:Nx
    %for j = 1:Ef_steps
        temp1 = H(i,:,:);
        temp2 = reshape(temp1,Ny,Ny);
        T22(i,:,:) = (inv(H_0_0p1))*(temp2 - Ef(25).*T_12);
        temp1 = T22(i,:,:);
        temp2 = reshape(temp1,Ny,Ny);
        T(i,:,:) = [Zero I ; T_21 temp2]; 
    %    T(j,i) = [Zero I ; T_21 T22]; 
    %end
end


T0 = [T_11, T_12 ; T_21, T_22];

[T0_new, Dsort] = getEigsAndSortModes(T0); % mode-to-slice conversion (new T0)

T11_new = T0_new(1:Ny, 1:Ny);
T21_new = T0_new(Ny+1:2*Ny, 1:Ny);
T12_new = T0_new(1:Ny, Ny+1:2*Ny);
T22_new = T0_new(Ny+1:2*Ny, Ny+1:2*Ny);

C_1_0 = eye(Ny);
C_2_0 = zeros(Ny);

P_2_0_1 = T21_new * C_2_0;
P_2_0_2 = T22_new;
P_2_0 = P_2_0_2^-1;
P_1_0 = -1*P_2_0 * T21_new*C_1_0;

P_mat_0 = [I Zero ; P_1_0 P_2_0];
C_mat_0 = [C_1_0 C_2_0 ; Zero I];

C_mat = zeros(Nx+2, 2*Ny, 2*Ny);

C_mat(1,:,:) = T0_new * C_mat_0 * P_mat_0;

P_2 = zeros(Nx, Ny, Ny);
P_1 = zeros(Nx, Ny, Ny);
C_1 = zeros(Nx, Ny, Ny);
C_2 = zeros(Nx, Ny, Ny);



for i = 1:Nx  
  fprintf('%d/%d\n',i,Nx);
  
  C_mattemp1 = C_mat(i,:,:);
  C_mattemp2 = reshape(C_mattemp1,2*Ny,2*Ny);
  
  C_mat_temp = C_mattemp2;
  C_1(i,:,:) = C_mat_temp(1:Ny,1:Ny);
  C_2(i,:,:) = C_mat_temp(Ny+1:2*Ny, Ny+1:2*Ny);
  
  C_1temp1 = C_1(i,:,:);
  C_1temp2 = reshape(C_1temp1,Ny,Ny);
  
  C_2temp1 = C_2(i,:,:);
  C_2temp2 = reshape(C_2temp1,Ny,Ny);
  
  %T22temp1 = T_22(i,:,:);
  %T22temp2 = reshape(T22temp1,Ny,Ny);
  %testP2 = reshape(P_2,Ny,Ny);
  P_2(i,:,:) = inv(T_21*C_1temp2 + T_22);
  P_1(i,:,:) = -1*P_2(i)*T_21*C_1temp2;
  
  P1temp1 = P_1(i,:,:);
  P1temp2 = reshape(P1temp1,Ny,Ny);
  
  P2temp1 = P_2(i,:,:);
  P2temp2 = reshape(P2temp1,Ny,Ny);
  
  P_mat = [I Zero ; P1temp2 P2temp2];
  C_mat(i,:,:) = [C_1temp2 C_2temp2 ; Zero I];
  
  C_mattemp1 = C_mat(i,:,:);
  C_mattemp2 = reshape(C_mattemp1,2*Ny,2*Ny);
  
  Ttemp1 = T(i,:,:);
  Ttemp2 = reshape(Ttemp1,2*Ny,2*Ny);
  
  C_mat(i+1,:,:) = Ttemp2*P_mat*C_mattemp2;
  
end




% Mode-to-Slice (0 -> 1)

% Slice-to-Slice (1 -> Nx)


% Slice-to-Mode (Nx+1, Nx+2, collect info at Nx+2)


% -----------------
% Problem 1
% -----------------

% a) Plot the conductance of the QPC as a fucntion of the Fermi Energy
% from 0 <= Ef <= 50 meV for the QPC being 50 nm wide.



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
