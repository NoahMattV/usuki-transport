% Top module for usuki-transport repo
% Noah Van Der Weide 4/22/2019
% UW Madison - ECE 845

% 2D transport, identifying slice-by-slice along x

% Problem 1: Quantum Point Contact
% GaAs Nanowire
clc;
close all;
format long;
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
Ef_steps = 10;
Ef = linspace(0,(50E-3)*e,Ef_steps); % start at 1
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

for i = 1:Nx
    H(i,:,:) = diag(4*t + V(:,i)) + H_0p1 + H_0n1;
    %H(i,:,:) = H(i,:,:) + H_0p1;
    %H(i,:,:) = H(i,:,:) + H_0n1;
end

H_0_0p1 = diag(t_vector);
H_0_0n1 = diag(t_vector);
G = zeros(Ef_steps,1);

for f = 10:Ef_steps

T_11 = zeros(Ny); % 0
T_12 = eye(Ny); % I
T_21 = -H_0_0p1\H_0_0n1; % No B = -I
T_22 = H_0_0p1\(Ef(f)*I - H_0); % for 0th slice

%T = zeros(Ef_steps, Nx, Ny, Ny);
T22 = zeros(Nx, Ny, Ny);
T = zeros(Nx, 2*Ny, 2*Ny);

% Defining T slices for 1 to Nx
%for i = 1:Nx
%        Htemp = reshape(H(i,:,:),Ny,Ny);
%        T22(i,:,:) = ((H_0_0p1)^-1)*(Ef(f)*I - Htemp); %Htemp - Ef(f)*I
%        T22temp = reshape(T22(i,:,:),Ny,Ny);
%        T(i,:,:) = [Zero I ; T_21 T22temp]; 
%end


T0 = [T_11 T_12 ; T_21 T_22];

% Eigenvalues and Eigenvectors for T0
%clc;
fprintf('Getting Eigs and Sorting Modes...');
[T0_new, Dsort, fp_modes] = getEigsAndSortModes(T0); % mode-to-slice conversion (new T0)

T11_new = T0_new(1:Ny, 1:Ny);
T12_new = T0_new(1:Ny, Ny+1:2*Ny);
T21_new = T0_new(Ny+1:2*Ny, 1:Ny);
T22_new = T0_new(Ny+1:2*Ny, Ny+1:2*Ny);

C_1_0 = eye(Ny);
C_2_0 = zeros(Ny);

%P_2_0_1 = T21_new * C_2_0;
%P_2_0_2 = T22_new;
P_2_0 = inv(T21_new * C_2_0 + T22_new);
P_1_0 = -P_2_0 * T21_new * C_1_0;

P_mat_0 = [I Zero ; P_1_0 P_2_0];
C_mat_0 = [C_1_0 C_2_0 ; Zero I];

C_mat = zeros(Nx+2, 2*Ny, 2*Ny);

C_mat(1,:,:) = T0_new * C_mat_0 * P_mat_0;

P_2 = zeros(Nx+2, Ny, Ny);
P_1 = zeros(Nx+2, Ny, Ny);
C_1 = zeros(Nx+2, Ny, Ny);
C_2 = zeros(Nx+2, Ny, Ny);

for i = 1:Nx  
  C_mat_temp = reshape(C_mat(i,:,:),2*Ny,2*Ny);
  %C_1(i,:,:) = C_mat_temp(1:Ny,1:Ny);
  %C_2(i,:,:) = C_mat_temp(1:Ny, Ny+1:2*Ny);
  C_1temp = C_mat_temp(1:Ny,1:Ny);
  C_2temp = C_mat_temp(1:Ny, Ny+1:2*Ny);
  
  %Htemp = reshape(H(i,:,:),Ny,Ny);
  Htemp = diag(4*t + V(:,i)) + H_0p1 + H_0n1;
  T22temp = H_0_0p1\(Ef(f)*I - Htemp); %Htemp - Ef(f)*I?
  Ttemp = [Zero I ; T_21 T22temp]; 

  P_2(i,:,:) = inv(T_21*C_2temp + T22temp);
  P2temp = reshape(P_2(i,:,:),Ny,Ny);
  P_1(i,:,:) = -P2temp*T_21*C_1temp;
  P1temp = reshape(P_1(i,:,:),Ny,Ny);  
  P_mat = [I Zero ; P1temp P2temp];
  
  C_mat(i+1,:,:) = Ttemp*C_mat_temp*P_mat;
  
end

% Slice-to-Mode (Nx+1, Nx+2, collect info at Nx+2)
C_mat_Nx1 = reshape(C_mat(Nx+1,:,:),2*Ny,2*Ny);
C_1(Nx+1,:,:) = C_mat_Nx1(1:Ny,1:Ny);
C_2(Nx+1,:,:) = C_mat_Nx1(1:Ny, Ny+1:2*Ny);
C_1_Nx1 = reshape(C_1(Nx+1,:,:),Ny,Ny);
C_2_Nx1 = reshape(C_2(Nx+1,:,:),Ny,Ny);

T12_Nx1 = inv(T21_new);
T22_Nx1 = -T11_new*T12_Nx1;
T_Nx1 = [Zero T12_Nx1 ; I T22_Nx1];
T11_Nx1 = T_Nx1(1:Ny, 1:Ny);
T21_Nx1 = T_Nx1(Ny+1:2*Ny, 1:Ny);

P_2_Nx1 = inv(T21_Nx1*C_2_Nx1 + T22_Nx1);
P_1_Nx1 = -P_2_Nx1*T21_Nx1*C_1_Nx1;

P_mat_Nx1 = [I Zero ; P_1_Nx1 P_2_Nx1]; 

C_mat_Nx2 = T_Nx1*C_mat_Nx1*P_mat_Nx1;

C_mat_Nx2temp = reshape(C_mat(Nx+1,:,:),2*Ny,2*Ny);
C_1(Nx+2,:,:) = C_mat_Nx2temp(1:Ny,1:Ny);
C_1_Nx2 = reshape(C_1(Nx+2,:,:),Ny,Ny); % t, matrix of transmission amplitude

Tmn = zeros(size(fp_modes,1)); % Transmission coefficient matrix (only forward propagating modes)
mn_length = size(fp_modes,1);
k = zeros(mn_length, mn_length);
v = zeros(mn_length, mn_length);
G = zeros(mn_length, mn_length);
tmn = zeros(mn_length, mn_length);
for i = 1:mn_length
    for j = 1:mn_length
        tmn(i,j) = C_1_Nx2(i,j); % transmission amplitude matrix with only forward prop modes
    end
end

Tn = inv(T0_new);
[T0_dummy, Dsort_dummy, fp_modes_out] = getEigsAndSortModes(Tn);

for m = 1:mn_length
    for n = 1:mn_length
        kn = log(fp_modes(n))/(1i*a);
        km = log(fp_modes_out(m))/(1i*a);
        vm = hbar*km/m;
        vn = hbar*kn/m;
        Tmn(m,n) = (vm/vn) * abs(tmn(m,n))^2;
        
    end
end

TmnSum = 0;
for m = 1:mn_length
    for n = 1:mn_length
        TmnSum = TmnSum + Tmn(m,n);
    end
end

G(f,1) = (2*e^2/h)*TmnSum;
clc;
fprintf('%d/%d\n',f,Ef_steps);
end %Ef loop


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
