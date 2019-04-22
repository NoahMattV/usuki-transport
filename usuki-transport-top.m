% Top module for usuki-transport repo
% Noah Van Der Weide 4/22/2019
% UW Madison - ECE 845

% 2D transport, identifying slice-by-slice along x

Nx = 100; % Points along x
Ny = 100; % Points along y

e = 1.602E-19;
h = 6.626E-34; % J-s
hbar = h/(2*pi); % Reduced planck's constant
m0 = 9.11E-31; % kg
m = 0.067*m0; % TODO: correct this

a = deltaX = deltaY; % TODO: correct this

t = (hbar^2)/(2*m*a^2); % Hopping Energy

% Mode-to-Slice (0 -> 1)


% Slice-to-Slice (1 -> Nx)


% Slice-to-Mode (Nx+1, Nx+2, collect info at Nx+2)
