%% Continental Deformation
%   William Davis, 25/11/17
%
%   Notes:
%   Solves equations for the thin viscous sheet approximation, England and
%   McKensie 1982 & 1983. 
%
%   Input files:
%   - 
%
%   Problems:
%   -
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
%% Inputs
L = 100; % Length scale, km
u0 = 50; % Initial velocity, mm/yr
g = 9.81; % Gravitational acceleration, m/s^2
pc = 2.95; % Crustal density, Mg/m^3
pm = 3.30; % Mantle density, Mg/m^3
s0 = 35; % Initial crustal thickness, km
n = 1; % Power law rheology
Ar = 1; % Argand number

Nx = 16; % Number of nodes in x-direction (Ny is the same)

% Poisson solving
max_steps = 5000; % Maximum iterations
alpha = 2.5E-3; % Stability parameter
beta = 1E-3; % Convergence criterion
% Time-solving
S_bound = 'const'; % South boundary condition: 'const' or 'neu'
dt = 0.1; % Non-dimensional time
nt = 10; % Number of time-steps

% Plotting
save_figure = 0; % Save figure?
filename = strcat('Ar',num2str(Ar),'n',num2str(n),'t',num2str(nt)); % File to save as
cbar_range1 = 1E-15*[-2,20]; % Minimum and maximum isotropic strain rates (chosen) [s^-1]
cbar_range2 = [20,65]; % Minimum and maximum thicknesses (chosen) [km]
x_section = [1,4,6,12]; % Where to plot cross sections (change this with grid size)
%x_section = [2,8,16,24]; % Where to plot cross sections (change this with grid size)

% Adding path
addpath('files')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup

% Non-dimensionalise initial setup
h = Nx/(Nx-1); % Grid step size, []
s0 = s0/L; % Crustal length scale, []
dt_dim = dt*u0/L; % Dimensional time, Ma

x = linspace(0,Nx,Nx); % Grid points x including boundaries, []
y = linspace(0,Nx,Nx); % Grid points y including boundaries, []

[X,Y] = meshgrid(x,y); % Position meshgrid
X = X'; % Transpose so that X(i,j),Y(i,j) are
Y = Y'; % coordinates of (i,j) point

S = s0*ones(Nx); % Crustal thickness meshgrid
Ux = zeros(Nx); % X-direction velocity meshgrid
Uy = zeros(Nx); % Y-direction velocity meshgrid

% Indenter function, velocities
x_f = linspace(0,Nx,Nx); % x Indices for F function
x_sec1 = floor(Nx/4); % End node for section 1
x_sec2 = floor(Nx/2); % End node for section 2
Fx = zeros(1,Nx); % Background velocity, []
Fx(1:x_sec1) = 1; % First section velocity, []
Fx(x_sec1+1:x_sec2) = cos(pi/2*(4*x_f(x_sec1+1:x_sec2)./Nx-1)).^2; % Second section velocity, []

% Initial conditions (boundary conditions)
Uy(:,1) = Fx; % Add to y-direction velocity meshgrid


%% Processing

[Ux_new,Uy_new,S_new] = time_solve(Ux,Uy,S,h,n,Ar,dt,S_bound,nt);
disp('Complete!')

%% Plotting

plot_cont;

