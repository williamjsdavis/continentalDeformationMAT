function [L,u0,g,pc,pm,s0,n,Ar,Nx,dt,nt,S_bound,poisson_set] = simulation_settings()
%Simulation settings
%   William Davis, 17/07/20
%
%   Notes:
%   Setup grids for fields. Assumes gris spaces dx=dy=h.
%
%   Inputs:
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Scales
L = 100; % Length scale, km
u0 = 50; % Initial velocity, mm/yr
g = 9.81; % Gravitational acceleration, m/s^2
pc = 2.95; % Crustal density, Mg/m^3
pm = 3.30; % Mantle density, Mg/m^3
s0 = 35; % Initial crustal thickness, km
n = 1; % Power law rheology
Ar = 1; % Argand number

% Grid nodes in x/y-direction
Nx = 16; 

% Poisson solver settings
poisson_set.max_steps = 5000; % Maximum iterations
poisson_set.alpha = 2.5E-3; % Stability parameter
poisson_set.beta = 1E-3; % Convergence criterion

% Time solver settings
S_bound = 'const'; % South boundary condition: 'const' or 'neu'
dt = 0.1; % Non-dimensional time
nt = 10; % Number of time-steps
end
