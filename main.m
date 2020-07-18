%% Continental Deformation
%   William Davis, 16/07/20
%
%   Notes:
%   Solves equations for the thin viscous sheet approximation, England and
%   McKensie 1982 & 1983. 
%   Updated from previous 25/11/17 work.
%
%   Input files:
%   - 
%
%   Problems:
%   -
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
addpath('src/')

% Settings
[L,u0,g,pc,pm,s0,n,Ar,Nx,dt,nt,S_bound,poisson_set] = simulation_settings();

% Setup grids
[Ux,Uy,S,s0,h,X,Y,x,y] = setup_grid(Nx,L,u0,s0,dt);

% Processing
[Ux_new,Uy_new,S_new] = time_solve(Ux,Uy,S,h,n,Ar,dt,S_bound,nt,poisson_set);
disp('Complete!')

% Plotting
plot_cont;

