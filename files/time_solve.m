function [ Ux,Uy,S ] = time_solve( Ux0,Uy0,S0,h,n,Ar,dt,S_bound,nt )
%Solves the solution in space and time
%   William Davis, 12/12/17
%
%   Notes:
%   Solves the Poisson equation u_{xx} + u_{yy} = f(x,y) and the equation 
%   s_{t} = -div(su) in relation to England and McKensie 1982 & 1983.
%   Assumes gris spaces dx=dy=h.
%
%   Inputs:
%   - "Ux"                      Initial velocity in x-direction, []
%   - "Uy"                      Initial velocity in y-direction, []
%   - "S"                       Initial crustal thickness, []
%   - "h"                       Spatial grid size, []
%   - "n"                       Power law rheology, []
%   - "Ar"                      Argand number, []
%   - "dt"                      Time-step, []
%   - "S_bound"                 South boundary type: 'const' or 'neu'
%   - "nt"                      Number of time-steps
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup
% Parameters
max_steps = 5000; % Maximum iterations
alpha = 2.5E-3; % Stability parameter
beta = 1E-3; % Convergence criterion

% Preallocating sizes
Ux = Ux0; Ux_old = Ux0; Ux_new = Ux0;
Uy = Uy0; Uy_old = Uy0; Uy_new = Uy0;
S = S0; S_old = S0; S_new = S0;

% Time-stepping

for i = 1:nt
    % Poisson solving (for velocity)
    [Ux_new,Uy_new,~] = poisson_vel(Ux,Uy,S,h,n,Ar,alpha,beta,max_steps); % New velocities

    % Upwind solving (for thickness)
    S_new = upwind_s(Ux_new,Uy_new,S,h,dt,S_bound); % New thicknesses

    % Change variables
    Ux_old = Ux; Uy_old = Uy; S_old = S;
    Ux = Ux_new; Uy = Uy_new; S = S_new;
    
    % Print message
    disp(['Progress: ',sprintf('%5.2f',100*i/nt),'%'])
end

end

