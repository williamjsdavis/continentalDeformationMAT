function [Ux,Uy,S] = time_solve(Ux0,Uy0,S0,h,n,Ar,dt,S_bound,nt,poisson_set)
%Solves the solution in space and time
%   William Davis, 12/12/17
%
%   Notes:
%   Solves the Poisson equation u_{xx} + u_{yy} = f(x,y) and the equation 
%   s_{t} = -div(su) in relation to England and McKensie 1982 & 1983.
%   Assumes gris spaces dx=dy=h.
%
%   Inputs:
%   - "Ux0"                     Initial velocity in x-direction, []
%   - "Uy0"                     Initial velocity in y-direction, []
%   - "S0"                      Initial crustal thickness, []
%   - "h"                       Spatial grid size, []
%   - "n"                       Power law rheology, []
%   - "Ar"                      Argand number, []
%   - "dt"                      Time-step, []
%   - "S_bound"                 South boundary type: 'const' or 'neu'
%   - "nt"                      Number of time-steps
%   - "poisson_set"             Poisson solver settings
%       - "alpha"                   Stability criterion, ~10E-2
%       - "beta"                    Convergence criterion, ~10E-3
%       - "max_steps"               Maximum number of iteration steps
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocating sizes
Ux = Ux0; 
Uy = Uy0; 
S = S0; 

% Set up Poisson stencils
Nx = length(Ux);
[Mx,My] = setup_poisson(Nx,h);

% Time-stepping
for i = 1:nt
    % Solve for velocity
    [Ux_new,Uy_new] = poisson_vel(Ux,Uy,Mx,My,S,h,n,Ar,poisson_set); % New velocities

    % Solve for thickness
    S_new = upwind_s(Ux_new,Uy_new,S,h,dt,S_bound); % New thicknesses

    % Update
    Ux = Ux_new; Uy = Uy_new; S = S_new;
    
    % Print message
    disp(['Progress: ',sprintf('%5.2f',100*i/nt),'%'])
end

end

