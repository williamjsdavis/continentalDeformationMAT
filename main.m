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

%mainFunctional()

mainOO()

%mainAnimation()

%% Functions
function mainFunctional()
% Functional example

% Settings
[L,u0,g,pc,pm,s0,n,Ar,Nx,dt,S_bound,poisson_set] = simulation_settings();

% Setup grids
[Ux,Uy,S,s0,h,X,Y,x,y] = setup_grid(Nx,L,u0,s0,dt);

% Solving
nt = 100; % Number of timesteps
[Ux_new,Uy_new,S_new] = time_solve(Ux,Uy,S,h,n,Ar,dt,S_bound,nt,poisson_set);
disp('Complete!')

% Plotting
plot_cont;
plot_3D;
end
function mainOO()
% Object oriented example

% Instantiate field object
thinViscousSheet = TVSfield();

% Setup grids and stencils
thinViscousSheet.setupGrids();

% Solve
nStep = 100; % Number of timesteps
thinViscousSheet.timeSolve(nStep);

% Plot
%thinViscousSheet.plot6();

figure('Position',[399,259,983,517])
thinViscousSheet.plot3D('default');
end
function mainAnimation()
% Animation example

% Instantiate field object and setup grids
thinViscousSheet = TVSfield();
thinViscousSheet.simSettings.Ar = 1;
thinViscousSheet.setupGrids();

% Solve and plot
nStep = 100; % Number of timesteps
i = 0;

figure('Position',[399,259,983,517])
while i < nStep
    thinViscousSheet.timeSolve(1);
    thinViscousSheet.plot3D('interp');
    drawnow
    i = i + 1;
end

end