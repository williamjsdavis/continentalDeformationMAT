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

mainExample()

%mainAnimation()

%% Functions
function mainExample()
% Object oriented example

% Instantiate field object
thinViscousSheet = TVSfield();

% Setup grids and stencils
thinViscousSheet.setupGrids();

% Solve
nTimeSteps = 100;
thinViscousSheet.timeSolve(nTimeSteps);

% Plots
thinViscousSheet.plot6();

figure('Position',[399,259,983,517])
thinViscousSheet.plot3D('default');
end
function mainAnimation()
% Animation example

% Instantiate field object and setup grids
thinViscousSheet = TVSfield();
thinViscousSheet.simSettings.Ar = 10;
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