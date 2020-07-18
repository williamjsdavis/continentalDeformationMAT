function [Ux,Uy,S,s0,h,X,Y,x,y] = setup_grid(Nx,L,u0,s0,dt)
%Set up grids
%   William Davis, 17/07/20
%
%   Notes:
%   Setup grids for fields. Assumes gris spaces dx=dy=h.
%
%   Inputs:
%   - "Nx"                      Number of grid points, []
%   - "L"                       Length scale, km
%   - "u0"                      Initial velocity, mm/yr
%   - "s0"                      Initial crustal thickness, km
%   - "dt"                      Non-dimensional time, []
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Non-dimensionalise initial setup
h = Nx/(Nx-1); % Grid step size, []
s0 = s0/L; % Crustal length scale, []
dt_dim = dt*u0/L; % Dimensional time, Ma

% Grid points x/y including boundaries, []
x = linspace(0,Nx,Nx); 
y = linspace(0,Nx,Nx);

% Position meshgrid
% Transpose so that X(i,j),Y(i,j) are coordinates of (i,j) point
[X,Y] = meshgrid(x,y); 
X = X'; 
Y = Y'; 

% Crustal thickness meshgrid
S = s0*ones(Nx); 

% X/Y-direction velocity meshgrid
Ux = zeros(Nx); 
Uy = zeros(Nx);

% Indenter function, velocities
x_f = linspace(0,Nx,Nx); 
x_sec1 = floor(Nx/4); 
x_sec2 = floor(Nx/2); 
Fx = zeros(1,Nx); 
Fx(1:x_sec1) = 1; 
Fx(x_sec1+1:x_sec2) = cos(pi/2*(4*x_f(x_sec1+1:x_sec2)./Nx-1)).^2; 

% Initial conditions (on South boundary)
Uy(:,1) = Fx;
end
