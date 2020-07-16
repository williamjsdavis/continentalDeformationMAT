function [ a_new ] = del__g( a,h,type )
%Finite difference derivatives
%   William Davis, 26/11/17
%
%   Notes:
%   Discretised (finite differences) partial derivative for 3 options:
%   x or y      - Partial derivative in one direction
%   xx or yy    - Double partial derivative in one direction
%   xy          - Double partial derivative in two directions
%   
%   Assumes gris spaces dx=dy=Delta.
%
%   Inputs:
%   - "a"                       Gridded input function
%   - "h"                       Grid spacing in x and y direction
%   - "type"                    Type of derivative to calculate
%
%   Problems:
%   - Maybe y axis is inverted?
%   - Fix boundary conditions?
%   - Look at units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set grid spacing
dx = h; % x-dir grid spacing
dy = h; % x-dir grid spacing
Nx = length(a); % Number of nodes on one direction
a_new = a; % Preallocating size

if (strcmp(type,'x') == 1) ||  (strcmp(type,'y') == 1) || strcmp(type,'xy')
    I = ones(Nx,1); % Identity array
    D = spdiags([-I,I],[-1,1],Nx,Nx);
    D(1,1) = -3;
    D(1,2) = 4;
    D(1,3) = -1;
    D(Nx,Nx) = 3;
    D(Nx,Nx-1) = -4;
    D(Nx,Nx-2) = 1;
    if strcmp(type,'x') == 1
        M = 1/(2*dx)*kron(diag(I),D);
    elseif strcmp(type,'y') == 1
        M = 1/(2*dy)*kron(D,diag(I));
    elseif strcmp(type,'xy') == 1
        M = 1/(4*dx*dy)*kron(D,D);
    end
elseif (strcmp(type,'xx') == 1) ||  (strcmp(type,'yy') == 1)
    I = ones(Nx,1); % Identity array
    D2 = spdiags([I,-2*I,I],[-1,0,1],Nx,Nx);
    D2(1,1) = 1;
    D2(1,2) = -2;
    D2(1,3) = 1;
    D2(Nx,Nx) = 1;
    D2(Nx,Nx-1) = -2;
    D2(Nx,Nx-2) = 1;
    if strcmp(type,'xx') == 1
        M = 1/(dx*dx)*kron(diag(I),D2);
    elseif strcmp(type,'yy') == 1
        M = 1/(dy*dy)*kron(D2,diag(I));
    end
end

% Forward problem
a_new(:) = M*a(:);

