function [ Mx,My ] = setup_poisson( Nx,h )
%Poisson equation solve for velocity
%   William Davis, 16/07/20
%
%   Notes:
%   Setup stencils for finite difference solution of Poisson equation 
%           
%               u_{xx} + u_{yy} = f(x,y) 
%
%   in relation to England and McKensie 1982 & 1983.
%   Assumes gris spaces dx=dy=h.
%
%   Inputs:
%   - "Ux"                      Velocity in x-direction, []
%   - "Uy"                      Velocity in y-direction, []
%   - "S"                       Crustal thickness, []
%   - "h"                       Spatial grid size
%   - "n"                       Power law rheology
%   - "Ar"                      Argand number
%   - "alpha"                   Stability criterion, ~10E-2
%   - "beta"                    Convergence criterion, ~10E-3
%   - "max_steps"               Maximum number of iteration steps
%
%   Problems:
%   - % Neumann condition on Uy West boundary (check line 135)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Form sparse matrices
% Identity matrices
I = ones(Nx,1); % Identity array
I_NS = I; % Block array for North and South stencil dependency
I_NS(Nx) = 0; % Remove East Dirichlet boundary dependency
I_EW = I_NS; % Block array for East and West stencil dependency
I_EW(1) = 0; % Removing South boundary (North was already removed)
I_NSB = I;  % Block array for North and South boundaries
I_NSB(2:Nx-1) = 0; % Removing interior nodes

% Block matrices
B_EW = diag(I_EW); % Block matrix for East and West stencil dependency
B_NSB = diag(I_NSB); % Block matrix for North and South boundaries
B_NS = spdiags([I,I],[-1,1],Nx,Nx); % Block matrix for North and South stencil dependency
B_NS(1,2) = 0; % Remove South boundary dependency
B_NS(Nx,Nx-1) = 0; % Remove North boundary dependency

% Stencil matrices
Sx_NS = diag(I_EW); % North and South stencil dependency (makes use of EW matrix already formed)
Sx_NSB = h^2*speye(Nx); % North and South boundary stencil dependency
Sx_EW = spdiags([I,-4*I,I],[-1,0,1],Nx,Nx); % East and West stencil dependency
Sx_EW(1,1) = h^2; % Add for West Dirichlet boundary (*)
Sx_EW(1,2) = 0; % Remove for West Dirichlet boundary
Sx_EW(Nx,Nx) = h^2; % Add for East Dirichlet boundary (*)
Sx_EW(Nx,Nx-1) = 0; % Remove for East Dirichlet boundary

Sy_NS = diag(I_NS); % North and South stencil dependency
Sy_NSB = h^2*speye(Nx); % North and South boundary stencil dependency (*)
Sy_EW = spdiags([I,-4*I,I],[-1,0,1],Nx,Nx); % East and West stencil dependency
Sy_EW(1,1) = -13/4; % Add West Neumann boundary
Sy_EW(1,3) = 1/4; % Add West Neumann boundary
Sy_EW(Nx,Nx) = h^2; % Add for East Dirichlet boundary (*)
Sy_EW(Nx,Nx-1) = 0; % Remove for East Dirichlet boundary

% (*) Will be replaced later so can be any number, but set to h^2 for ease 
% of identification.

% Kronecker tensor products
Mx_NS = kron(B_NS,Sx_NS); % Matrix for North and South stencil dependency
Mx_EW = kron(B_EW,Sx_EW); % Matrix for East and West stencil dependency
Mx_NSB = kron(B_NSB,Sx_NSB); % Matrix for North and South boundaries

My_NS = kron(B_NS,Sy_NS); % Matrix for North and South stencil dependency
My_EW = kron(B_EW,Sy_EW); % Matrix for East and West stencil dependency
My_NSB = kron(B_NSB,Sy_NSB); % Matrix for North and South boundaries

% Full sparse matrices
Mxx = (Mx_NS+Mx_EW+Mx_NSB)/h^2; % Add all to form sparse matrix for Ux (spA_x)
Myy = (My_NS+My_EW+My_NSB)/h^2; % Add all to form sparse matrix for Uy (spA_y)

% LU decomposition
[Mx.L,Mx.U,Mx.P] = lu(Mxx);
[My.L,My.U,My.P] = lu(Myy);
end
