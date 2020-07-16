function [ S_new ] = upwind_s( Ux,Uy,S,h,dt,S_bound )
%Upwind scheme to solve for thickness
%   William Davis, 08/12/17
%
%   Notes:
%   Solves the equation s_{t} = -div(su) in relation to
%   England and McKensie 1982 & 1983.
%   Assumes gris spaces dx=dy=h.
%
%   Inputs:
%   - "Ux"                      Velocity in x-direction, []
%   - "Uy"                      Velocity in y-direction, []
%   - "S"                       Crustal thickness, []
%   - "h"                       Spatial grid size
%   - "dt"                      Time-step, []
%   - "S_bound"                 South boundary type: 'const' or 'neu'
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Iteration independent steps
% Parameters
Nx = length(Ux); % Number of nodes on one direction

% Preallocating
UxE = zeros(Nx);
UxW = UxE;
UyN = UxE;
UyS = UxE;

% Cardinal direction nodes
Ux_cent = Ux(2:Nx-1,2:Nx-1); % Center nodes
Ux_east = Ux(3:Nx,2:Nx-1); % East nodes
Ux_west = Ux(1:Nx-2,2:Nx-1); % West nodes
Uy_cent = Uy(2:Nx-1,2:Nx-1); % Center nodes
Uy_nort = Uy(2:Nx-1,3:Nx); % North nodes
Uy_sout = Uy(2:Nx-1,1:Nx-2); % South nodes

% Half step approximation
UxE_cent = 0.5*(Ux_east + Ux_cent); % East velocity
UxW_cent = 0.5*(Ux_west + Ux_cent); % West velocity
UyN_cent = 0.5*(Uy_nort + Uy_cent); % East velocity
UyS_cent = 0.5*(Uy_sout + Uy_cent); % West velocity

% Filling central nodes
% UxE(2:Nx-1,2:Nx-1) = UxE_cent; % Needed?
% UxW(2:Nx-1,2:Nx-1) = UxW_cent;
% UyN(2:Nx-1,2:Nx-1) = UyN_cent;
% UyS(2:Nx-1,2:Nx-1) = UyS_cent;

% Cardinal thicknesses
S_cent = S(2:Nx-1,2:Nx-1); % Center nodes
S_east = S(3:Nx,2:Nx-1); % East nodes
S_nort = S(2:Nx-1,3:Nx); % North nodes
S_west = S(1:Nx-2,2:Nx-1); % West nodes
S_sout = S(2:Nx-1,1:Nx-2); % South nodes

% UxE_sig = sign(UxE(2:Nx-1,2:Nx-1));
% UxW_sig = sign(UxW(2:Nx-1,2:Nx-1));
% UyN_sig = sign(UyN(2:Nx-1,2:Nx-1));
% UyS_sig = sign(UyS(2:Nx-1,2:Nx-1));

% Determine signs of velocities for upwind scheme
UxE_sig = sign(UxE_cent);
UxW_sig = sign(UxW_cent);
UyN_sig = sign(UyN_cent);
UyS_sig = sign(UyS_cent);

% Thicknesses are assigned by signs of velocities
SE_cent = 0.5*((1 + UxE_sig).*S_cent + (1 - UxE_sig).*S_east);
SW_cent = 0.5*((1 + UxW_sig).*S_west + (1 - UxW_sig).*S_cent);
SN_cent = 0.5*((1 + UyN_sig).*S_cent + (1 - UyN_sig).*S_nort);
SS_cent = 0.5*((1 + UyS_sig).*S_sout + (1 - UyS_sig).*S_cent);

%% Solving for new thickness
% Preallocating size
S_new = S; % New thickness

% Stability criterion (Courant-Friedrichs-Lewy)
stability = dt/(max(abs(Ux(:)))+max(abs(Uy(:))));
if stability > 1
    dt_suggest = (max(abs(Ux(:)))+max(abs(Uy(:))));
    error(['Timestep too large. Does not satisfy Courant-Friedrichs-Lewy criterion. Suggestion: Change time-step to ',...
        num2str(dt_suggest)])
end
%coeff = 1; % Stability criterion (Courant-Friedrichs-Lewy)
%dt = coeff*(max(abs(Ux(:)))+max(abs(Uy(:)))); % Evaluate timestep

% Center nodes
S_new(2:Nx-1,2:Nx-1) = S_cent - dt/h*(SE_cent.*UxE_cent - SW_cent.*UxW_cent + SN_cent.*UyN_cent - SS_cent.*UyS_cent);

% Boundary nodes (matrix method)
I = speye(Nx);
I1 = diag([1;zeros(Nx-1,1)]);
I2 = I1; I2(1,1) = 0; I2(1,2) = 4/3; I2(1,3) = -1/3;
B = I-I1+I2; B(Nx,Nx) = 0; B(Nx,Nx-1) = 4/3; B(Nx,Nx-2) = -1/3;
if strcmp(S_bound,'const') == 1
    M = full(kron(B-I2,B) + kron(I1,I)); % Constant on South boundary
elseif strcmp(S_bound,'neu') == 1
    M = kron(B,B); % Neumann on South boundary
end

S_new(:) = M*S_new(:); % Forward problem to find all nodes
end