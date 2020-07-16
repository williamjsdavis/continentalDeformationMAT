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

%% Upwind setup
[ S_cent,SE_cent,UxE_cent,SW_cent,UxW_cent,SN_cent,UyN_cent,SS_cent,...
    UyS_cent,B,I,I1,I2 ] = setup_upwind( Ux,Uy,S );

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
S_new(2:end-1,2:end-1) = S_cent - dt/h*(SE_cent.*UxE_cent - SW_cent.*UxW_cent + SN_cent.*UyN_cent - SS_cent.*UyS_cent);

if strcmp(S_bound,'const') == 1
    M = full(kron(B-I2,B) + kron(I1,I)); % Constant on South boundary
elseif strcmp(S_bound,'neu') == 1
    M = kron(B,B); % Neumann on South boundary
end

S_new(:) = M*S_new(:); % Forward problem to find all nodes
end