function [ eij_dot ] = strain_rate( Ux,Uy,h,scale )
%Calculates strain rate tensor
%   William Davis, 14/12/17
%
%   Notes:
%   
%
%   Inputs:
%   - "Ux"                      Velocity in x-direction, []
%   - "Uy"                      Velocity in y-direction, []
%   - "h"                       Spatial grid size, []
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocating
UiXj = cell(2); 
eij_dot = cell(2);

% Derivatives of velocity
UiXj{1,1} = gradient(Ux',h)';   % i=1, j=1, du/dx component
UiXj{1,2} = gradient(Ux,h);     % i=1, j=2, du/dy component
UiXj{2,1} = gradient(Uy',h)';   % i=2, j=1, dv/dx component
UiXj{2,2} = gradient(Uy,h);     % i=2, j=2, dv/dy component

% Strain rate tensor
for i = 1:2
    for j = 1:2
        eij_dot{i,j} = 0.5*(UiXj{i,j}+UiXj{j,i})*scale;
    end
end

end

