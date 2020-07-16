function [  ] = plot_elipsiod( tensor,X,Y,L,scale )
%Plots elipsiods for tensors
%   William Davis, 15/12/17
%
%   Notes: Designed for 2x2 matrix tensor
%   
%
%   Inputs:
%   - "tensor"                  Tensor for plotting (2x2 matrix)
%   - "Uy"                      Velocity in y-direction, []
%   - "S"                       Crustal thickness, []
%   - "h"                       Spatial grid size
%   - "dt"                      Time-step, []
%   - "S_bound"                 South boundary type: 'const' or 'neu'
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = length(X);

%% Eigenvalue decomposition
lambda1 = 0.5*(tensor{1,1} + tensor{2,2} + sqrt((tensor{1,1} + tensor{2,2}).^2 + 4*tensor{1,2}.^2));
lambda2 = 0.5*(tensor{1,1} + tensor{2,2} - sqrt((tensor{1,1} + tensor{2,2}).^2 + 4*tensor{1,2}.^2));
ang = 0.5*atan2(2*tensor{1,2},tensor{1,1}-tensor{2,2});

t = 0:0.1:2*pi;

hold on
for i = 1:(Nx/16):Nx
    for j = 1:(Nx/16):Nx
        x = X(i,j) + scale*lambda1(i,j)*cos(t)*cos(ang(i,j)) - scale*lambda2(i,j)*sin(t)*sin(ang(i,j));
        y = Y(i,j) + scale*lambda1(i,j)*cos(t)*sin(ang(i,j)) + scale*lambda2(i,j)*sin(t)*cos(ang(i,j));
        fill(L*x,L*y,'k')
    end
end

end

