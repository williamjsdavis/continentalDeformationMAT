function [ out ] = myGradient( in )
%Custom gradient function 
%   William Davis, 16/07/20
%
%   Notes:
%   Solves the Poisson equation u_{xx} + u_{yy} = f(x,y) and the equation 
%   s_{t} = -div(su) in relation to England and McKensie 1982 & 1983.
%   Assumes gris spaces dx=dy=h.
%
%   Inputs:
%   - "F"                       Initial velocity in x-direction, []
%   - "Uy"                      Initial velocity in y-direction, []
%   - "S"                       Initial crustal thickness, []
%   - "h"                       Spatial grid size, []
%   - "n"                       Power law rheology, []
%   - "Ar"                      Argand number, []
%   - "dt"                      Time-step, []
%   - "S_bound"                 South boundary type: 'const' or 'neu'
%   - "nt"                      Number of time-steps
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end