function [ e11_dot,e12_dot,e21_dot,e22_dot ] = strain_rate( Ux,Uy,h,scale )
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
%   - "scale"                   Scale for plotting, []
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Derivatives of velocity
[UxDy,UxDx] = gradient(Ux,h);
[UyDy,UyDx] = gradient(Uy,h);
% U1X1 = gradient(Ux',h)';   % i=1, j=1, du/dx component
% U1X2 = gradient(Ux,h);     % i=1, j=2, du/dy component
% U2X1 = gradient(Uy',h)';   % i=2, j=1, dv/dx component
% U2X2 = gradient(Uy,h);     % i=2, j=2, dv/dy component

% Strain rate tensor
% e11_dot = 0.5*(U1X1+U1X1)*scale;
% e12_dot = 0.5*(U1X2+U2X1)*scale;
% e21_dot = 0.5*(U2X1+U1X2)*scale;
% e22_dot = 0.5*(U2X2+U2X2)*scale;
e11_dot = UxDx*scale;
[e12_dot,e21_dot] = deal(0.5*(UxDy+UyDx)*scale);
e22_dot = UyDy*scale;

end

