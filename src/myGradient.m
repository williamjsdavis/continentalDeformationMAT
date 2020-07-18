function [Fy,Fx] = myGradient(F,h)
%Custom gradient function 
%   William Davis, 16/07/20
%
%   Notes:
%   Solves x and y derivatives by finite differences.
%
%   Inputs:
%   - "F"                       Initial velocity in x-direction, []
%   - "h"                       Initial velocity in y-direction, []
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x direction
Fx = zeros(size(F));
Fx(1,:) = (F(2,:) - F(1,:))/h;
Fx(end,:) = (F(end,:) - F(end-1,:))/h;
Fx(2:end-1,:) = (F(3:end,:) - F(1:end-2,:)) ./ (2*h);

% y direction
Fy = zeros(size(F));
Fy(:,1) = (F(:,2) - F(:,1))/h;
Fy(:,end) = (F(:,end) - F(:,end-1))/h;
Fy(:,2:end-1) = (F(:,3:end) - F(:,1:end-2)) ./ (2*h);
end