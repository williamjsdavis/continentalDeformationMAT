function [ Ux_new,Uy_new,beta_arr ] = poisson_velint( steps,Nx,h,Ux,Uy,S,n,Ar,alpha )
%Poisson equation solve for velocity
%   William Davis, 27/11/17
%
%   Notes:
%   Solves the Poisson equation u_{xx} + u_{yy} = f(x,y) in relation to
%   England and McKensie 1982 & 1983.
%   Interior points method.
%
%   Different settingsstep
%   x or y      - Partial derivative in one direction
%   xx or yy    - Double partial derivative in one direction
%   xy          - Double partial derivative in two directions
%
%   Assumes gris spaces dx=dy=h.
%
%   Inputs:
%   - "steps"                   Number of iteration steps
%   - "Nx"                      Number of interior nodes in one dimension
%   - "h"                       Spatial grid size
%   - "Ux"                      Velocity in x-direction, []
%   - "Uy"                      Velocity in y-direction, []
%   - "S"                       Crustal thickness, []
%   - "n"                       Power law rheology
%   - "Ar"                      Argand number
%   - "alpha"                   Stability criterion, ~10E-2
%
%   Problems:
%   - Maybe y axis is inverted?
%   - Fix boundary conditions?
%   - Look at units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Iteration independent steps
% Preallocate convergence array
beta_arr = zeros(steps,2); % Array to put convergence metric for each step

%%% Form sparse matrices
spI = speye(Nx); % Sparse identity matrix
e = ones(Nx,1); % Column of ones
spS = spdiags([e,e],[-1,1],Nx,Nx); % Forming block matrix for edge part

% X-direction
Uxint_p = zeros(Nx); % Preallocating dimensions
spT_x = spdiags([e,-4*e,e],[-1,0,1],Nx,Nx); % Forming block matrix for central part
spA_x = (kron(spI,spT_x) + kron(spS,spI))/h^2; % Use Kronecker tensor product to form matrix

% Y-direction
Uyint_p = zeros(Nx); % Preallocating dimensions
eN_y = e; % Column of ones
eN_y(1:2) = 2; % Neumann boundary condition
spT_y = spdiags([e,-4*e,eN_y],[-1,0,1],Nx,Nx); % Forming block matrix for central part
spA_y = (kron(spI,spT_y) + kron(spS,spI))/h^2; % Use Kronecker tensor product to form matrix

%%% Other
% Interior points
Iint = 2:Nx+1; % Indices of interior points in x
Jint = 2:Nx+1; % Indices of interior points in y
Uxint = Ux(Iint,Jint); % Interior velocity in x-direction
Uyint = Uy(Iint,Jint); % Interior velocity in y-direction
Sint = S(Iint,Jint); % Interior crustal thickness

% Exterior points
Ux_bound = Ux; % X-direction
Ux_bound(2:Nx+1,2:Nx+1) = 0; % Setting boundaries
Uy_bound = Uy; % Y-direction
Uy_bound(2:Nx+1,2:Nx+1) = 0; % Setting boundaries

UiXj = cell(2); % Preallocating and preassigning
eij_dot = cell(2);
%{
% % Derivatives of velocity
% UiXj = cell(2); % Preallocating
% UiXj{1,1} = gradient(Ux',h)';   % i=1, j=1, du/dx component
% UiXj{1,2} = gradient(Ux,h); % i=1, j=2, du/dy component
% UiXj{2,1} = gradient(Uy',h)';   % i=2, j=1, dv/dx component
% UiXj{2,2} = gradient(Uy,h); % i=2, j=2, dv/dy component
% for k = 1:4
%     UiXj{k} = UiXj{k}(2:Nx+1,2:Nx+1); % Interior points
% end
% 
% % Strain rate tensor
% eij_dot = cell(2); % Preallocating
% for i = 1:2
%     for j = 1:2
%         eij_dot{i,j} = 0.5*(UiXj{i,j}+UiXj{j,i});
%     end
% end
% 
% % Second invariant of the strain rate tensor
% E_dot = sqrt(2)*sqrt(eij_dot{1,1}.^2 + eij_dot{2,2}.^2 + eij_dot{1,2}.^2 + eij_dot{1,1}.*eij_dot{2,2});
%}

%% Iteration steps
for p_step = 1:steps
    %%% Parameters from velocity
    % Derivatives of velocity
    UiXj{1,1} = gradient((Ux_bound+padarray(Uxint,[1,1],0))',h)'; % i=1, j=1, du/dx component
    UiXj{1,2} = gradient(Ux_bound+padarray(Uxint,[1,1],0),h);     % i=1, j=2, du/dy component
    UiXj{2,1} = gradient((Uy_bound+padarray(Uyint,[1,1],0))',h)'; % i=2, j=1, dv/dx component
    UiXj{2,2} = gradient(Uy_bound+padarray(Uyint,[1,1],0),h);     % i=2, j=2, dv/dy component
    for k = 1:4
        UiXj{k} = UiXj{k}(2:Nx+1,2:Nx+1); % Interior points
    end
    
    % Strain rate tensor
    for i = 1:2
        for j = 1:2
            eij_dot{i,j} = 0.5*(UiXj{i,j}+UiXj{j,i});
        end
    end
    
    % Second invariant of the strain rate tensor
    E_dot = sqrt(2)*sqrt(eij_dot{1,1}.^2 + eij_dot{2,2}.^2 + eij_dot{1,2}.^2 + eij_dot{1,1}.*eij_dot{2,2});
    
    %%% Setting up RHS (note: RHS only contains interior points)
    if p_step == 1
        RHS_x = zeros(Nx); % X-direction
        RHS_y = zeros(Nx); % Y-direction
    else
        % X-direction
        RHS_x = - del__g(Uxint,h,'xx') - del__g(Uyint,h,'xy') + 2*(1-1/n)*E_dot.^(-1).*(eij_dot{1,1}.*del__g(E_dot,h,'x') ...
            + eij_dot{1,2}.*del__g(E_dot,h,'y')) + Ar*E_dot.^(1-1/n).*Sint.*del__g(Sint,h,'x');
        % Y-direction
        RHS_y = - del__g(Uyint,h,'yy') - del__g(Uxint,h,'xy') + 2*(1-1/n)*E_dot.^(-1).*(eij_dot{2,1}.*del__g(E_dot,h,'x') ...
            + eij_dot{2,2}.*del__g(E_dot,h,'y')) + Ar*E_dot.^(1-1/n).*Sint.*del__g(Sint,h,'y');
    end
    
    % Dirichlet boundary conditions
    % X-direction (should all be 0)
    RHS_x(:,Nx) = RHS_x(:,Nx) - Ux(Iint,Nx+2)/h^2; % North boundary
    RHS_x(Nx,:) = RHS_x(Nx,:) - Ux(Nx+2,Jint)/h^2; % East boundary
    RHS_x(:,1) = RHS_x(:,1) - Ux(Iint,1)/h^2; % South boundary
    RHS_x(1,:) = RHS_x(1,:) - Ux(1,Jint)/h^2; % West boundary
    
    % Y-direction
    RHS_y(:,Nx) = RHS_y(:,Nx) - Uy(Iint,Nx+2)/h^2; % North boundary
    RHS_y(Nx,:) = RHS_y(Nx,:) - Uy(Nx+2,Jint)/h^2; % East boundary
    RHS_y(:,1) = RHS_y(:,1) - Uy(Iint,1)/h^2; % South boundary
    %RHS_y(1,:) = RHS_y(1,:) - Uy(1,Jint)/h^2; % West boundary (not altered -> Neumann)
    
    %%% Solving inverse problems
    Uxint_p(:) = spA_x\RHS_x(:); % Solutions of the poisson equation
    Uyint_p(:) = spA_y\RHS_y(:); % 
    
    % Combined with previous solutions for velocity
    if (p_step == steps) && (steps == 1)
        beta_arr(p_step,1) = max(Uxint_p(:)-Uxint(:))/max(Uxint(:)); % Convergence metric, X velocity
        beta_arr(p_step,2) = max(Uyint_p(:)-Uyint(:))/max(Uyint(:)); % Y velocity
        Uxint_new = Uxint_p; % If only one step
        Uyint_new = Uyint_p; % I.e. first and final step
    elseif p_step == steps
        Uxint_new = alpha.*Uxint_p + (1 - alpha).*Uxint; % Stability condition
        Uyint_new = alpha.*Uyint_p + (1 - alpha).*Uyint; % Final step
        beta_arr(p_step,1) = max(Uxint_new(:)-Uxint(:))/max(Uxint(:)); % Convergence metric, X velocity
        beta_arr(p_step,2) = max(Uyint_new(:)-Uyint(:))/max(Uyint(:)); % Y velocity
    elseif p_step == 1
        beta_arr(p_step,1) = max(Uxint_p(:)-Uxint(:))/max(Uxint(:)); % Convergence metric, X velocity
        beta_arr(p_step,2) = max(Uyint_p(:)-Uyint(:))/max(Uyint(:)); % Convergence metric
        Uxint = Uxint_p; % Initial step
        Uyint = Uyint_p;
    else
        Uxint_new = alpha.*Uxint_p + (1 - alpha).*Uxint; % Stability condition
        Uyint_new = alpha.*Uyint_p + (1 - alpha).*Uyint; % Intermediate step
        beta_arr(p_step,1) = max(Uxint_new(:)-Uxint(:))/max(Uxint(:)); % Convergence metric, X velocity
        beta_arr(p_step,2) = max(Uyint_new(:)-Uyint(:))/max(Uyint(:)); % Convergence metric
        Uxint = Uxint_new; % Reassigning velocities
        Uyint = Uyint_new;
    end
    
    %figure(1),surf(50*Uyint),drawnow
end

% Adding exterior points
% X-direction
Ux_new = padarray(Uxint_new,[1,1],0) + Ux_bound; % Dirichlet

% Y-direction
Uy_new = padarray(Uyint_new,[1,1],0) + Uy_bound; % Dirichlet
Uy_new(1,2:Nx+1) = Uy_new(3,2:Nx+1); % Neumann











end

