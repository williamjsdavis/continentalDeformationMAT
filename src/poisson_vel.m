function [Ux_new,Uy_new] = poisson_vel(Ux,Uy,Mx,My,S,h,n,Ar,poisson_set)
%Poisson equation solve for velocity
%   William Davis, 27/11/17
%
%   Notes:
%   Solves the Poisson equation u_{xx} + u_{yy} = f(x,y) in relation to
%   England and McKensie 1982 & 1983.
%   Assumes gris spaces dx=dy=h.
%
%   Inputs:
%   - "Ux"                      Velocity in x-direction, []
%   - "Uy"                      Velocity in y-direction, []
%   - "Mx"                      Velocity stencil in x-direction, []
%   - "My"                      Velocity stencil in y-direction, []
%   - "S"                       Crustal thickness, []
%   - "h"                       Spatial grid size
%   - "n"                       Power law rheology
%   - "Ar"                      Argand number
%   - "poisson_set"             Poisson solver settings
%       - "alpha"                   Stability criterion, ~10E-2
%       - "beta"                    Convergence criterion, ~10E-3
%       - "max_steps"               Maximum number of iteration steps
%
%   Problems:
%   - % Neumann condition on Uy West boundary (check line 135)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate
beta_arr = nan(poisson_set.max_steps,2); % Convergence metric
Ux_old = Ux;
Uy_old = Uy;
Ux_p  = Ux;
Uy_p  = Uy;

%% Iteration steps
%%% Step 1
p_step = 1; % Initialize step counter

% Setting up RHS
RHS_x = Ux; % X-direction
RHS_y = Uy; % Y-direction
RHS_x(2:end-1,2:end-1) = 0; % Set as 0
RHS_y(2:end-1,2:end-1) = 0; % Y-direction

%%% Solving inverse problems
[Ux_p,Uy_p] = solveVelocity(Ux_p,Uy_p,Mx,My,RHS_x,RHS_y,Ux_old,Uy_old);

% Combined with previous solutions for velocity
Ux_new = Ux_p; % If only one step
Uy_new = Uy_p; % I.e. first and final step

% Reassigning velocities
Ux_old = Ux; 
Uy_old = Uy;
Ux = Ux_new;
Uy = Uy_new;

[beta_arr(p_step,:),beta_curr] = convergenceMetric(Ux_new,Uy_new,Ux_old,Uy_old);

while beta_curr > poisson_set.beta % Run until solution converges
    
    p_step = p_step + 1; % Step counter
    
    % Strain rate tensor
    [e11_dot,e12_dot,e21_dot,e22_dot] = strain_rate(Ux,Uy,h,1);
    
    %%% Setting up RHS (note: RHS only contains interior points)
    [RHS_x,RHS_y] = calculateRHS(Ux,Uy,S,e11_dot,e12_dot,e21_dot,...
        e22_dot,h,n,Ar);
    
    %%% Solving inverse problems
    [Ux_p,Uy_p] = solveVelocity(Ux_p,Uy_p,Mx,My,RHS_x,RHS_y,Ux_old,Uy_old);

    % Combined with previous solutions for velocity
    Ux_new = poisson_set.alpha.*Ux_p + (1 - poisson_set.alpha).*Ux;
    Uy_new = poisson_set.alpha.*Uy_p + (1 - poisson_set.alpha).*Uy;
    
    % Reassigning velocities
    Ux_old = Ux;
    Uy_old = Uy;
    Ux = Ux_new;
    Uy = Uy_new;
    
    % Convergence metric
    [beta_arr(p_step,:),beta_curr] = convergenceMetric(Ux_new,Uy_new,Ux_old,Uy_old);
    
    % Iteration limit
    if p_step > poisson_set.max_steps % Setting limit on number of steps
        figure % Plot failure of convergence
        hold on
        plot(1:max_steps,log10(max(beta_arr,[],2)))
        plot(xlim,log10(poisson_set.beta.*[1,1]),'r--')
        title(['Convergence metric, \beta = 10^{',num2str(log10(beta)),'}'])
        xlabel('Iteration number')
        ylabel('Convergence parameter, log(\beta)')
        box on,grid on
        error(['Solution did not converge in ',num2str(max_steps),' iteration(s)'])
    end
end
end
function [RHS_x,RHS_y] = calculateRHS(Ux,Uy,S,e11_dot,e12_dot,e21_dot,...
                                      e22_dot,h,n,Ar)
%% Setting up RHS (note: RHS only contains interior points)
% Second invariant of the strain rate tensor
E_dot = sqrt(2)*sqrt(e11_dot.^2 + e22_dot.^2 + ...
    e12_dot.^2 + e11_dot.*e22_dot);

% Gradients
[UyDy,~] = myGradient(Uy,h); % UyDx
[~,UxDx] = myGradient(Ux,h); % UxDy

[UyDyy,UyDyx] = myGradient(UyDy,h);
%[UyDxy,UyDxx] = gradient(UyDx,h);
[UxDxy,UxDxx] = myGradient(UxDx,h);
%[UxDyy,UxDyx] = gradient(UxDy,h);

[E_dotDy,E_dotDx] = myGradient(E_dot,h);
[SDy,SDx] = myGradient(S,h);

% X-direction
RHS_x = - 3*UxDxx - 3*UyDyx + 2*(1-1/n)*E_dot.^(-1).*(e11_dot.*E_dotDx ...
        + e12_dot.*E_dotDy + (UxDx + UyDy).*(E_dotDx + E_dotDy)) ...
        + 2*Ar*E_dot.^(1-1/n).*S.*SDx;
% Y-direction
RHS_y = - 3*UyDyy - 3*UxDxy + 2*(1-1/n)*E_dot.^(-1).*(e21_dot.*E_dotDx ...
        + e22_dot.*E_dotDy + (UxDx + UyDy).*(E_dotDx + E_dotDy)) ...
        + 2*Ar*E_dot.^(1-1/n).*S.*SDy;
end
function [Ux_p,Uy_p] = solveVelocity(Ux_p,Uy_p,Mx,My,RHS_x,RHS_y,...
                                     Ux_old,Uy_old)
%% Solving inverse problem for velocity field
% Solutions of the poisson equation (LU)
Ux_p(:) = Mx.U\(Mx.L\(Mx.P*RHS_x(:))); 
Uy_p(:) = My.U\(My.L\(My.P*RHS_y(:))); 

% Boundaries (all Dirichlet)
Ux_p(:,end) = Ux_old(:,end); % North
Ux_p(:,1) = Ux_old(:,1); % South
Ux_p(end,:) = Ux_old(end,:); % East 
Ux_p(1,:) = Ux_old(1,:); % West 

Uy_p(:,end) = Uy_old(:,end); % North 
Uy_p(:,1) = Uy_old(:,1); % South 
Uy_p(end,:) = Uy_old(end,:); % East 

% West Boundary (Neumann)
Uy_p(1,:) = 1/3*(4*Uy_p(2,:)-Uy_p(3,:)); 
end
function [beta_arr,beta_curr] = convergenceMetric(Ux_new,Uy_new,Ux_old,Uy_old)
% Calculate convergence metric
beta_arr(2) = max(Uy_new(:)-Uy_old(:))/max(Uy_old(:)); % Y velocity
beta_arr(1) = max(Ux_new(:)-Ux_old(:))/max(Ux_old(:)); % X velocity
beta_curr = max(beta_arr); % Maximum of either convergence metric
end