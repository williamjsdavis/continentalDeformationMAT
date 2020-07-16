function [ Ux_new,Uy_new,beta_arr ] = poisson_vel( Ux,Uy,Mx,My,S,h,n,Ar,alpha,beta,max_steps )
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
%   - "alpha"                   Stability criterion, ~10E-2
%   - "beta"                    Convergence criterion, ~10E-3
%   - "max_steps"               Maximum number of iteration steps
%
%   Problems:
%   - % Neumann condition on Uy West boundary (check line 135)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate convergence metrics
beta_arr = nan(max_steps,2); % Array to put convergence metric for each step
beta_curr = 1; % Current convergence value

Ux_old = Ux;
Uy_old = Uy;
Ux_p  = Ux;
Uy_p  = Uy;
Ux_new = Ux;
Uy_new = Uy;

%% Iteration steps
p_step = 0; % Initialize step counter
while beta_curr > beta % Run until solution converges
    
    p_step = p_step + 1; % Step counter
    
    %%% Parameters from velocity
    % Strain rate tensor
    eij_dot = strain_rate(Ux,Uy,h,1);
    
    % Second invariant of the strain rate tensor
    E_dot = sqrt(2)*sqrt(eij_dot{1,1}.^2 + eij_dot{2,2}.^2 + eij_dot{1,2}.^2 + eij_dot{1,1}.*eij_dot{2,2});
    
    %%% Setting up RHS (note: RHS only contains interior points)
    if p_step == 1
        RHS_x = Ux; % X-direction
        RHS_y = Uy; % Y-direction
        RHS_x(2:end-1,2:end-1) = 0; % Set as 0
        RHS_y(2:end-1,2:end-1) = 0; % Y-direction
    else
        % X-direction
        RHS_x = - 3*gradient(gradient(Ux',h),h)' - 3*gradient(gradient(Uy',h)',h) + 2*(1-1/n)*E_dot.^(-1).*...
            (eij_dot{1,1}.*gradient(E_dot',h)' + eij_dot{1,2}.*gradient(E_dot,h) + ...
            (gradient(Ux',h)' + gradient(Uy,h)).*(gradient(E_dot',h)' + gradient(E_dot,h))) + ...
            Ar*E_dot.^(1-1/n).*S.*gradient(S',h)';
        % Y-direction
        RHS_y = - 3*gradient(gradient(Uy,h),h) - 3*gradient(gradient(Ux',h)',h) + 2*(1-1/n)*E_dot.^(-1).*...
            (eij_dot{2,1}.*gradient(E_dot',h)' + eij_dot{2,2}.*gradient(E_dot,h) + ...
            (gradient(Ux',h)' + gradient(Uy,h)).*(gradient(E_dot',h)' + gradient(E_dot,h))) + ...
            Ar*E_dot.^(1-1/n).*S.*gradient(S,h);
    end
    
    %%% Solving inverse problems
    Ux_p(:) = Mx\RHS_x(:); % Solutions of the poisson equation
    Uy_p(:) = My\RHS_y(:); % 
    
    Ux_p(:,end) = Ux_old(:,end); % North boundary (Dirichlet)
    Ux_p(:,1) = Ux_old(:,1); % South boundary (Dirichlet)
    Ux_p(end,:) = Ux_old(end,:); % East Boundary (Dirichlet)
    Ux_p(1,:) = Ux_old(1,:); % West Boundary (Dirichlet)
    
    Uy_p(1,:) = 1/3*(4*Uy_p(2,:)-Uy_p(3,:)); % West Boundary (Neumann)
    Uy_p(:,end) = Uy_old(:,end); % North boundary (Dirichlet)
    Uy_p(:,1) = Uy_old(:,1); % South boundary (Dirichlet)
    Uy_p(end,:) = Uy_old(end,:); % East Boundary (Dirichlet)
    
    Ux_old = Ux; % Reassigning velocities
    Uy_old = Uy;

    % Combined with previous solutions for velocity
    if p_step == 1
        Ux_new = Ux_p; % If only one step
        Uy_new = Uy_p; % I.e. first and final step
    else
        Ux_new = alpha.*Ux_p + (1 - alpha).*Ux_old; % Stability condition
        Uy_new = alpha.*Uy_p + (1 - alpha).*Uy_old; % Intermediate step
    end
    
    % Reassigning velocities
    Ux = Ux_new;
    Uy = Uy_new;
    
    % Convergence metric
    beta_arr(p_step,1) = max(Ux_new(:)-Ux_old(:))/max(Ux_old(:)); % Convergence metric, X velocity
    beta_arr(p_step,2) = max(Uy_new(:)-Uy_old(:))/max(Uy_old(:)); % Y velocity
    beta_curr = max(beta_arr(p_step,:)); % Maximum of either convergence metric
    
    % Iteration limit
    if p_step > max_steps % Setting limit on number of steps
        beta_curr = beta;
        figure % Plot failure of convergence
        hold on
        plot(1:max_steps,log10(max(beta_arr,[],2)))
        plot(xlim,log10([beta,beta]),'r--')
        title(['Convergence metric, \beta = 10^{',num2str(log10(beta)),'}'])
        xlabel('Iteration number')
        ylabel('Convergence parameter, log(\beta)')
        box on,grid on
        error(['Solution did not converge in ',num2str(max_steps),' iteration(s)'])
    end
end
end

%         RHS_x = - 3*del__g(Ux,h,'xx') - 3*del__g(Uy,h,'xy') + 2*(1-1/n)*E_dot.^(-1).*...
%             (eij_dot{1,1}.*del__g(E_dot,h,'x') + eij_dot{1,2}.*del__g(E_dot,h,'y') + ...
%             (del__g(Ux,h,'x') + del__g(Uy,h,'y')).*(del__g(E_dot,h,'x') + del__g(E_dot,h,'y'))) + ...
%             Ar*E_dot.^(1-1/n).*S.*del__g(S,h,'x');
%         % Y-direction
%         RHS_y = - 3*del__g(Uy,h,'yy') - 3*del__g(Ux,h,'xy') + 2*(1-1/n)*E_dot.^(-1).*...
%             (eij_dot{2,1}.*del__g(E_dot,h,'x') + eij_dot{2,2}.*del__g(E_dot,h,'y') + ...
%             (del__g(Ux,h,'x') + del__g(Uy,h,'y')).*(del__g(E_dot,h,'x') + del__g(E_dot,h,'y'))) + ...
%             Ar*E_dot.^(1-1/n).*S.*del__g(S,h,'y');