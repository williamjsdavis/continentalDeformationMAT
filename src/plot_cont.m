%Plot continental model
%   William Davis, 12/12/17
%
%   Notes:
%   Plots results from continental deformation modelling
%
%   Inputs:
%   - 
%
%   Problems:
%   - No scale for tensor plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings
save_figure = 0; % Save figure?
filename = strcat('Ar',num2str(Ar),'n',num2str(n),'t',num2str(nt)); % File to save as
cbar_range1 = 1E-15*[-2,20]; % Minimum and maximum isotropic strain rates (chosen) [s^-1]
cbar_range2 = [20,65]; % Minimum and maximum thicknesses (chosen) [km]
x_section = [1,4,6,12]; % Where to plot cross sections (change this with grid size)
%x_section = [2,8,16,24]; % Where to plot cross sections (change this with grid size)

% Colourbars
cbar = redblue(); % Set colourbar to red/blue

% Isotropic strain rate colourbar
cbar_coeff = (0-cbar_range1(1))/(cbar_range1(2)-cbar_range1(1)); % Colourbar coeficient
cbar1 = cbar(floor(64*(1-2*cbar_coeff)/(2-2*cbar_coeff)-1)+3:64,:); % Scale so that original thickness is white.

% Thickness colourbar
cbar_coeff = (L*s0-cbar_range2(1))/(cbar_range2(2)-cbar_range2(1)); % Colourbar coeficient
cbar2 = cbar(floor(64*(1-2*cbar_coeff)/(2-2*cbar_coeff)-1)+2:64,:); % Scale so that original thickness is white.

% Velocity plot
logmax = log10(u0);
logmin = log10(1E-10);
npts = 11;
colourbar_log = logspace(logmin,logmax,npts);
colourbar_lin = linspace(u0,0,npts);
cbar3 = cbar(32:64,:);

%% Dimensionalise parameters
% Properties

% Strain rate tensor
[e11_dot,e12_dot,e21_dot,e22_dot] = strain_rate(Ux_new,Uy_new,h,1);

% Isotropic strain rate
eIso_dot = -(e11_dot + e22_dot); 

% Isotropic strain rate [s^-1]
eIso_dot_dim = 1E-6/(365.25*24*60*60)*u0/L*eIso_dot; 

% Second invariant of the strain rate tensor
E_dot = sqrt(2)*sqrt(e11_dot.^2 + e22_dot.^2 + e12_dot.^2 + ...
    e11_dot.*e22_dot); 

% Deviatoric stress tensor [MPa]
[t11_dim,t12_dim,t21_dim,t22_dim] = ...
    strain_rate(Ux_new,Uy_new,h,1*g*pc*L*(1-pc/pm)/Ar*E_dot.^(1/n-1)); 
tij_dim = cell(2);
tij_dim{1,1} = t11_dim;
tij_dim{1,2} = t12_dim;
tij_dim{2,1} = t21_dim;
tij_dim{2,2} = t22_dim;

% Strain rate tensor [s^-1]
[e11_dot_dim,e12_dot_dim,e21_dot_dim,e22_dot_dim] = ...
    strain_rate(Ux_new,Uy_new,h,1E-6/(365.25*24*60*60)*u0/L); 
eij_dot_dim = cell(2);
eij_dot_dim{1,1} = e11_dot_dim;
eij_dot_dim{1,2} = e12_dot_dim;
eij_dot_dim{2,1} = e21_dot_dim;
eij_dot_dim{2,2} = e22_dot_dim;

%% Plotting
figure('Position',[204,140,1237,665])
fig_dims = [2,3];

sp1 = subplot(fig_dims(1),fig_dims(2),1);
hold on
contourf(L*X,L*Y,eIso_dot_dim,'LineColor','none')
contour(L*X,L*Y,eIso_dot_dim,[cbar_range1(1):5E-15:cbar_range1(2)],'k')
%surf(L*Xint,L*Yint,u0*Uxint)
title('Isotropic strain rate')
xlabel('X distance, km')
ylabel('Y distance, km')
colormap(sp1,cbar1)
cb_sp1 = colorbar('eastoutside');
cb_sp1.Label.String = 'Isotropic strain rate, s^{-1}';
cb_sp1.TickLength = 0.1;
caxis(cbar_range1)
axis equal
axis([0,L*Nx,0,L*Nx])
box on, grid on
plot([0,0,L*Nx,L*Nx],[0,L*Nx,L*Nx,0],'k')
set(gca,'TickDir','out')

subplot(fig_dims(1),fig_dims(2),4)
plot_elipsiod(eij_dot_dim,X,Y,L,1E15)
title('Principle strain rates')
xlabel('X distance, km')
ylabel('Y distance, km')
axis equal
axis([0,L*Nx,0,L*Nx])
box on, grid on
plot([0,0,L*Nx,L*Nx],[0,L*Nx,L*Nx,0],'k')
set(gca,'TickDir','out')

subplot(fig_dims(1),fig_dims(2),2)
plot_elipsiod(tij_dim,X,Y,L,1E-3) % -3
title('Horizontal deviatoric stress')
xlabel('X distance, km')
ylabel('Y distance, km')
axis equal
axis([0,L*Nx,0,L*Nx])
box on, grid on
plot([0,0,L*Nx,L*Nx],[0,L*Nx,L*Nx,0],'k')
set(gca,'TickDir','out')

vector = [1:(Nx/16):Nx,1:(Nx/16):Nx];
sp5 = subplot(fig_dims(1),fig_dims(2),5);
hold on
contourf(L*X,L*Y,u0*sqrt(Ux_new.^2+Uy_new.^2),'LineColor','none');
quiver(L*X(vector,vector),L*Y(vector,vector),...
    u0*Ux_new(vector,vector)./sqrt(Ux_new(vector,vector).^2+Uy_new(vector,vector).^2),...
    u0*Uy_new(vector,vector)./sqrt(Ux_new(vector,vector).^2+Uy_new(vector,vector).^2),'LineWidth',1.5,'Color','k')
title('Crustal velocity')
xlabel('X distance, km')
ylabel('Y distance, km')
colormap(sp5,cbar3)
cb_sp5 = colorbar('eastoutside');
cb_sp5.Label.String = 'Velocity, mm yr^{-1}';
axis equal
axis([0,L*Nx,0,L*Nx])
box on, grid on
plot([0,0,L*Nx,L*Nx],[0,L*Nx,L*Nx,0],'k')
set(gca,'TickDir','out')

sp3 = subplot(fig_dims(1),fig_dims(2),3);
hold on
contourf(L*X,L*Y,L*S_new,'LineColor','none')
contour(L*X,L*Y,L*S_new,cbar_range2(1):5:cbar_range2(2),'k')
for i = 1:4
    plot(L*X(x_section(i),:),L*Y(x_section(i),:),'k--')
    text(L*(0.5+x(x_section(i))),0.95*Nx*L,strcat(char(64+i),'-',char(64+i),char(39)))
end
title('Crustal thickness')
xlabel('X distance, km')
ylabel('Y distance, km')
colormap(sp3,cbar2)
cb_sp3 = colorbar('eastoutside');
cb_sp3.Label.String = 'Crustal thickness, km';
cb_sp3.TickLength = 0.1;
caxis(cbar_range2)
view(0,90)
axis equal
axis([0,L*Nx,0,L*Nx])
box on, grid on
plot([0,0,L*Nx,L*Nx],[0,L*Nx,L*Nx,0],'k')
set(gca,'TickDir','out')

subplot(fig_dims(1),fig_dims(2),6)
plot_section(x_section,y,S_new,pc/pm,L)
%title(['Y velocity, mm/yr (Time = ',sprintf('%0.1f',nt*dt*u0/L),' Ma)'])
title('Crustal section profiles')
xlabel('X distance, km')
ylabel('Elevation above undisturbed crust, km')

% Overall title
dim = [0.42,0.68,0.3,0.3];
str = ['Properties at time = ',sprintf('%0.1f',nt*dt*u0/L),' Ma'];
annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none',...
    'FontSize',20,'FontWeight','bold');

%% Save figure
if save_figure == 1
    addpath('/Users/williamdavis/Documents/masters_project/matlab/altmany-export_fig-2763b78')
    set(gcf,'Color','white')
    export_fig(gcf,filename,'-painters','-pdf' );
    movefile(strcat(filename,'.pdf'),'writeup/figures')
    close all
end

% Old plot
%{
figure('Position',[204,140,1237,665])
fig_dims = [2,3];

subplot(fig_dims(1),fig_dims(2),1)
surf(L*X,L*Y,u0*Uy)
%surf(L*Xint,L*Yint,u0*Uxint)
title('Y velocity initial')
xlabel('X distance, km')
ylabel('Y distance, km')
cb_vel1 = colorbar('eastoutside');
cb_vel1.Label.String = 'Velocity, mm/yr';
view(0,90)
caxis([0,u0])
axis equal
axis([0,L*Nx,0,L*Nx])

sp2 = subplot(fig_dims(1),fig_dims(2),4);
hold on
contourf(L*X,L*Y,L*S_new,'LineColor','none')
contour(L*X,L*Y,L*S_new,[cbar_range(1):5:cbar_range(2)],'k')
%surf(L*Xint,L*Yint,L*Sint)
title(['Crustal thickness, km (Time = ',sprintf('%0.1f',nt*dt*u0/L),' Ma)'])
xlabel('X distance, km')
ylabel('Y distance, km')
colormap(sp2,cbar)
cb_thi = colorbar('eastoutside');
cb_thi.Label.String = 'Crustal thickness, km';
cb_thi.TickLength = 0.1;
caxis(cbar_range)
view(0,90)
axis equal
axis([0,L*Nx,0,L*Nx])

subplot(fig_dims(1),fig_dims(2),2)
surf(L*X,L*Y,u0*Ux_1)
%surf(L*Xint,L*Yint,u0*Uyint)
title('X velocity')
xlabel('X distance, km')
ylabel('Y distance, km')
view(0,90)
caxis([0,u0])
axis equal
axis([0,L*Nx,0,L*Nx])

subplot(fig_dims(1),fig_dims(2),5)
surf(L*X,L*Y,u0*Uy_1)
%surf(L*Xint,L*Yint,u0*vvec)
title('Y velocity')
xlabel('X distance, km')
ylabel('Y distance, km')
cb_vel2 = colorbar('eastoutside');
cb_vel2.Label.String = 'Velocity, mm/yr';
view(0,90)
axis equal
axis([0,L*Nx,0,L*Nx])

subplot(fig_dims(1),fig_dims(2),3)
surf(L*X,L*Y,u0*Ux_2)
%surf(L*Xint,L*Yint,u0*Uyint)
%title('X velocity')
title(['X velocity, mm/yr (Time = ',sprintf('%0.1f',nt*dt*u0/L),' Ma)'])
xlabel('X distance, km')
ylabel('Y distance, km')
view(0,90)
caxis([0,u0])
axis equal
axis([0,L*Nx,0,L*Nx])

subplot(fig_dims(1),fig_dims(2),6)
surf(L*X,L*Y,u0*Uy_2)
%surf(L*Xint,L*Yint,u0*vvec)
%title('Y velocity')
title(['Y velocity, mm/yr (Time = ',sprintf('%0.1f',nt*dt*u0/L),' Ma)'])
xlabel('X distance, km')
ylabel('Y distance, km')
cb_vel3 = colorbar('eastoutside');
cb_vel3.Label.String = 'Velocity, mm/yr';
view(0,90)
axis equal
axis([0,L*Nx,0,L*Nx])
%}