%Plot 3D continental model
%   William Davis, 20/07/20
%
%   Notes:
%   Plots 3D results from continental deformation modelling
%
%   Inputs:
%   - 
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mirror domain
thick2elev = (1 - pc/pm);
Escaled = L*thick2elev*[flipud(S_new);S_new];
Xscaled = L*[X;X+max(X(:))+X(2,2)];
Yscaled = L*[flipud(Y);Y];

% Plot
figure('Position',[399,259,983,517])
surf(Xscaled,Yscaled,Escaled)
zlim([2.5,7])
zlabel('Elevation, km')
xlabel('X distance, km')
ylabel('Y distance, km')

% View settings
load('DEMcolormap.mat')
view(145,13)
colormap(mymap)
caxis([2.5 6.5])
set(gca,'DataAspectRatio',[300 300 1])

% Overall title
dim = [0.24,0.68,0.3,0.3];
timeMa = nt*dt*u0/L;
str = ['Topography at time = ',sprintf('%0.1f',timeMa),' Ma'];
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
           'EdgeColor','none','FontSize',20,'FontWeight','bold');

