function [  ] = plot_section( x,y,S,pcpm,L )
%Plots elipsiods for tensors
%   William Davis, 16/12/17
%
%   Notes: Designed for 2x2 matrix tensor
%   
%
%   Inputs:
%   - "x"                       Index of section to plot
%   - "y"                       Array of y direction distances, []
%   - "S"                       Crustal thickness matrix, []
%   - "pcpm"                    Ratio of crustal density to mantle density, []
%   - "dt"                      Time-step, []
%   - "S_bound"                 South boundary type: 'const' or 'neu'
%
%   Problems:
%   - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = length(y);
bottom = [0.75,0.50,0.25,0.00];
scale = 5;

hold on
%fill(L*[y,fliplr(y)],centre+scale*L*[S(x,:)*(1-pcpm),fliplr(-S(x,:)*(pcpm))],'k')
for i = 1:4
    %fill(L*[y,fliplr(y)],bottom(i)+scale*L*[S(x(i),:)*(1-pcpm),zeros(size(y))],'r')
    fill(L*[y,fliplr(y)],bottom(i)+[scale*(S(x(i),:)-0.35)*(1-pcpm)+0.125,zeros(size(y))],'r')
    plot(L*y,bottom(i)+scale*L*0.35*(1-pcpm)*ones(size(y)),'k--')
end

axis([0,L*Nx,0,1])
axis square
box on
lnwidth = 1.5;
plot([0,0,L*Nx,L*Nx,0],[0.00,0.25,0.25,0.00,0.00],'k','LineWidth',lnwidth)
plot([0,0,L*Nx,L*Nx,0],[0.25,0.50,0.50,0.25,0.25],'k','LineWidth',lnwidth)
plot([0,0,L*Nx,L*Nx,0],[0.50,0.75,0.75,0.50,0.50],'k','LineWidth',lnwidth)
plot([0,0,L*Nx,L*Nx,0],[0.75,1.00,1.00,0.75,0.75],'k','LineWidth',lnwidth)
x_offset = 0.015;
y_offset = 0.04;

for i = 1:4
    text(L*Nx*x_offset,0.25+bottom(i)-y_offset,strcat(char(64+i),'-',char(64+i),char(39)),'FontSize',14,'FontWeight','bold')
end
set(gca,'TickDir','out')
set(gca,'ytick',[])

end

