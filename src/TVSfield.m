classdef TVSfield < handle
    %Thin Viscous Sheet field
    %   Detailed explanation goes here
    
    properties
        Ux
        Uy
        S
        X
        Y
        stencil
        simSettings
        otherProp
    end
    
    methods
        function obj = TVSfield()
            %Field constructor
            %   Setup grids and stencils
            
            % Settings
            obj.simulationSettings()
            
            % Setup grids
            obj.setupGrid()
            
            % Set up Poisson stencils
            obj.setupPoisson()
        end
        
        function simulationSettings(obj)
            %Getting simulation settings
            %   Parse input and pass to functions
            [L,u0,g,pc,pm,s0,n,Ar,Nx,dt,nt,S_bound,poisson_set] = ...
                simulation_settings();
            obj.simSettings.L = L;
            obj.simSettings.u0 = u0;
            obj.simSettings.g = g;
            obj.simSettings.pc = pc;
            obj.simSettings.pm = pm;
            obj.simSettings.s0 = s0;
            obj.simSettings.n = n;
            obj.simSettings.Ar = Ar;
            obj.simSettings.Nx = Nx;
            obj.simSettings.dt = dt;
            obj.simSettings.nt = nt;
            obj.simSettings.S_bound = S_bound;
            obj.simSettings.poisson_set = poisson_set;
        end
        function setupGrid(obj)
            %Set up field grids
            %   Parse input and pass to functions
            [obj.Ux,obj.Uy,obj.S,s0,h,obj.X,obj.Y,x,y] = ...
                setup_grid(obj.simSettings.Nx,...
                           obj.simSettings.L,...
                           obj.simSettings.u0,...
                           obj.simSettings.s0,...
                           obj.simSettings.dt);
            obj.otherProp.s0 = s0;
            obj.otherProp.h = h;
            obj.otherProp.x = x;
            obj.otherProp.y = y;
        end
        function setupPoisson(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            [Mx,My] = setup_poisson(obj.simSettings.Nx,obj.otherProp.h);
            obj.stencil.Mx = Mx;
            obj.stencil.My = My;
        end
        function timeSolve(obj)
            %Solves the solution in space and time
            %   Parse input and pass to functions
            
            % Time-stepping
            for i = 1:obj.simSettings.nt
                % Solve for velocity
                obj.poissonVel(obj.stencil.Mx,obj.stencil.My);
                
                % Solve for thickness
                obj.upwindS();
                
                % Print message
                pctTime = 100*i/obj.simSettings.nt;
                disp(['Progress: ',sprintf('%5.2f',pctTime),'%'])
            end
        end
        function poissonVel(obj,Mx,My)
            %Poisson equation solve for velocity
            %   Parse input and pass to function

            % New velocities
            [obj.Ux,obj.Uy] = poisson_vel(obj.Ux,obj.Uy,...
                                          Mx,My,...
                                          obj.S,obj.otherProp.h,...
                                          obj.simSettings.n,...
                                          obj.simSettings.Ar,...
                                          obj.simSettings.poisson_set); 
        end
        function upwindS(obj)
            %Upwind scheme to solve for thickness
            %   Parse input and pass to function

            % New thickness
            obj.S = upwind_s(obj.Ux,obj.Uy,obj.S,...
                             obj.otherProp.h,...
                             obj.simSettings.dt,...
                             obj.simSettings.S_bound); 
        end
        function plot6(obj)
            %Plot diagnostic results to six subplots
            %   Passes object to plotting function, showing various parameters.
            plot_cont;
        end
        function plot3D(obj)
            %Plot elevation
            %   Plots 3D results from continental deformation modelling
            
            % Mirror domain
            thick2elev = (1 - obj.simSettings.pc/obj.simSettings.pm);
            Lscale = obj.simSettings.L;
            Xoffset = max(obj.X(:))+obj.X(2,2);
            
            Xscaled = Lscale*[obj.X;obj.X+Xoffset];
            Yscaled = Lscale*[flipud(obj.Y);obj.Y];
            Escaled = Lscale*thick2elev*[flipud(obj.S);obj.S];
            
            % Plot
            figure('Position',[399,259,983,517])
            surf(Xscaled,Yscaled,Escaled)
            zlim([2.5,7])
            zlabel('Elevation, km')
            xlabel('X distance, km')
            ylabel('Y distance, km')
            
            % View settings
            load 'DEMcolormap.mat' mymap
            view(145,13)
            colormap(mymap)
            caxis([2.5 6.5])
            set(gca,'DataAspectRatio',[300 300 1])
            
            % Overall title
            dim = [0.24,0.68,0.3,0.3];
            timeMa = obj.simSettings.nt*obj.simSettings.dt*...
                obj.simSettings.u0/Lscale;
            str = ['Topography at time = ',sprintf('%0.1f',timeMa),' Ma'];
            annotation('textbox',dim,'String',str,'FitBoxToText','on',...
                'EdgeColor','none','FontSize',20,'FontWeight','bold');
        end
    end
end

