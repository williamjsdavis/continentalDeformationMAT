classdef TVSfield < handle
    %Thin Viscous Sheet field
    %   Detailed explanation goes here
    
    properties
        Ux
        Uy
        S
        X
        Y
        simSettings
        otherProp
    end
    
    methods
        function obj = TVSfield()
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            
            % Settings
            simulationSettings(obj)
            
            
            % Setup grids
            setupGrid(obj)
        end
        
        function simulationSettings(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
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
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
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
        function timeSolve(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            % Set up Poisson stencils
            Nx = obj.simSettings.Nx;
            [Mx,My] = setup_poisson(Nx,obj.otherProp.h);
            
            % Time-stepping
            for i = 1:obj.simSettings.nt
                % Solve for velocity
                obj.poissonVel(Mx,My);
                
                % Solve for thickness
                obj.upwindS();
                
                % Update
%                obj.Ux = Ux_new; obj.Uy = Uy_new; obj.S = S_new;
%                 Ux = Ux_new; Uy = Uy_new; S = S_new;
                
                % Print message
                pctTime = 100*i/obj.simSettings.nt;
                disp(['Progress: ',sprintf('%5.2f',pctTime),'%'])
            end
        end
        function poissonVel(obj,Mx,My)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            % New velocities
            [obj.Ux,obj.Uy] = poisson_vel(obj.Ux,obj.Uy,...
                                          Mx,My,...
                                          obj.S,obj.otherProp.h,...
                                          obj.simSettings.n,...
                                          obj.simSettings.Ar,...
                                          obj.simSettings.poisson_set); 
        end
        function upwindS(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            % New thickness
            obj.S = upwind_s(obj.Ux,obj.Uy,obj.S,...
                             obj.otherProp.h,...
                             obj.simSettings.dt,...
                             obj.simSettings.S_bound); 
        end
        function [] = plot6(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            plot_cont;
        end
    end
end

