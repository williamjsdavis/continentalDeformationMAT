%% Benchmark performance
clearvars,close all
addpath('src/')

Ar_vector = linspace(0,10,11);
n_vector = linspace(1,5,5);
Nsamples = 5;
timeSamples = zeros(Nsamples,length(Ar_vector),length(n_vector));

for nn = 1:numel(n_vector)
    for nAr = 1:numel(Ar_vector)
        for n = 1:Nsamples
            disp(['nSamp = ',num2str(n)])
            disp(['n = ',num2str(n_vector(nn))])
            disp(['Ar = ',num2str(Ar_vector(nAr))])
            tic
            callSolve(Ar_vector(nAr),n_vector(nn));
            timeSamples(n,nAr,nn) = toc;
        end
    end
end

%% Functions
function callSolve2(Ar,n)
inv(magic(100));
end
function callSolve(ArSet,nSet)

% Settings
[L,u0,g,pc,pm,s0,n,Ar,Nx,dt,nt,S_bound,poisson_set] = simulation_settings();
Ar = ArSet;
Nx = 32;
n = nSet;
nt = 10;

% Setup grids
[Ux,Uy,S,s0,h,X,Y,x,y] = setup_grid(Nx,L,u0,s0,dt);

% Processing
[Ux_new,Uy_new,S_new] = time_solve(Ux,Uy,S,h,n,Ar,dt,S_bound,nt,poisson_set);
disp('Complete!')

end
