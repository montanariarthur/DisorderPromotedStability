%% Disorder-promoted stability in several dynamical systems models,
% including first- and second-order oscillators.
clear all; close all; clc;
addpath([pwd,'/Models/'])
warning('off')

poolobj = gcp('nocreate');
delete(poolobj);
ncores = input('Number of cores: ')
parpool(ncores);

% Choose one of the dynamical systems listed below
system = 'power-grid'
% van-der-pol
% power-grid (second-order Kuramoto oscillator)
% neuron (FitzHugh-Nagumo oscillators)
% metamaterial (spring-mass-damper networks)
% multi-agent (consensus, flocking, opinion)
% josephson-junction
% phase-amplitude
% first-order (leaky integrator)

Nparam = 100;       % number of parameter steps
Nmc = 100;          % number of optimization realizations (chooses best)

%% Default parameters

% Network
N = 100                         % number of oscillators
A = ones(N,N) - eye(N);
A = 1*A.*rand(N,N);
A = A.*(rand(N,N) <= 1);
param.A = A;
param.L = diag(sum(A,2)) - A;   % graph Laplacian 

% System parameters
switch system
    case 'van-der-pol'
        K = 1;                                 % coupling strength
        param_list = linspace(0,N,Nparam);    % nonlinearity factor
        param.yeq = zeros(N,1);                % equilibrium
        param.xeq = zeros(N,1);
        zeroLE = 0;

    case 'power-grid'
        param_list = linspace(0,10*N,Nparam);  % damping
        zeroLE = 1;

    case 'neuron'
        param_list = linspace(0,N,Nparam); % damping
        param.a = 1.5;                      % pump
        param.eps = 1;                      % gain
        param.K = 1;                        % coupling strength
        zeroLE = 1;

    case 'metamaterial'
        param_list = linspace(0,N,Nparam);  % damping
        param.m = ones(N,1);                % mass
        param.k = 0.5;                      % spring stiffness
        zeroLE = 0;

    case 'multi-agent'
        param_list = linspace(0,N,Nparam); % damping
        param.gamma = 1;                    % factor
        param.k = 0.1;                      % couplign strength
        zeroLE = 0;

    case 'josephson-junction'
        param_list = linspace(0,N,Nparam); % resistance
        param.C = ones(N,1);                % capacitance
        param.Ic = 1.5*ones(N,1);           % critical current
        param.Ib = 1.0*ones(N,1);           % bias current
        zeroLE = 0;

    case 'phase-amplitude'
        param_list = linspace(0,N,Nparam);  % damping
        param.eps = 0.1;                    
        param.gamma = 1;                    
        zeroLE = 0;

    case 'first-order'
        param_list = linspace(0,N,Nparam);
        param.k = 1;
        zeroLE = 1;
end

%% Stability optimization for homogeneous system
for k = 1:Nparam
    phom(:,k) = param_list(k)/2 * ones(N,1);
    [LyapExp_hom(k),xeq_hom(:,k)] = findLyapExp(system,param,N,phom(:,k),zeroLE);
end

%% Stability optimization for heterogeneous system

% Optimization over many trials
ptrial = zeros(N,Nmc);
minLyapExp = zeros(Nmc,1);
[~,argmin] = min(LyapExp_hom); 
p0 = param_list(argmin)/2;
parfor i = 1:Nmc            % optimization
    % Nelder-Mead simplex method
    % [ptrial(:,i),minLyapExp(i)] = fminsearch(@(p)findLyapExp(system,param,N,p,zeroLE), p0+0.5*randn(N,1));

    % Quasi-Newton (BFGS) method
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'MaxIterations', 1e4, 'Display','off');
    [ptrial(:,i), minLyapExp(i)] = fminunc(@(p)findLyapExp(system,param,N,p,zeroLE), p0+0.5*randn(N,1), options); 
end

% Picks best trial
[~,minindex] = min(minLyapExp);
pvec = ptrial(:,minindex);

% Calculates stability over wide range
for k = 1:Nparam            % computation over parameter range
    phet(:,k) = param_list(k) * pvec / norm(pvec);
    [LyapExp_het(k),xeq_het(:,k)] = findLyapExp(system,param,N,phet(:,k),zeroLE);
end

%% Plot
figure(1);

subplot(131);  hold on;
plot(param_list, LyapExp_hom, 'red')
plot(param_list, LyapExp_het, 'blue')
xlabel('parameter')
ylabel('\lambda_m_a_x')

subplot(132);  hold on;
plot(param_list, vecnorm(xeq_hom), 'red')
plot(param_list, vecnorm(xeq_het), 'blue')
xlabel('parameter')
ylabel('equilibrium')

subplot(133);
histogram(pvec/norm(pvec)); hold on;
[f, xi] = ksdensity(pvec/norm(pvec));  % Estimate the density
plot(xi, f, 'LineWidth', 2);           % Plot the smooth curve
xlabel('pvec');
ylabel('Density');

fontsize(18,"points")



