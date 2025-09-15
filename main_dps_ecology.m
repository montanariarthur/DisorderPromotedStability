%% Disorder-promoted stability in ecological models
clear all; close all; clc;
addpath([pwd,'/Models/'])
warning('off')

poolobj = gcp('nocreate');
delete(poolobj);
ncores = input('Number of cores: ')
parpool(ncores);

system = input('Choose model (implicit; explicit): ')
p_edge = input('Choose edge probability: ')
p_sign = input('Choose sign probability: ')
parameter = input('Choose control parameter (network; connectivity; sign; growth; none): ')

%% Default parameters

N = 100;            % number of nodes

Nhet = 100;         % number of heterogeneity steps
sigma = linspace(0,2,Nhet);     % heterogeneity interval
growth_het = 10;    % growth heterogeneity

if N <= 100         % number of Monte Carlo realizations
    Nmc = 1000; 
else
    Nmc = 100;
end

% Control parameter
switch parameter
    case 'network'
        param_range = [2 5 10 20 50 100 200 500 1000];
    case 'connectivity'
        param_range = [0.01 0.02 0.05 0.1 0.2 0.5 1];
    case 'sign'
        param_range = 0:0.05:1;
    case 'growth'
        param_range = 0:1:20;
    case 'none'
        param_range = 1;
end

%% Stability analysis

% Parameter loop
for j = 1:length(param_range)
    disp(['Control parameter (', parameter, '): ', num2str([param_range(j)])])

    switch parameter
        case 'network'
            N = param_range(j);
        case 'connectivity'
            p_edge = param_range(j);
        case 'sign'
            p_sign = param_range(j);
        case 'growth'
            growth_het = param_range(j);
    end

    % Independent realizations loop
    LyapExp_run = zeros(Nmc,Nhet);
    feasible_run = ones(Nmc,Nhet);
    parfor i = 1:Nmc
        
        % Adjacency matrix
        Adj = econetwork(N,p_edge,p_sign)

        % Self-regulation
        b = 1 + growth_het*rand(N,1);      % growth rate
        alpha = - 1*ones(N,1);             % self interaction

        % Disorder level
        for k = 1:Nhet

            switch system

                % Implicit model (linear ODE xdot = Ax)
                case 'implicit'
                    A = (1/N)*sigma(k)*Adj;
                    B = - diag(b);
                    J = B + A;
                
                % Explicit model (Lotka-Volterra ODE)
                case 'explicit'
                    % Adjacency matrix
                    A = (1/N)*(sigma(k)*Adj) + alpha.*eye(N);

                    % Finds a feasible equilibrium of the LK model
                    xeq = - inv(A)*b;
                    X = diag(xeq);

                    % Feasible equilibrium?
                    if sum(xeq >= 0) == N
                        feasible_run(i,k) = 1;
                    else
                        feasible_run(i,k) = NaN;
                    end

                    % Jacobian
                    J = X*A;
            end

            % Eigenvalues
            eigJ = eig(J);

            % Largest Lyapunov exponent
            LyapExp_run(i,k) = max(real(eigJ));
        end
    end

    % Save data
    feasible{j} = feasible_run;
    LyapExp{j} = LyapExp_run;
end

%% Performance analysis
for j = 1:length(param_range)
    % Stability improvement
    for i = 1:Nmc
        stabilityimprovement{j}(i,:) = ( LyapExp{j}(i,:) - LyapExp{j}(i,1) );% / LyapExp{j}(iter,1);
    end

    % Fraction of stabilized systems by disorder
    stabilizedfraction{j} = sum(LyapExp{j}(:,:)'.*feasible{j}(:,:)' < LyapExp{j}(:,1)' ,2)/Nmc;
end

%% Plot
Nsamps = 10;  % number of plotted curves
j = 1;

figure(1)
subplot(2,2,1)
plot(sigma,LyapExp{j}(1:Nsamps,:)'.*feasible{j}(1:Nsamps,:)','LineWidth',2)
xlim([0 2]); %ylim([-2 .2]);
xlabel('heterogeneity \sigma')
ylabel('\lambda_m_a_x')

subplot(2,2,2)
plot(sigma,nanmean(stabilityimprovement{j},1),'LineWidth',2)
xlim([0 2]); %ylim([-2 .2]);
xlabel('heterogeneity \sigma')
ylabel('relative improvement')

subplot(2,2,3)
plot(sigma,stabilizedfraction{j},'LineWidth',2)
xlim([0 2]); %ylim([-2 .2]);
xlabel('heterogeneity \sigma')
ylabel('\lambda_m_a_x')

fontsize(18,"points")
