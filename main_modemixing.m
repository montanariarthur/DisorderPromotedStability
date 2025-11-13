%% Mode mixing for second-order systems.
clear all; close all; clc;
addpath([pwd,'/Models/'])

%% Second-order Kuramoto model
system = 'power-grid'

% Network parameters
N = 3;
A = ones(N,N) - eye(N);
A = A.*rand(N,N);
A = (A+A')/2

% Matrix used in the plot
% A = [      0    0.9096    0.2027;
%       0.9096         0    0.3222;
%       0.2027    0.3222         0];

L = diag(sum(A,2)) - A;
param.L = L;
[VL,DL] = eig(L);

% Heterogeneity over parameters
datapoints = 1000;
param_list = linspace(0,5,datapoints);

% Homogeneous optimum damping
for k = 1:datapoints
    bhom(:,k) = param_list(k) * ones(N,1);
    [LyapExp_hom(k),xeq_hom(:,k)] = findLyapExp(system,param,N,bhom(:,k),1);
end
[~,argmin] = min(LyapExp_hom); 
bhom_opt = param_list(argmin)*ones(N,1);

% Heterogeneous optimum damping
Nmc = 10;
minLyapExp = zeros(Nmc,1);
btrial = zeros(N,Nmc);
parfor i = 1:Nmc            % optimization
    % Quasi-Newton (BFGS) method
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'MaxIterations', 1e4, 'Display','off');
    [btrial(:,i), minLyapExp(i)] = fminunc(@(p)findLyapExp(system,param,N,p,1), bhom_opt+0.5*randn(N,1), options); 
end
[~,minindex] = min(minLyapExp); % picks best trial
bhet = btrial(:,minindex);

% Disordered damping
disorder = randn(N,1);


%% Jacobian construction
datapoints = 200;
epsilon = linspace(0,2,datapoints);
eigJ = [];
alignment = [];

for i = 1:datapoints

    % Optimized parameters
    b = bhom_opt + epsilon(i)*(bhet - bhom_opt);

    % Jacobian matrix
    Jseq(:,:,i) = [zeros(N,N) eye(N); -L -diag(b)];
end

% Compute eigenvalues and eigenvectors
[Vseq,Dseq] = eigenshuffle(Jseq);

% Eigenvalues
for i = 1:datapoints
    eigJ(:,i) = real(Dseq(:,i));
    if i == 1
        [LyapExp(i),minind] = max(real(Dseq(2:end,i)));
    end

    % Dot products
    for k = 1:2*N
        alignment(k,i) = Vseq(:,minind+1,i)'*Vseq(:,k,i);
    end
end

%% Plot
figure(1); subplot(121)
plot(epsilon,eigJ);
xlabel('path length')
ylabel('real part of eigenvalues')
fontsize(16,"points")
% ylim([-2.05 0.05])
legend

figure(1); subplot(122)
plot(epsilon,abs(alignment));
xlabel('disorder level')
ylabel('dot products')
fontsize(16,"points")
ylim([-0.05 1.05])
legend
