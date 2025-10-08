%% Arnold tongues for conventional and non-conventional oscillators.
clear all; close all; clc;

%% Classical Kuramoto model
L = [-1 1; 1 -1];
pcenter = 0.0;

datapoints = 500;
K = linspace(0,1,datapoints);
deltaP = linspace(-2,2,datapoints);
Kuramoto = zeros(datapoints,datapoints);

for i = 1:datapoints
    for j = 1:datapoints
        p = [pcenter + deltaP(j)/2, pcenter - deltaP(j)/2];
        deltaW = abs(p(1) - p(2));

        LyapExp(i,j) = max(K);
        if abs(deltaW/(2*K(i))) <= 1
            Kuramoto(i,j) = 1;
            phi = asin(deltaW/(2*K(i)));
            J = K(i) * cos(phi) * L;
            eigJ = sort(real(eig(J)));
            LyapExp(i,j) = eigJ(1);
        end
    end
end

% Plot
figure(1); subplot(221)
imagesc(LyapExp(end:-1:1,1:end))
% colormap(slanCM('coolwarm'))
colorbar;
m = max(abs(LyapExp(:)));
caxis([-m m])  

axis square
hold on
contour(LyapExp(end:-1:1,1:end), [0 0], 'k', 'LineWidth', 1);

yticks([1 datapoints/2 datapoints])
yticklabels({'1','0.5','0'})
xticks([1 datapoints/2 datapoints])
xticklabels({'-1','0','1'})
xlabel('\Delta p')
ylabel('K')
fontsize(16,"points")

hold off

%% Leaky Kuramoto model
L = [-1 1; 1 -1];
pcenter = 0.6;

datapoints = 500;
K = linspace(0,0.1,datapoints);
deltaP = linspace(-2,2,datapoints);
DampedConsensus = zeros(datapoints,datapoints);

for i = 1:datapoints
    for j = 1:datapoints       
        p = [pcenter + deltaP(j)/2, pcenter - deltaP(j)/2];
        J =  K(i)*L - diag(p);
        eigJ = eig(J);
        LyapExp(i,j) = max(real(eigJ));
    end
end

% Plot
figure(1); subplot(222)
imagesc(LyapExp(end:-1:1,1:end))
% colormap(slanCM('coolwarm'))
colorbar;
m = max(abs(LyapExp(:)));
caxis([-m m])  

axis square
hold on
contour(LyapExp(end:-1:1,1:end), [0 0], 'k', 'LineWidth', 1);

yticks([1 datapoints/2 datapoints])
yticklabels({'1','0.5','0'})
xticks([1 datapoints/2 datapoints])
xticklabels({'-2','0','2'})
xlabel('\Delta p')
ylabel('K')
fontsize(16,"points")

hold off

%% Second-order Kuramoto
L = [-1 1; 1 -1];
pcenter = 0.75;

datapoints = 500;
K = linspace(0,0.1,datapoints);
deltaP = linspace(-2,2,datapoints);
Kuramoto2ndOrder = zeros(datapoints,datapoints);

for i = 1:datapoints
    for j = 1:datapoints        
        p = [pcenter + deltaP(j)/2, pcenter - deltaP(j)/2];
        J = [zeros(2,2) eye(2); K(i)*L -diag(p)];
        eigJ = real(eig(J));
        eigJ = sort(eigJ);
        eigJ(end) = [];
        LyapExp(i,j) = max(real(eigJ));
    end
end

% Plot
figure(1); subplot(223)
imagesc(LyapExp(end:-1:1,1:end))
% colormap(slanCM('coolwarm'))
colorbar;
m = max(abs(LyapExp(:)));
caxis([-m m])  

axis square
hold on
contour(LyapExp(end:-1:1,1:end), [0.01 0.01], 'k', 'LineWidth', 1);

yticks([1 datapoints/2 datapoints])
yticklabels({'1','0.5','0'})
xticks([1 datapoints/2 datapoints])
xticklabels({'-2','0','2'})
xlabel('\Delta p')
ylabel('K')
fontsize(16,"points")

hold off



%% Phase-amplitude model
eps = 0.2;
alpha = 0.7;
c = 0.3
L = [-1 1; 1 -1];
pcenter = 1.5;

datapoints = 500;
K = linspace(0,0.035,datapoints);
deltaP = linspace(-2,2,datapoints);
ArnoldPhaseAmp = zeros(datapoints,datapoints);

for i = 1:datapoints
    for j = 1:datapoints        
        p = [pcenter + deltaP(j)/2, pcenter - deltaP(j)/2];
        J = [-K(i)*L-c*diag(p)+alpha*eye(2) eye(2); -eps*eye(2), -diag(p)+alpha*eye(2)] ;
        eigJ = eig(J);
        LyapExp(i,j) = max(real(eigJ));
    end
end

figure(1); subplot(224)
imagesc(LyapExp(end:-1:1,1:end))
% colormap((slanCM('coolwarm')))
colorbar;
m = max(abs(LyapExp(:)));
caxis([-m m])  

axis square
hold on
contour(LyapExp(end:-1:1,1:end), [0 0], 'k', 'LineWidth', 1);

yticks([1 datapoints/2 datapoints])
yticklabels({'1','0.5','0'})
xticks([1 datapoints/2 datapoints])
xticklabels({'-2','0','2'})
xlabel('\Delta p')
ylabel('K')
fontsize(16,"points")

hold off
