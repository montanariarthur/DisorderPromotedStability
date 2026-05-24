%% Disorder scan for phase-amplitude oscillators
clear; clc; close all;

%% Parameters
N = 10;
A = 0.3 * directed_ring_circulant(N, 1, 2, 1);

b0 = 0.5;
epsVal = 1;
omega = 0;

sigmaVals = linspace(0, 0.5, 31);
nReal = 1000;

nIC = 1000;
tFinal = 100;

ampTol = 1e-2;
phaseTol = 1e-2;

solverTol = 1e-10;
resTol = 1e-8;
maxIter = 1000;

%% Homogeneous locked-state initial guess
kMean = mean(sum(A, 2));
rho0 = max(1e-6, 1 + epsVal * kMean / b0) * ones(N, 1);
theta0 = zeros(N, 1);
Omega0 = omega + mean(rho0) - 1;

deltaB = randn(nReal, N);

lambdaTrans = NaN(nReal, numel(sigmaVals));
basinFraction = NaN(nReal, numel(sigmaVals));

%% Main scan
parfor rr = 1:nReal
    [lambdaTrans(rr,:), basinFraction(rr,:)] = run_realization( ...
        A, deltaB(rr,:).', sigmaVals, b0, epsVal, omega, ...
        rho0, theta0, Omega0, nIC, tFinal, ...
        ampTol, phaseTol, solverTol, resTol, maxIter);
end

%% Summary
lambdaTrans_mean = mean(lambdaTrans, 1, 'omitnan');
lambdaTrans_std  = std(lambdaTrans, 0, 1, 'omitnan');

basin_mean = mean(basinFraction, 1, 'omitnan');
basin_std  = std(basinFraction, 0, 1, 'omitnan');

% Compare every disorder point with the homogeneous sigma = 0 case.
tol = 1e-12;

validLyap = isfinite(lambdaTrans) & isfinite(lambdaTrans(:,1));
validBasin = isfinite(basinFraction) & isfinite(basinFraction(:,1));

betterLyap = lambdaTrans < lambdaTrans(:,1) - tol;
betterBasin = basinFraction > basinFraction(:,1) + tol;

fracBetterLyap = sum(betterLyap & validLyap, 1) ./ sum(validLyap, 1);
fracBetterBasin = sum(betterBasin & validBasin, 1) ./ sum(validBasin, 1);

%% Plot
figure(1)

subplot(221)
shaded_curve(sigmaVals, lambdaTrans_mean, lambdaTrans_std);
yline(0, '--');
xlim([0 0.5])
ylim([-0.7 0.7])
xlabel('\sigma')
ylabel('\Lambda_{max}')
box on

subplot(223)
shaded_curve(sigmaVals, basin_mean, basin_std);
ylim([-0.05 1.05])
xlim([0 0.5])
xlabel('\sigma')
ylabel('basin size')
box on

subplot(222)
plot(sigmaVals, fracBetterLyap, 'LineWidth', 1.2);
ylim([0 1])
xlim([0 0.5])
xlabel('\sigma')
ylabel('fraction of systems with higher stability')
box on

subplot(224)
plot(sigmaVals, fracBetterBasin, 'LineWidth', 1.2);
ylim([0 1])
xlim([0 0.5])
xlabel('\sigma')
ylabel('fraction of systems with larger basin')
box on

fontsize(16, "points")

%% Functions
function [lambdaTrans, basinFraction] = run_realization( ...
    A, deltaB, sigmaVals, b0, epsVal, omega, ...
    rhoGuess, thetaGuess, OmegaGuess, nIC, tFinal, ...
    ampTol, phaseTol, solverTol, resTol, maxIter)
% Runs numerical continuation to track the linear stability and basin sizes
% for the frequency-synchronized state.

N = size(A, 1);
nSigma = numel(sigmaVals);

lambdaTrans = NaN(1, nSigma);
basinFraction = NaN(1, nSigma);

for ss = 1:nSigma
    b = b0 + sigmaVals(ss) * deltaB;

    [rho, theta, Omega, residual, ok] = solve_locked_state( ...
        A, b, epsVal, omega, rhoGuess, thetaGuess, OmegaGuess, ...
        solverTol, resTol, maxIter);

    if ~ok || residual > resTol
        continue
    end

    % Continuation: use the previous solution as the next initial guess.
    rhoGuess = rho;
    thetaGuess = theta;
    OmegaGuess = Omega;

    J = jacobian_pa(rho, theta, A, b, epsVal);
    ev = eig(J);

    % Remove the neutral global phase-shift eigenvalue.
    [~, idxNeutral] = min(abs(ev));
    ev(idxNeutral) = [];

    lambdaTrans(ss) = max(real(ev));
    basinFraction(ss) = estimate_basin( ...
        A, b, epsVal, omega, Omega, rho, theta, ...
        nIC, tFinal, ampTol, phaseTol);
end
end

function basin = estimate_basin( ...
    A, b, epsVal, omega, Omega, rho, theta, ...
    nIC, tFinal, ampTol, phaseTol)
% Estimate sizes of the basin of attraction by sampling random initial
% conditions

N = size(A, 1);
success = false(nIC, 1);

rhs = @(t, x) rhs_rotating(x, A, b, epsVal, omega, Omega);
opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);

for q = 1:nIC
    r0 = 1 + 2 * rand(N, 1);
    th0 = -pi + 2*pi * rand(N, 1);

    try
        [~, X] = ode45(rhs, [0 tFinal], [r0; th0], opts);
    catch
        continue
    end

    xEnd = X(end, :).';
    rEnd = xEnd(1:N);
    thEnd = wrap_pi(xEnd(N+1:end));

    ampErr = norm(rEnd - rho) / sqrt(N);
    phaseErr = phase_pattern_error(thEnd, theta);

    success(q) = ampErr < ampTol && phaseErr < phaseTol;
end

basin = mean(success);
end

function [rho, theta, Omega, residual, ok] = solve_locked_state( ...
    A, b, epsVal, omega, rhoGuess, thetaGuess, OmegaGuess, ...
    solverTol, resTol, maxIter)
% Finds stationary solution

N = size(A, 1);

thetaGuess = wrap_pi(thetaGuess(:) - thetaGuess(1));
y0 = [rhoGuess(:); thetaGuess(2:end); OmegaGuess];

F = @(y) locked_equations(y, A, b, epsVal, omega);

opts = optimoptions('fsolve', ...
    'Display', 'off', ...
    'FunctionTolerance', solverTol, ...
    'StepTolerance', solverTol, ...
    'OptimalityTolerance', solverTol, ...
    'MaxIterations', maxIter, ...
    'MaxFunctionEvaluations', max(5000, 1000*numel(y0)));

[y, fval, exitflag] = fsolve(F, y0, opts);

rho = y(1:N);
theta = [0; y(N+1:2*N-1)];
theta = wrap_pi(theta - theta(1));
Omega = y(end);

residual = norm(fval) / sqrt(numel(fval));
ok = exitflag > 0 && residual < resTol && all(rho > 0) && all(isfinite(y));
end

function F = locked_equations(y, A, b, epsVal, omega)
N = size(A, 1);

rho = y(1:N);
theta = [0; y(N+1:2*N-1)];
Omega = y(end);

D = theta.' - theta;
C = sum(A .* cos(D), 2);
S = sum(A .* sin(D), 2);

Fr = b .* rho .* (1 - rho) + epsVal .* rho .* C;
Ftheta = omega + rho - 1 + rho .* S - Omega;

F = [Fr; Ftheta];
end

function dx = rhs_rotating(x, A, b, epsVal, omega, Omega)
% ODEs of the coupled phase-amplitude oscillators

N = size(A, 1);

r = x(1:N);
theta = x(N+1:end);

D = theta.' - theta;
C = sum(A .* cos(D), 2);
S = sum(A .* sin(D), 2);

rdot = b .* r .* (1 - r) + epsVal .* r .* C;
thetadot = omega + r - 1 + r .* S - Omega;

dx = [rdot; thetadot];
end

function J = jacobian_pa(r, theta, A, b, epsVal)
% Calculates Jacobian matrix

N = numel(r);

D = theta.' - theta;
SinD = sin(D);
CosD = cos(D);

C = sum(A .* CosD, 2);
S = sum(A .* SinD, 2);

Jrr = diag(b .* (1 - 2*r) + epsVal .* C);
Jrth = zeros(N);
Jthr = diag(1 + S);
Jthth = zeros(N);

for i = 1:N
    for j = 1:N
        if i == j
            Jrth(i,j) = epsVal * r(i) * sum(A(i,:) .* SinD(i,:));
            Jthth(i,j) = -r(i) * sum(A(i,:) .* CosD(i,:));
        else
            Jrth(i,j) = -epsVal * r(i) * A(i,j) * SinD(i,j);
            Jthth(i,j) = r(i) * A(i,j) * CosD(i,j);
        end
    end
end

J = [Jrr, Jrth; Jthr, Jthth];
end

function err = phase_pattern_error(thetaObserved, thetaTarget)
d = wrap_pi(thetaObserved(:) - thetaTarget(:));
shift = angle(mean(exp(1i*d)));
err = norm(wrap_pi(d - shift)) / sqrt(numel(d));
end

function y = wrap_pi(x)
y = mod(x + pi, 2*pi) - pi;
end

function shaded_curve(x, y, err)
x = x(:).';
y = smoothdata(y(:).', 'movmean', 3);
err = smoothdata(err(:).', 'movmean', 3);

hold on
fill([x fliplr(x)], ...
    [y - err fliplr(y + err)], ...
    [0.7 0.7 0.7], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.35);

plot(x, y, 'LineWidth', 1.5);
hold off
end

function A = directed_ring_circulant(N, W1, W2, K)
% Generates a circulant network.    
A = zeros(N);

for i = 1:N
    for k = 1:K
        j1 = mod(i - 1 + k, N) + 1;
        j2 = mod(i - 1 - k, N) + 1;

        A(i, j1) = W1;
        A(i, j2) = W2;
    end
end
end