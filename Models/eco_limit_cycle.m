function results = main()

    % ---------------------------------------------------------------------
    % User input
    % ---------------------------------------------------------------------
    p_edge = input('Choose edge probability: ');
    p_sign = input('Choose sign probability: ');

    parameter = lower(strtrim(input( ...
        'Choose control parameter (network; connectivity; sign; growth; none): ', ...
        's')));

    % ---------------------------------------------------------------------
    % Default parameters
    % ---------------------------------------------------------------------
    N = 10;

    Nhet = 400;

    sigma = linspace(0,1,Nhet);

    growth_het = 10;

    Nmc = 100;

    % ---------------------------------------------------------------------
    % Control parameter range
    % ---------------------------------------------------------------------
    switch parameter

        case 'network'
            param_range = [2 5 10 20 50 100 200 500 1000];

        case 'connectivity'
            param_range = 1;

        case 'sign'
            param_range = 0:0.05:1;

        case 'growth'
            param_range = 0:1:20;

        case 'none'
            param_range = 1;

        otherwise
            error('Unknown control parameter.');

    end

    % ---------------------------------------------------------------------
    % Stability analysis
    % ---------------------------------------------------------------------
    LyapExp = cell(length(param_range),1);

    feasible = cell(length(param_range),1);

    for j = 1:length(param_range)

        p = param_range(j);

        fprintf('Control parameter (%s): %g\n', parameter, p);

        % Defaults
        N_j = N;
        p_edge_j = p_edge;
        p_sign_j = p_sign;
        growth_het_j = growth_het;

        % Update controlled parameter
        switch parameter

            case 'network'
                N_j = round(p);

            case 'connectivity'
                p_edge_j = p;

            case 'sign'
                p_sign_j = p;

            case 'growth'
                growth_het_j = p;

        end

        % Monte Carlo realizations
        LyapExp_run = zeros(Nmc,Nhet);

        feasible_run = zeros(Nmc,Nhet);

        parfor mc = 1:Nmc

            seed = 10000*j + mc;

            [lyap, feas] = mc_realization( ...
                seed, N_j, p_edge_j, p_sign_j, ...
                growth_het_j, sigma);

            LyapExp_run(mc,:) = lyap;

            feasible_run(mc,:) = feas;

        end

        LyapExp{j} = LyapExp_run;

        feasible{j} = feasible_run;

    end

    % ---------------------------------------------------------------------
    % Performance analysis
    % ---------------------------------------------------------------------
    stabilityimprovement = cell(length(param_range),1);

    stabilizedfraction = cell(length(param_range),1);

    for j = 1:length(param_range)

        baseline = LyapExp{j}(:,1);

        stabilityimprovement{j} = LyapExp{j} - baseline;

        feas = feasible{j};

        lam = LyapExp{j};

        improved = (lam .* feas) < baseline;

        stabilizedfraction{j} = ...
            nansum(improved,1) ./ nansum(~isnan(feas),1);

    end

    % ---------------------------------------------------------------------
    % Plot
    % ---------------------------------------------------------------------
    Nsamps = 10;

    j = 1;

    figure('Position',[100 100 1200 900])

    subplot(2,2,1)

    plot(sigma, LyapExp{j}(1:Nsamps,:)','LineWidth',2)

    xlabel('heterogeneity \sigma')

    ylabel('\lambda_{max}')

    set(gca,'FontSize',14)

    subplot(2,2,2)

    plot(sigma, nanmean(stabilityimprovement{j},1), ...
        'LineWidth',2)

    xlabel('heterogeneity \sigma')

    ylabel('relative improvement')

    set(gca,'FontSize',14)

    subplot(2,2,3)

    plot(sigma, stabilizedfraction{j}, ...
        'LineWidth',2)

    xlabel('heterogeneity \sigma')

    ylabel('stabilized fraction')

    set(gca,'FontSize',14)

    subplot(2,2,4)

    axis off

    % ---------------------------------------------------------------------
    % Output
    % ---------------------------------------------------------------------
    results.sigma = sigma;

    results.param_range = param_range;

    results.LyapExp = LyapExp;

    results.feasible = feasible;

    results.stabilityimprovement = stabilityimprovement;

    results.stabilizedfraction = stabilizedfraction;

end

% =========================================================================
% Monte Carlo realization
% =========================================================================
function [lyap, feas] = mc_realization( ...
    seed, N, p_edge, p_sign, growth_het, sigma)

    rng(seed);

    Nhet = length(sigma);

    lyap = zeros(1,Nhet);

    feas = ones(1,Nhet);

    % ---------------------------------------------------------------------
    % Adjacency matrix
    % ---------------------------------------------------------------------
    A = econetwork(N,p_edge,p_sign);

    d_A = rand(N,N);

    for i = 1:N
        for j = 1:N

            if j ~= mod(i,N)+1
                d_A(i,j) = 0;
            end

        end
    end

    % ---------------------------------------------------------------------
    % Heterogeneity loop
    % ---------------------------------------------------------------------
    for k = 1:Nhet

        A_fin = A + sigma(k)/N * d_A;

        x0 = (1:N)'/N;

        x0 = x0/sum(x0);

        lyapvals = largest_lyapunov( ...
            A_fin, x0, ...
            2000, 20000, 1, 2);

        lyap(k) = lyapvals(2);

        if k == 1
            disp(lyap(k))
        end

    end

end

% =========================================================================
% Largest Lyapunov exponent (Benettin algorithm)
% =========================================================================
function lyap = largest_lyapunov( ...
    A, x0, t_transient, t_total, renorm_dt, k)

    N = size(A,1);

    % ---------------------------------------------------------------------
    % Burn transient
    % ---------------------------------------------------------------------
    opts = odeset('RelTol',1e-6,'AbsTol',1e-7);

    [~,X] = ode45(@(t,x) replicator_rhs(t,x,A), ...
        [0 t_transient], x0, opts);

    x = X(end,:)';

    % ---------------------------------------------------------------------
    % Initial tangent vectors
    % ---------------------------------------------------------------------
    [Q,~] = qr(randn(N,k),0);

    % ---------------------------------------------------------------------
    % Benettin loop
    % ---------------------------------------------------------------------
    n_steps = floor(t_total/renorm_dt);

    log_sums = zeros(k,1);

    for step = 1:n_steps

        y0 = [x; Q(:)];

        [~,Y] = ode45(@(t,y) variational_rhs(t,y,A,k), ...
            [0 renorm_dt], y0, opts);

        y = Y(end,:)';

        x = y(1:N);

        Q = reshape(y(N+1:end),N,k);

        [Q,R] = qr(Q,0);

        log_sums = log_sums + log(abs(diag(R)));

    end

    lyap = log_sums/(n_steps*renorm_dt);

end

% =========================================================================
% Variational equations
% =========================================================================
function dy = variational_rhs(~, y, A, k)

    N = size(A,1);

    x = y(1:N);

    Q = reshape(y(N+1:end),N,k);

    dx = replicator_rhs(0,x,A);

    J = jacobian_re(x,A);

    dQ = J*Q;

    dy = [dx; dQ(:)];

end

% =========================================================================
% Replicator RHS
% =========================================================================
function dx = replicator_rhs(~, x, A)

    Ax = A*x;

    avg_payoff = x' * Ax;

    dx = x .* (Ax - avg_payoff);

end

% =========================================================================
% Replicator Jacobian
% =========================================================================
function J = jacobian_re(x, A)

    Ax = A*x;

    avg_payoff = x' * Ax;

    J = diag(Ax - avg_payoff) + ...
        x .* (A - ((A + A')*x));

end

% =========================================================================
% Lotka-Volterra RHS
% =========================================================================
function dx = lotka_volterra_rhs(~, x, A, b)

    dx = x .* (b + A*x);

end

% =========================================================================
% Lotka-Volterra Jacobian
% =========================================================================
function J = jacobian_lv(x, A, b)

    J = diag(b + A*x) + x .* A;

end

% =========================================================================
% Ecological network
% =========================================================================
function A = econetwork(N, p_edge, p_sign)

    a = ones(N,1);

    A = zeros(N,N);

    for i = 1:N

        A(i,mod(i,N)+1) = a(i);

        A(i,i) = -0.5;

    end

end