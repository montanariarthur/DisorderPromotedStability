function results = main()

    % ---------------------------------------------------------------------
    % Main parameters
    % ---------------------------------------------------------------------
    N = 10;
    Nhet = 10;

    sigma = linspace(0,0.1,Nhet);

    Nmc = 10;

    b = 0.1;
    c = 14.0;
    eps = 0.221;

    % Ring Laplacian
    L = zeros(N,N);

    for i = 1:N
        L(i,i) = -2;
        L(i,mod(i,N)+1) = 1;
        L(i,mod(i-2,N)+1) = 1;
    end

    % ---------------------------------------------------------------------
    % Monte Carlo realizations
    % ---------------------------------------------------------------------
    LyapExp_run = zeros(Nmc,Nhet);

    parfor mc = 1:Nmc
        LyapExp_run(mc,:) = mc_realization(mc,N,sigma,b,c,L,eps);
    end

    % ---------------------------------------------------------------------
    % Performance analysis
    % ---------------------------------------------------------------------
    baseline = LyapExp_run(:,1);

    stabilityimprovement = LyapExp_run - baseline;

    improved = LyapExp_run < baseline;

    stabilizedfraction = sum(improved,1);

    % ---------------------------------------------------------------------
    % Plot
    % ---------------------------------------------------------------------
    Nsamps = min(10,Nmc);

    figure('Position',[100 100 1200 900])

    subplot(2,2,1)
    plot(sigma, LyapExp_run(1:Nsamps,:)','LineWidth',2)
    xlabel('\sigma')
    ylabel('\lambda_{max}')
    set(gca,'FontSize',14)

    subplot(2,2,2)
    plot(sigma, mean(stabilityimprovement,1),'LineWidth',2)
    xlabel('\sigma')
    ylabel('relative improvement')
    set(gca,'FontSize',14)

    subplot(2,2,3)
    plot(sigma, stabilizedfraction,'LineWidth',2)
    xlabel('\sigma')
    ylabel('stabilized fraction')
    set(gca,'FontSize',14)

    subplot(2,2,4)
    axis off

    % ---------------------------------------------------------------------
    % Output structure
    % ---------------------------------------------------------------------
    results.sigma = sigma;
    results.LyapExp = LyapExp_run;
    results.stabilityimprovement = stabilityimprovement;
    results.stabilizedfraction = stabilizedfraction;

end

% =========================================================================
% Monte Carlo realization
% =========================================================================
function lyap = mc_realization(seed,N,sigma,b,c,L,eps)

    rng(seed);

    Nhet = length(sigma);

    lyap = zeros(1,Nhet);

    a = 0.1*ones(N,1);

    d_A = rand(N,1);

    for k = 1:Nhet

        a_fin = a + sigma(k)/N * d_A;

        x0 = (1:(3*N))'/(3*N);
        x0 = x0/sum(x0);

        lyap_3 = largest_lyapunov( ...
            a_fin, N, b, c, L, eps, x0, ...
            2000, 20000, 1, 15);

        lyap(k) = lyap_3(2);

        disp(lyap_3')
        fprintf('\n');

    end

end

% =========================================================================
% Largest Lyapunov exponents (Benettin algorithm)
% =========================================================================
function lyap = largest_lyapunov( ...
    a, N, b, c, L, eps, x0, ...
    t_transient, t_total, renorm_dt, k)

    dt = 0.01;

    X = x0;

    % ---------------------------------------------------------------------
    % Burn transient
    % ---------------------------------------------------------------------
    ntrans = floor(t_transient/dt);

    for i = 1:ntrans
        X = rk4_Rossler_step(X,dt,a,b,c,N,L,eps);
    end

    % ---------------------------------------------------------------------
    % Initial orthonormal tangent vectors
    % ---------------------------------------------------------------------
    [Q,~] = qr(randn(3*N,k),0);

    n_steps = floor(t_total/renorm_dt);

    log_sums = zeros(k,1);

    % ---------------------------------------------------------------------
    % Benettin loop
    % ---------------------------------------------------------------------
    for step = 1:n_steps

        Y = [X; Q(:)];

        nsmall = floor(renorm_dt/dt);

        for j = 1:nsmall
            Y = rk4_variational_step(Y,dt,a,b,c,N,L,eps,k);
        end

        X = Y(1:3*N);

        Q = reshape(Y(3*N+1:end),3*N,k);

        [Q,R] = qr(Q,0);

        log_sums = log_sums + log(abs(diag(R)));

    end

    lyap = log_sums/(n_steps*renorm_dt);

end

% =========================================================================
% RK4 step for Rossler system
% =========================================================================
function Xnew = rk4_Rossler_step(X,dt,a,b,c,N,L,eps)

    k1 = Rossler_rhs(X,a,b,c,N,L,eps);

    k2 = Rossler_rhs(X + 0.5*dt*k1,a,b,c,N,L,eps);

    k3 = Rossler_rhs(X + 0.5*dt*k2,a,b,c,N,L,eps);

    k4 = Rossler_rhs(X + dt*k3,a,b,c,N,L,eps);

    Xnew = X + dt*(k1 + 2*k2 + 2*k3 + k4)/6;

end

% =========================================================================
% RK4 step for variational equations
% =========================================================================
function Ynew = rk4_variational_step(Y,dt,a,b,c,N,L,eps,k)

    k1 = variational_rhs(Y,a,b,c,N,L,eps,k);

    k2 = variational_rhs(Y + 0.5*dt*k1,a,b,c,N,L,eps,k);

    k3 = variational_rhs(Y + 0.5*dt*k2,a,b,c,N,L,eps,k);

    k4 = variational_rhs(Y + dt*k3,a,b,c,N,L,eps,k);

    Ynew = Y + dt*(k1 + 2*k2 + 2*k3 + k4)/6;

end

% =========================================================================
% Variational RHS
% =========================================================================
function dY = variational_rhs(Y,a,b,c,N,L,eps,k)

    X = Y(1:3*N);

    Q = reshape(Y(3*N+1:end),3*N,k);

    dX = Rossler_rhs(X,a,b,c,N,L,eps);

    J = jacobian_ross(X,a,c,eps,L,N);

    dQ = J*Q;

    dY = [dX; dQ(:)];

end

% =========================================================================
% Rossler RHS
% =========================================================================
function dX = Rossler_rhs(X,a,b,c,N,L,eps)

    x = X(1:N);

    y = X(N+1:2*N);

    z = X(2*N+1:3*N);

    dx = -y - z;

    dy = x + a.*y + eps*(L*y);

    dz = b + z.*x - c*z;

    dX = [dx; dy; dz];

end

% =========================================================================
% Jacobian
% =========================================================================
function J = jacobian_ross(X,a,c,eps,L,N)

    x = X(1:N);

    z = X(2*N+1:3*N);

    J = zeros(3*N,3*N);

    I = eye(N);

    % dx/dy and dx/dz
    J(1:N, N+1:2*N) = -I;
    J(1:N, 2*N+1:3*N) = -I;

    % dy/dx and dy/dy
    J(N+1:2*N,1:N) = I;
    J(N+1:2*N,N+1:2*N) = diag(a) + eps*L;

    % dz/dx and dz/dz
    for i = 1:N
        J(2*N+i,i) = z(i);
        J(2*N+i,2*N+i) = x(i) - c;
    end

end