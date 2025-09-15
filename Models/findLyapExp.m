function [LyapExp, xeq] = findLyapExp(system,param,N,p,zeroLE)
%% Calculates the largest Lyapunov exponent for different dynamical systems
%
% LyapExp   -  largest (transversal) Lyapunov exponent
% xeq       -  equilibrium point (calculated only for certain systems)
% system    -  dynamical model under consideration
% param     -  tuple of parameters of the considered model
% N         -  network size
% p         -  optimization parameter
% zeroLE    -  0 if the system does not have identically null exponents,
%              1 otherwise

% Jacobian matrix
switch system
    case 'van-der-pol'
        % Parameters
        yeq = param.yeq;
        xeq = param.xeq;
        L = param.L;
        
        % Jacobian matrix
        J = [zeros(N,N), eye(N); ...
            - L - eye(N), - diag(p)];

    case 'power-grid'
        % Parameters
        B = diag(p);
        L = param.L;
        xeq = zeros(N,1);

        % Jacobian
        J11 = zeros(N,N);
        J12 = eye(N);
        J21 = - L;
        J22 = - B;
        J = [J11 J12; J21 J22];

    case 'neuron'
            
        % Parameters
        b = p;
        a = param.a;
        eps = param.eps;
        K = param.K;
        L = param.L;
        A = param.A;
        
        % Find equilibrium
        x0 = rand(2*N,1);
        options = optimoptions('fsolve', 'TolFun', 1e-2);
        xeq = fsolve(@(x) fhn_system(x, N, a, eps, b, 0, K, A), x0, options);
        veq = xeq(1:N,1);

        % Jacobian
        Jvv = (eye(N) - diag(veq).^2) - K*L;
        Jvw = - eye(N);
        Jwv = eps*eye(N);
        Jww = - eps*diag(b);

        J = [Jvv Jvw; Jwv Jww];

    case 'metamaterial'
        % Parameters
        B = diag(p);
        Minv = diag(1./param.m);
        k = param.k;
        L = param.L;
        xeq = zeros(N,1);

        % Jacobian
        J11 = zeros(N,N);
        J12 = eye(N);
        J21 = - Minv * (L + k*eye(N));
        J22 = - Minv * B;
        J = [J11 J12; J21 J22];
        
    case 'multi-agent'
        % Parameters
        B = diag(p);
        k = param.k;
        gamma = param.gamma;
        L = param.L;
        xeq = zeros(N,1);

        % Jacobian
        J11 = zeros(N,N);
        J12 = eye(N);
        J21 = - (k*L + B);
        J22 = - gamma*(k*L + B);

        J = [J11 J12; J21 J22];

    case 'josephson-junction'

        % Parameters
        R = p; 
        Cinv = diag(1./param.C);
        Ic = param.Ic;
        Ib = param.Ib;
        L = param.L;

        % Find equilibrium
        xeq = asin(Ib./Ic);

        % Jacobian matrix
        J11 = zeros(N,N);
        J12 = eye(N);
        J21 = - Cinv*(diag(Ic.*cos(xeq)));
        J22 = - Cinv *(L + diag(R));

        J = [J11 J12; J21 J22];

    case 'phase-amplitude'

        % Parameters
        b = p;
        eps = param.eps;
        gamma = param.gamma;
        A = param.A;
        L = param.L;
        xeq = zeros(2*N,1);

        % Jacobian
        Jrr = - diag(b);
        Jrt = - eps*A;
        Jtr = eye(N);
        Jtt = - gamma*L;
        J = [Jrr Jrt; Jtr Jtt];
        
    case 'first-order'

        % Parameters
        b = p;
        L = param.L;
        k = param.k;
        xeq = zeros(N,1);

        % Jacobian
        J = -k*L - diag(b);

end


% Eigenvalues around equilibrium point
eigJ = eig(J);
real_eig = real(eigJ);
if zeroLE == 1
    threshold = 1e-10;
    [~, idx] = min(abs(real_eig));
    if abs(real_eig(idx)) < threshold  % check if it is within the threshold (close to zero)
        real_eig(idx) = [];            % remove the eigenvalue closest to zero
    end
end
LyapExp = max(real_eig);

end