%% Calculates the equilibrium of FitzHugh-Nagumo neuron models
%   x       -   equilibrium state
%   N       -   network size
%   a,b,I,K -   intrinsic parameters
%   A       -   adjacency matrix

function F = fhn_system(x, N, a, epsilon, b, I, K, A)
    v = x(1:N);
    w = x(N+1:end);

    Fv = zeros(N,1);
    Fw = zeros(N,1);
    
    for i = 1:N
        coupling = K * sum(A(i,:) .* (v' - v(i)));
        Fv(i) = v(i) - (v(i)^3)/3 - w(i) + I + coupling;
        Fw(i) = epsilon * (v(i) + a - b(i)*w(i));
    end

    F = [Fv; Fw];
end