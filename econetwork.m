function Adj = econetwork(N,p_edge,p_sign)
%% Generates a random network topology  
%  with the desired iteraction type for ecological networks.
%
%  Adj          -   weighted adjacency matrix
%  N            -   number of species
%  p_edge       -   connection probability
%  p_sign       -   sign probability, 
%                   from competitive (p=0) to mutualistic (p=1)

% Random network weights
Adj = abs(randn(N,N));       % only positive interactions thus far

for i = 1:N
    for j = 1:N

        % No self-edge
        if i == j                
            Adj(i,j) = 0;
        else
        
            % Edge probability (independent removal for directed edge)
            if rand > p_edge   
                Adj(i,j) = 0;
            else
                % Edge sign
                if rand > p_sign 
                    Adj(i,j) = - Adj(i,j);
                end
            end
        end
    end
end


end


