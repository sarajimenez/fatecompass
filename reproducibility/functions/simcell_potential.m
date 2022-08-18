function states = simcell_potential(E,K,s0,numiter,numcells)
    %state = s0;
    %energy = E(state);
    states = zeros(numiter,numcells);
    
    for f = 1:numcells
        
        fprintf('%i\n',f);
        % state and energy should be here so they are reinicialized for
        % each simulated cell.
        state = s0;
        energy = E(state);
        states(1,f) = s0;
        
        for i = 1:numiter-1
            
            numnn = length(K{state});
            newstate = K{state}(randi(numnn));%K(state,randi(numnn));
            newenergy = E(newstate);
            
            if (exp(newenergy-energy) > rand())
                state = newstate;
                energy = newenergy; 
            end
            
            states(i+1,f) = state;
            
        end
    end
end
