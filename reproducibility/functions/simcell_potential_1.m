function states = simcell_potential_1(E,K,s0,numiter,numcells)
    numnn = size(K,2);    
    states = zeros(numiter,numcells);
    for f = 1:numcells
        fprintf('%i\n',f);
        states(1,f) = s0;
        state = s0;
        energy = E(state);
        for i = 1:numiter-1
            newstate = K(state,randi(numnn));
            newenergy = E(newstate);
            
            if (exp(newenergy-energy) > rand())
                state = newstate;
                energy = newenergy; 
            end
            
            states(i+1,f) = state;
        end
    end
end
