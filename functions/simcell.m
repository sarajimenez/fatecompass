function states = simcell(W,N,s0,numiter,numcells)
    states = zeros(numiter,numcells);
    for f = 1:numcells
        fprintf('%i\n',f);
        states(1,f) = s0;
        for i = 1:numiter-1
            n = randp(W{states(i,f)});
            states(i+1,f) = N{states(i,f)}(n);
        end
    end
end

function i = randp(p)
    c = cumsum(p);
    r = c(end)*rand();
    i = 1;
    while (r > c(i))
        i = i + 1;
    end
end