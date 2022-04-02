function qave = absorbing_probabilities(indA,P,numabs)
    
    [numS,numA] = size(indA);
    numcells = length(P);
    Ia = eye(numA);
    Iq = eye(numcells-numA);
    indF = zeros(1,numA);
    qave = zeros(numcells,numA);
    
    for i = 1:numabs
        fprintf('%i\n',i)
        
        % Selection absorving/transient estates 
        for f = 1:numA
            indF(f) = indA(randi(numS),f);
        end
        indT = true(numcells,1);
        indT(indF) = 0;
        
        Q = P(indT,indT)';
        R = P(indF,indT)';
        q = (Iq-Q)\R;
        qexp = zeros(numcells,4);
        qexp(indT,:) = q;
        qexp(indF,:) = Ia;

        qave = qave + qexp;
    end
    qave = qave/numabs;
end