function n = average_number_connections(Dv,D,cutoff)
    Pv = exp(-(Dv.^2)/D);
    Pv = Pv./sum(Pv);
    n = mean(sum(Pv>cutoff));   
end