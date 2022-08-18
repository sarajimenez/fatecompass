function i = randp(p)
    c = cumsum(p);
    r = c(end)*rand();
    i = 1;
    while (r > c(i))
        i = i + 1;
    end
end