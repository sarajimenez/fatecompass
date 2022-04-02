function c = colorscatter(vector,d)
    if (d == 2)
        c = linspace(1,10,length(vector));
    else 
        c = zeros(length(vector),3);
        map = colormap('parula');
        for i = 1:length(c)
            n = floor(length(map)*(i-1)/length(c))+1;
            c(i,:) = map(n,:);
        end
    end
end