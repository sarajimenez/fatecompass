function featuresKeyTFs_h(motif,E2,G,A,mapgeneID,meanG,semG,ffG,mapTFID,meanTF,semTF,ffTF)

tf = split(motif,"_");

EC_c    = [238 132 124]/255;
A_c     = [53 122 197]/255;
B_c     = [165 225 127]/255;

k = 1; 

while k <= length(tf)
    
    if isKey(mapgeneID,tf{k}) == 1
    
    %figure('position',[0 0 700 300])
    figure('position',[0 0 300 700])
    
    % 1st panel umap with activity profile
    subplot(3,2,1)
    scatter(E2(:,1),E2(:,2),10,A(:,mapTFID(motif)),'filled')
    str = strcat(motif, {' activity'});
    title(str, 'Interpreter', 'none');
    colormap(redbluecmap)
    axis off
    
    % 2nd panel activity over trajectories 
    subplot(3,2,3)
    box on
    x = 10:1000;
    y = meanTF{7}(10:1000,mapTFID(motif))';
    sd = semTF{7}(10:1000,mapTFID(motif))';
    plot(x, y, 'color', A_c)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], A_c, 'EdgeColor','none'), alpha(0.5)
    %ylim([-0.15 0.15])
    axis tight

    hold on 
    y = meanTF{8}(10:1000,mapTFID(motif))';
    sd = semTF{8}(10:1000,mapTFID(motif))';
    plot(x, y, 'color', B_c)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], B_c, 'EdgeColor','none'), alpha(0.5)
    %ylim([-0.05 0.06])
    axis tight
    str = strcat({'Profile over stochastic trajectories'});
    %title(str);
    %xlabel('Simulated time');
    %ylabel({motif, 'activity'}, 'Interpreter', 'none');
    hold off 
    
    
    % 4th panel umap with the expression profile 
    subplot(3,2,2)
    scatter(E2(:,1),E2(:,2),10,G(:,mapgeneID(tf{k})),'filled')
    str = strcat(tf{k}, {'  expression'});
    title(str);
    colormap(redbluecmap)
    axis off
    
    % 5ft panel expression over trajectories 
    subplot(3,2,4)
    box on
    x = 10:1000;
    y = meanG{7}(10:1000,mapgeneID(tf{k}))';
    sd = semG{7}(10:1000,mapgeneID(tf{k}))';
    plot(x, y, 'color', A_c)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], A_c, 'EdgeColor','none'), alpha(0.5)
    %ylim([-0.15 0.15])
    axis tight

    hold on 
    y = meanG{8}(10:1000,mapgeneID(tf{k}))';
    sd = semG{8}(10:1000,mapgeneID(tf{k}))';
    plot(x, y, 'color', B_c)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], B_c, 'EdgeColor','none'), alpha(0.5)
    %ylim([-0.05 0.06])
    axis tight
    str = strcat({'Profile over stochastic trajectories'});
    %title(str);
    %xlabel('Simulated time');
    %ylabel({tf{k},'expression'});
    hold off 

    
    % 6th panel cross-correlation
    [Ca,laga] = xcorr((meanTF{7}(:,mapTFID(motif))-mean(meanTF{7}(:,mapTFID(motif)))),...
        (meanG{7}(:,mapgeneID(tf{k}))-mean(meanG{7}(:,mapgeneID(tf{k})))),'coeff');
    
    [Cb,lagb] = xcorr((meanTF{8}(:,mapTFID(motif))-mean(meanTF{8}(:,mapTFID(motif)))),...
        (meanG{8}(:,mapgeneID(tf{k}))-mean(meanG{8}(:,mapgeneID(tf{k})))),'coeff');
    
    subplot(3,2,[5 6])
    plot(laga,Ca,'color',A_c,'LineWidth',2),hold on
    axis tight
    hold on
    plot(lagb,Cb,'color',B_c,'LineWidth',2),hold on
    str = strcat({'Expression vs. activity'});
    title(str);
    xlabel('Time lag');
    ylabel('pearson correlation');
    axis tight
    
    
    end
    

    k = k + 1;

end
end 