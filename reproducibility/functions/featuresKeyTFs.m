function featuresKeyTFs(motif,E2,G,A,mapgeneID,meanG,semG,ffG,mapTFID,meanTF,semTF,ffTF)

tf = split(motif,"_");
Acol = [0 121 185]/255;
Bcol = [165 225 127]/255;

k = 1; 

while k <= length(tf)
    
    if isKey(mapgeneID,tf{k}) == 1
    
    %figure('position',[0 0 700 300])
    figure('position',[0 0 300 700])
    
    % 1st panel umap with activity profile
    subplot(3,2,1)
    %subplot(2,3,1)
    scatter(E2(:,1),E2(:,2),10,A(:,mapTFID(motif)),'filled')
    str = strcat(motif, {' activity'});
    %title(str, 'Interpreter', 'none');
    colormap(redbluecmap)
    axis off
    
    % 2nd panel activity over trajectories 
    subplot(3,2,3)
    %subplot(2,3,2)
    box on
    x = 10:1500;
    y = meanTF{6}(10:1500,mapTFID(motif))';
    sd = semTF{6}(10:1500,mapTFID(motif))';
    plot(x, y, 'color', Acol)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], Acol, 'EdgeColor','none'), alpha(0.5)
    %ylim([-0.15 0.15])
    axis tight

    hold on 
    y = meanTF{5}(10:1500,mapTFID(motif))';
    sd = semTF{5}(10:1500,mapTFID(motif))';
    plot(x, y, 'color', Bcol)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], Bcol, 'EdgeColor','none'), alpha(0.5)
    %ylim([-0.05 0.06])
    axis tight
    str = strcat({'Profile over stochastic trajectories'});
    %title(str);
    %xlabel('Simulated time');
    %ylabel({motif, 'activity'}, 'Interpreter', 'none');
    hold off 
    
    % 3rd panel  fano factor activities 
%     subplot(2,3,3)
%     box on
%     hold on
%     plot(log(ffTF{6}(10:1500,mapTFID(motif))),'color', Acol,'LineWidth',2);
%     plot(log(ffTF{5}(10:1500,mapTFID(motif))),'color', Bcol,'LineWidth',2);
%     str = strcat(motif, {' fano factor'});
%     title(str, 'Interpreter', 'none');
%     xlabel('Simulated time');
%     ylabel('Log(Fano factor)');
%     hold off
    
    % 4th panel umap with the expression profile 
    subplot(3,2,2)
    %subplot(2,3,4)
    scatter(E2(:,1),E2(:,2),10,G(:,mapgeneID(tf{k})),'filled')
    str = strcat(tf{k}, {'  expression'});
    %title(str);
    colormap(redbluecmap)
    axis off
    
    % 5ft panel expression over trajectories 
    subplot(3,2,4)
    %subplot(2,3,5)
    box on
    x = 10:1500;
    y = meanG{6}(10:1500,mapgeneID(tf{k}))';
    sd = semG{6}(10:1500,mapgeneID(tf{k}))';
    plot(x, y, 'color', Acol)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], Acol, 'EdgeColor','none'), alpha(0.5)
    %ylim([-0.15 0.15])
    axis tight

    hold on 
    y = meanG{5}(10:1500,mapgeneID(tf{k}))';
    sd = semG{5}(10:1500,mapgeneID(tf{k}))';
    plot(x, y, 'color', Bcol)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], Bcol, 'EdgeColor','none'), alpha(0.5)
    %ylim([-0.05 0.06])
    axis tight
    str = strcat({'Profile over stochastic trajectories'});
    %title(str);
    %xlabel('Simulated time');
    %ylabel({tf{k},'expression'});
    hold off 

    
    % 6th panel cross-correlation
    [Ca,laga] = xcorr((meanTF{6}(:,mapTFID(motif))-mean(meanTF{6}(:,mapTFID(motif)))),...
        (meanG{6}(:,mapgeneID(tf{k}))-mean(meanG{6}(:,mapgeneID(tf{k})))),'coeff');
    
    [Cb,lagb] = xcorr((meanTF{5}(:,mapTFID(motif))-mean(meanTF{5}(:,mapTFID(motif)))),...
        (meanG{5}(:,mapgeneID(tf{k}))-mean(meanG{5}(:,mapgeneID(tf{k})))),'coeff');
    
    subplot(3,2,[5 6])
    %subplot(2,3,6)
    plot(laga,Ca,'color',Acol,'LineWidth',2),hold on
    axis tight
    hold on
    plot(lagb,Cb,'color',Bcol,'LineWidth',2),hold on
    str = strcat({'Expression vs. activity'});
    %title(str);
    %xlabel('Time lag');
    %ylabel('pearson correlation');
    axis tight
    
    
    end
    

    k = k + 1;

end
end 