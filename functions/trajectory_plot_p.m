function trajectory_plot_p(var,mode,mapgeneID,meanG,semG,mapTFID,meanTF,semTF)

Acol = [0 121 185]/255;
Bcol = [165 225 127]/255;

if strcmp(mode, 'gene')
    
    % expression over trajectories 
    %figure()
    box on
    x = 10:20000;
    y = meanG{6}(10:20000,mapgeneID(var))';
    sd = semG{6}(10:20000,mapgeneID(var))';
    plot(x, y, 'color', Acol)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], Acol, 'EdgeColor','none'), alpha(0.5)
    axis tight

    hold on 
    y = meanG{5}(10:20000,mapgeneID(var))';
    sd = semG{5}(10:20000,mapgeneID(var))';
    plot(x, y, 'color', Bcol)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], Bcol, 'EdgeColor','none'), alpha(0.5)
    axis tight
    title({var}, 'Interpreter', 'none');
    xlabel('Simulated time');
    ylabel('Expression');
    hold off 

elseif strcmp(mode, 'tf')
    
    % activity over trajectories 
    %figure()
    box on
    x = 10:20000;
    y = meanTF{6}(10:20000,mapTFID(var))';
    sd = semTF{6}(10:20000,mapTFID(var))';
    plot(x, y, 'color', Acol)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], Acol, 'EdgeColor','none'), alpha(0.5)
    axis tight

    hold on 
    y = meanTF{5}(10:20000,mapTFID(var))';
    sd = semTF{5}(10:20000,mapTFID(var))';
    plot(x, y, 'color', Bcol)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], Bcol, 'EdgeColor','none'), alpha(0.5)
    axis tight
    title({var}, 'Interpreter', 'none');
    xlabel('Simulated time');
    ylabel('Activity');
    hold off 
    
else
    
    figure()
    
    
end
    
    % 6th panel cross-correlation
%     [Ca,laga] = xcorr((meanTF{6}(:,mapTFID(motif))-mean(meanTF{6}(:,mapTFID(motif)))),...
%         (meanG{6}(:,mapgeneID(tf{k}))-mean(meanG{6}(:,mapgeneID(tf{k})))),'coeff');
%     
%     [Cb,lagb] = xcorr((meanTF{5}(:,mapTFID(motif))-mean(meanTF{5}(:,mapTFID(motif)))),...
%         (meanG{5}(:,mapgeneID(tf{k}))-mean(meanG{5}(:,mapgeneID(tf{k})))),'coeff');
%     
%     subplot(2,3,6)
%     plot(laga,Ca,'color',Acol,'LineWidth',2),hold on
%     axis tight
%     hold on
%     plot(lagb,Cb,'color',Bcol,'LineWidth',2),hold on
%     str = strcat({'Expression vs. activity'});
%     title(str);
%     xlabel('Time lag');
%     ylabel('pearson correlation');
%     axis tight
%     
  
end 