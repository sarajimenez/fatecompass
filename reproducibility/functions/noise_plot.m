function noise_plot(var,mode,mapgeneID,meanG,cvG,ffG,semG,mapTFID,meanTF,cvTF,ffTF,semTF)

Acol = [0 121 185]/255;
Bcol = [165 225 127]/255;

if strcmp(mode, 'gene')

    figure('position',[0 0 650 150]) %this works for (1,3,x)
    
    %figure('position',[0 0 650 500])
    
    % 1st panel Expression profile over stochastic trajectories
    subplot(1,3,1)
    %subplot(2,2,1)
    box on
    x = 10:1500;
    y = meanG{6}(10:1500,mapgeneID(var))';
    sd = semG{6}(10:1500,mapgeneID(var))';
    plot(x, y, 'color', Acol)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], Acol, 'EdgeColor','none'), alpha(0.5)
    axis tight

    hold on 
    y = meanG{5}(10:1500,mapgeneID(var))';
    sd = semG{5}(10:1500,mapgeneID(var))';
    plot(x, y, 'color', Bcol)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], Bcol, 'EdgeColor','none'), alpha(0.5)
    axis tight
    xlabel('Simulated time');
    ylabel({var,'expression'});
    hold off 
    
    yl = ylim;
    
    % 2nd panel  fano factor over trajectories 
%     subplot(2,2,2)
%     box on
%     hold on
%     plot(log(ffG{6}(10:1500,mapgeneID(var))),'color', Acol,'LineWidth',2);
%     plot(log(ffG{5}(10:1500,mapgeneID(var))),'color', Bcol,'LineWidth',2);
%     xlabel('Simulated time');
%     ylabel('Log(Fano factor)');
%     hold off
 
    % 3rd panel noise over alpha trajectories 
    subplot(1,3,2)
    %subplot(2,2,3)
    box on
    x = 10:1500;
    scatter(x,meanG{6}(10:1500,mapgeneID(var)),...
    cvG{6}(10:1500,mapgeneID(var)).*10,...
    ffG{6}(10:1500,mapgeneID(var)),'o','filled');
    ylim(yl);
    title({var, 'noise over Alpha trajectories'}, 'Interpreter', 'none');
    xlabel('Simulated time');
    ylabel('Expression');
    alpha(.5)
    colormap(redbluecmap)
    c = colorbar;
    c.Label.String = 'Fano factor';
    box off
    
    % 4th panel noise over beta trajectories 
    subplot(1,3,3)
    %subplot(2,2,4)
    box on
    x = 10:1500;
    scatter(x,meanG{5}(10:1500,mapgeneID(var)),...
    cvG{5}(10:1500,mapgeneID(var)).*10,...
    ffG{5}(10:1500,mapgeneID(var)),'o','filled');
    ylim(yl);
    title({var, 'noise over Beta trajectories'}, 'Interpreter', 'none');
    xlabel('Simulated time');
    ylabel('Expression');
    alpha(.5)
    colormap(redbluecmap)
    c = colorbar;
    c.Label.String = 'Fano factor';
    box off

elseif strcmp(mode, 'tf')

    figure('position',[0 0 650 150])
    
    %figure('position',[0 0 650 500])
    
    % 1st panel Activity profile over stochastic trajectories
    subplot(1,3,1)
    %subplot(2,2,1)
    box on
    x = 10:1500;
    y = meanTF{6}(10:1500,mapTFID(var))';
    sd = semTF{6}(10:1500,mapTFID(var))';
    plot(x, y, 'color', Acol)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], Acol, 'EdgeColor','none'), alpha(0.5)
    axis tight

    hold on 
    y = meanTF{5}(10:1500,mapTFID(var))';
    sd = semTF{5}(10:1500,mapTFID(var))';
    plot(x, y, 'color', Bcol)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], Bcol, 'EdgeColor','none'), alpha(0.5)
    axis tight
    xlabel('Simulated time');
    ylabel({var, 'activity'}, 'Interpreter', 'none');
    hold off 
    
    yl = ylim;
    
    % 2nd panel  fano factor over trajectories 
%     subplot(2,2,2)
%     box on
%     hold on
%     plot(log(ffTF{6}(10:1500,mapTFID(var))),'color', Acol,'LineWidth',2);
%     plot(log(ffTF{5}(10:1500,mapTFID(var))),'color', Bcol,'LineWidth',2);
%     xlabel('Simulated time');
%     ylabel('Log(Fano factor)');
%     hold off
%     
    %I did log of cv and ff to have nicer scales, I'd to take the abs value
    %of cv bc size can only be positive numbers. 
    
    % 3rd panel noise over alpha trajectories 
    subplot(1,3,2)
    %subplot(2,2,3)
    box on
    x = 10:1500;
    scatter(x,meanTF{6}(10:1500,mapTFID(var)),...
    cvTF{6}(10:1500,mapTFID(var))*200,...
    cvTF{6}(10:1500,mapTFID(var)),'o','filled');
    ylim(yl);
    title({var, 'noise over Alpha trajectories'}, 'Interpreter', 'none');
    xlabel('Simulated time');
    ylabel('Activity');
    alpha(.5)
    colormap(redbluecmap)
    c = colorbar;
    c.Label.String = 'C.V.';
    box off
    
    % 4th panel noise over beta trajectories 
    subplot(1,3,3)
    %subplot(2,2,4)
    box on
    x = 10:1500;
    scatter(x,meanTF{5}(10:1500,mapTFID(var)),...
    cvTF{5}(10:1500,mapTFID(var))*200,...
    cvTF{5}(10:1500,mapTFID(var)),'o','filled');
    ylim(yl);
    title({var, 'noise over Beta trajectories'}, 'Interpreter', 'none');
    xlabel('Simulated time');
    ylabel('Activity');
    alpha(.5)
    colormap(redbluecmap)
    c = colorbar;
    c.Label.String = 'C.V.';
    box off
    
     
else
    
    figure()
    
    
end

end 