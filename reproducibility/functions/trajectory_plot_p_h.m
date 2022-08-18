function trajectory_plot_p_h(var,mode,mapgeneID,meanG,semG,mapTFID,meanTF,semTF)

Acol = [53 122 197]/255;
Bcol = [165 225 127]/255;
ECcol = [238 132 124]/255;

if strcmp(mode, 'gene')
    
    % expression over trajectories 
    %figure()
    box on
    x = 10:1000;
    y = meanG{7}(10:1000,mapgeneID(var))';
    sd = semG{7}(10:1000,mapgeneID(var))';
    plot(x, y, 'color', Acol)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], Acol, 'EdgeColor','none'), alpha(0.5)
    axis tight

    hold on 
    y = meanG{8}(10:1000,mapgeneID(var))';
    sd = semG{8}(10:1000,mapgeneID(var))';
    plot(x, y, 'color', Bcol)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], Bcol, 'EdgeColor','none'), alpha(0.5)
    axis tight
    
    hold on 
    y = meanG{9}(10:1000,mapgeneID(var))';
    sd = semG{9}(10:1000,mapgeneID(var))';
    plot(x, y, 'color', ECcol)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], ECcol, 'EdgeColor','none'), alpha(0.5)
    axis tight
    title({var}, 'Interpreter', 'none');
    xlabel('Simulated time');
    ylabel('Expression');
    hold off 

elseif strcmp(mode, 'tf')
    
    % activity over trajectories 
    %figure()
    box on
    x = 10:1000;
    y = meanTF{7}(10:1000,mapTFID(var))';
    sd = semTF{7}(10:1000,mapTFID(var))';
    plot(x, y, 'color', Acol)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], Acol, 'EdgeColor','none'), alpha(0.5)
    axis tight

    hold on 
    y = meanTF{8}(10:1000,mapTFID(var))';
    sd = semTF{8}(10:1000,mapTFID(var))';
    plot(x, y, 'color', Bcol)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], Bcol, 'EdgeColor','none'), alpha(0.5)
    axis tight
    
    hold on 
    y = meanTF{9}(10:1000,mapTFID(var))';
    sd = semTF{9}(10:1000,mapTFID(var))';
    plot(x, y, 'color', ECcol)
    hold on 
    patch([x fliplr(x)], [y-sd  fliplr(y+sd)], ECcol, 'EdgeColor','none'), alpha(0.5)
    axis tight
    title({var}, 'Interpreter', 'none');
    xlabel('Simulated time');
    ylabel('Activity');
    hold off 
    
else
    
    figure()
    
    
end
  
end 