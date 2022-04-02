function plotdiffmap(X,A,colcluster)

    figure
    for i = 1:length(A)
        for j = (i+1):length(A)
            if (A(i,j) > 0) 
                % gray = [208/255 206/255 206/255];
                gray = [0.2 0.2 0.2];
                line([X(1,i) X(1,j)], [X(2,i) X(2,j)],'Color',gray,'LineWidth',0.5);hold on; 
            end
        end        
    end
%     scatter(X(1,:),X(2,:),20,'MarkerFaceColor',colcluster,...
%         'MarkerEdgeColor',colcluster,'MarkerFaceAlpha',.2,...
%         'MarkerEdgeAlpha',.2);
    scatter(X(1,:),X(2,:),20,colcluster,'filled');
    axis off  
end