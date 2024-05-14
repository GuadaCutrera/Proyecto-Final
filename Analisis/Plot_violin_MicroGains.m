function Plot_violin_MicroGains(TL_mean,mogs_mean,mongs_mean,TL_group,mogs_group,mongs_group,titulo)

 
Matriz_Ganancias(:,1)=TL_group(:,end);
Matriz_Ganancias(:,2)=mogs_group(:,end);
Matriz_Ganancias(:,3)=mongs_group(:,end);

figure;set(gcf,'Color','white'); box OFF; hold on; sgtitle(titulo)
set(gcf,'Position',get(0,'ScreenSize'));

    %TL
    subplot(1,4,3)
    plot(TL_mean,'k.','MarkerSize', 15); hold on;
    for i=1:size(TL_group,1)
        s=scatter(1:size(TL_group,2),TL_group(i,:),'k','filled','MarkerEdgeAlpha',0.1);
        s.MarkerFaceAlpha = 0.1;
        hold on;
    end
    xlabel('Bloque','FontName','Arial','FontSize',12)
    ylabel('Suma de deltas (s)','FontName','Arial','FontSize',12); 
    ylim([-2 2])
    xlim([0.7 10.3])
    yline(0)
    title('Aprendizaaje Total')
    
    
    %MOGS
    subplot(1,4,1)
    plot(mogs_mean,'r.','MarkerSize', 15); hold on;
    for i=1:size(mogs_group,1)
        s=scatter(1:size(mogs_group,2),mogs_group(i,:),'r','filled','MarkerEdgeAlpha',0.1);
        s.MarkerFaceAlpha = 0.1;
        hold on;
    end
    xlabel('Bloque','FontName','Arial','FontSize',12)
    ylabel('Suma de deltas (s)','FontName','Arial','FontSize',12); 
    ylim([-2 2])
    xlim([0.7 10.3])
    title('MOGS')
    yline(0)
    
    
    %MONGS
    subplot(1,4,2)
    plot(mongs_mean,'b.','MarkerSize', 15); hold on;
    for i=1:size(mongs_group,1)
        s=scatter(1:size(mongs_group,2),mongs_group(i,:),'b','filled','MarkerEdgeAlpha',0.1);
        s.MarkerFaceAlpha = 0.1;
        hold on;
    end
    xlabel('Bloque','FontName','Arial','FontSize',12)
    ylabel('Suma de deltas (s)','FontName','Arial','FontSize',12); 
    ylim([-2 2])
    xlim([0.7 10.3])
    title('MONGS')
    yline(0)
    

    subplot(1,4,4)
    violin(Matriz_Ganancias,'xlabel',{'TL','MOGS','MONGS'},'facecolor',[.3 .3 .3; 1 0 0; 0 0 1],'edgecolor','k','mc','k','ylabel','Sum of deltas (s)','plotlegend',0,'medc',[]);
    plot(Matriz_Ganancias','k.','Markersize',15)
    yline(0)
    ylim([-2 2])
