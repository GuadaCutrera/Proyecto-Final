function Plot_violin_MicroMicroGains(mogs_media_mean,mongs_media_mean,TL_media_mean,mogs_media_median,mongs_media_median,...
    mogs_mediana_mean,mongs_mediana_mean,TL_mediana_mean,mogs_mediana_median,mongs_mediana_median,...
    mogs_group_mediaB,mongs_group_mediaB,TL_group_mediaB,mogs_group_medianaB,mongs_group_medianaB,TL_group_medianaB,titulo)

Matriz_Ganancias_media(:,1)=mogs_group_mediaB(:,end);
Matriz_Ganancias_media(:,2)=mongs_group_mediaB(:,end);
Matriz_Ganancias_media(:,3)=TL_group_mediaB(:,end);

Matriz_Ganancias_mediana(:,1)=mogs_group_medianaB(:,end);
Matriz_Ganancias_mediana(:,2)=mongs_group_medianaB(:,end);
Matriz_Ganancias_mediana(:,3)=TL_group_medianaB(:,end);

figure;set(gcf,'Color','white'); box OFF; hold on; sgtitle(titulo)
set(gcf,'Position',get(0,'ScreenSize'));

%% media x bloque
%mogs
subplot(2,4,1)
for i=1:size(mogs_group_mediaB,1)
    s=scatter(1:size(mogs_group_mediaB,2),mogs_group_mediaB(i,:),'r','filled','MarkerEdgeAlpha',0.25);
    s.MarkerFaceAlpha = 0.25;
    hold on;
end
p_mean=plot(1:1:length(mogs_media_mean),mogs_media_mean,'k.','MarkerSize',12);
%p_median=plot(1:1:length(mogs_media_median),mogs_media_median,'ko','MarkerSize',6);
hold on;
xlabel('Bloque','FontName','Arial','FontSize',12);
ylabel(['Suma de deltas (s)'],'FontName','Arial','FontSize',12); 
title('MicroMOGS')
xlim([0 length(mogs_media_mean)+1]);
ylim([-2.5 2.5])
yline(0)
%legend([p_mean p_median],'mean','median')

%mongs
subplot(2,4,2)
for i=1:size(mongs_group_mediaB,1)
    s=scatter(1:size(mongs_group_mediaB,2),mongs_group_mediaB(i,:),'b','filled','MarkerEdgeAlpha',0.25);
    s.MarkerFaceAlpha = 0.25;
    hold on;
end
p_mean=plot(1:1:length(mongs_media_mean),mongs_media_mean,'k.','MarkerSize',12);
%p_median=plot(1:1:length(mongs_media_median),mongs_media_median,'ko','MarkerSize',6);
hold on;
xlabel('Bloque','FontName','Arial','FontSize',12);
ylabel(['Suma de deltas (s)'],'FontName','Arial','FontSize',12); 
title('MicroMONGS')
xlim([0 length(mongs_media_mean)+1]);
ylim([-2.5 2.5])
yline(0)
%legend([p_mean p_median],'mean','median')

subplot(2,4,3)
for i=1:size(TL_group_mediaB,1)
    s=scatter(1:size(TL_group_mediaB,2),TL_group_mediaB(i,:),'k','filled','MarkerEdgeAlpha',0.25);
    s.MarkerFaceAlpha = 0.25;
    hold on;
end
p_mean=plot(1:1:length(TL_media_mean),TL_media_mean,'k.','MarkerSize',12);
%p_median=plot(1:1:length(mogs_media_median),mogs_media_median,'ko','MarkerSize',6);
hold on;
xlabel('Bloque','FontName','Arial','FontSize',12);
ylabel(['Suma de deltas (s)'],'FontName','Arial','FontSize',12); 
title('MicroTL')
xlim([0 length(TL_media_mean)+1]);
ylim([-2.5 2.5])
yline(0)

% violines
subplot(2,4,4)
violin(Matriz_Ganancias_media,'xlabel',{'Micro MOGS','Micro MONGS','MicroTL'},'facecolor',[1 0 0; 0 0 1;0.2 0.2 0.2],'edgecolor','k','mc','k','ylabel','Sum of deltas (s)','plotlegend',0,'medc',[])
yline(0)
plot(Matriz_Ganancias_media','k.','MarkerSize',12)
ylim([-2.5 2.5])
title('Grupo con ITI')
%legend('','Mean','Median')

%% mediana x bloque
%mogs
subplot(2,4,5)
for i=1:size(mogs_group_medianaB,1)
    s=scatter(1:size(mogs_group_medianaB,2),mogs_group_medianaB(i,:),'r','filled','MarkerEdgeAlpha',0.25);
    s.MarkerFaceAlpha = 0.25;
    hold on;
end
p_mean=plot(1:1:length(mogs_mediana_mean),mogs_mediana_mean,'k.','MarkerSize',12);
%p_median=plot(1:1:length(mogs_mediana_median),mogs_mediana_median,'ko','MarkerSize',6);
hold on;
xlabel('Bloque','FontName','Arial','FontSize',12);
ylabel(['Suma de deltas (s)'],'FontName','Arial','FontSize',12); 
title('MicroMOGS')
xlim([0 length(mogs_media_mean)+1]);
ylim([-2.5 2.5])
yline(0)
%legend([p_mean p_median],'mean','median')

%mongs
subplot(2,4,6)
for i=1:size(mongs_group_medianaB,1)
    s=scatter(1:size(mongs_group_medianaB,2),mongs_group_medianaB(i,:),'b','filled','MarkerEdgeAlpha',0.25);
    s.MarkerFaceAlpha = 0.25;
    hold on;
end
p_mean=plot(1:1:length(mongs_mediana_mean),mongs_mediana_mean,'k.','MarkerSize',12);
%p_median=plot(1:1:length(mongs_mediana_median),mongs_mediana_median,'ko','MarkerSize',6);
hold on;
xlabel('Bloque','FontName','Arial','FontSize',12);
ylabel(['Suma de deltas (s)'],'FontName','Arial','FontSize',12); 
title('MicroMONGS')
xlim([0 length(mongs_media_mean)+1]);
ylim([-2.5 2.5])
yline(0)
%legend([p_mean p_median],'mean','median')

%Micro TL
subplot(2,4,7)
for i=1:size(TL_group_medianaB,1)
    s=scatter(1:size(TL_group_medianaB,2),TL_group_medianaB(i,:),'k','filled','MarkerEdgeAlpha',0.25);
    s.MarkerFaceAlpha = 0.25;
    hold on;
end
p_mean=plot(1:1:length(TL_mediana_mean),TL_mediana_mean,'k.','MarkerSize',12);
%p_median=plot(1:1:length(mogs_mediana_median),mogs_mediana_median,'ko','MarkerSize',6);
hold on;
xlabel('Bloque','FontName','Arial','FontSize',12);
ylabel(['Suma de deltas (s)'],'FontName','Arial','FontSize',12); 
title('MicroTL')
xlim([0 length(TL_mediana_mean)+1]);
ylim([-2.5 2.5])
yline(0)

% violines
subplot(2,4,8)
violin(Matriz_Ganancias_mediana,'xlabel',{'Micro MOGS','Micro MONGS','MicroTL'},'facecolor',[1 0 0; 0 0 1;0.2 0.2 0.2],'edgecolor','k','mc','k','ylabel','Sum of deltas (s)','plotlegend',0,'medc',[])
yline(0)
plot(Matriz_Ganancias_mediana','k.','MarkerSize',12)
ylim([-2.5 2.5])
title('Grupo sin ITI')
%legend('','Mean','Median')