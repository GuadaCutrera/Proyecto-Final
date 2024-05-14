function MicroGains_Plot_Grupal(MOGS_mean,MOGS_std,MONGS_mean,MONGS_std,TL_mean,TL_std,MOGS_subjects,MONGS_subjects,TL_subjects,titulo,Group_Parameters)
%Guadalupe cutrera 5/8/23                       
figure; set(gcf,'Color','white'); box OFF; hold on; sgtitle(titulo)
set(gcf,'Position',get(0,'ScreenSize'));


%% MOGS 
subplot(1,4,3)
plot(MOGS_mean,'r.','MarkerSize', 15); hold on;

for i=1:size(MOGS_subjects,1)
    s=scatter(1:size(MOGS_subjects,2),MOGS_subjects(i,:),'r','filled','MarkerEdgeAlpha',0.1);
    s.MarkerFaceAlpha = 0.1;
    hold on;
end
xlabel('Block','FontName','Arial','FontSize',12)
ylabel('MOGs','FontName','Arial','FontSize',12);
if Group_Parameters.flag_norm==0 || contains(titulo,'crudo') 
    ylim([-1 1])
else
    if strcmp(Group_Parameters.flag_tipo_norm,'Z')
        ylim([-10 10])
    else %norm 01
        %ylim([-3 3])
        ylim([-1 1])
    end
end
xlim([1 length(MOGS_mean)]);
yline(0);

%% MONGS 

subplot(1,4,2)
plot(MONGS_mean,'b.','MarkerSize',15); hold on;

for i=1:size(MONGS_subjects,1)
    s=scatter(1:size(MONGS_subjects,2),MONGS_subjects(i,:),'b','filled','MarkerEdgeAlpha',0.1);
    s.MarkerFaceAlpha = 0.1;
    hold on;
end

xlabel('Block','FontName','Arial','FontSize',12)
ylabel('MONGs','FontName','Arial','FontSize',12);
if Group_Parameters.flag_norm==0 || contains(titulo,'crudo') 
    ylim([-1 1])
else
    if strcmp(Group_Parameters.flag_tipo_norm,'Z')
        ylim([-10 10])
    else %norm 01
        %ylim([-3 3])
        ylim([-1 1])
    end
end
xlim([1 length(MONGS_mean)-1]);
yline(0);

%% TOTAL LEARNING

subplot(1,4,1)
plot(TL_mean,'k.','MarkerSize',15); hold on;

for i=1:size(TL_subjects,1)
    s=scatter(1:size(TL_subjects,2),TL_subjects(i,:),'k','filled','MarkerEdgeAlpha',0.1);
    s.MarkerFaceAlpha = 0.1;
    hold on;
end

xlabel('Block','FontName','Arial','FontSize',12)
ylabel('Total Learning','FontName','Arial','FontSize',12);
if Group_Parameters.flag_norm==0 || contains(titulo,'crudo') 
    ylim([-1 1])
else
    if strcmp(Group_Parameters.flag_tipo_norm,'Z')
        ylim([-10 10])
    else %norm 01
        %ylim([-3 3])
        ylim([-1 1])
    end
end
xlim([1 length(TL_mean)]);
yline(0);

%% Grafico completo
subplot(1,4,4)

%mogs
    %error estandar
    x_vector = [1:length(MOGS_mean), fliplr(1:length(MOGS_mean))];
    patch = fill(x_vector, [MOGS_mean + MOGS_std , fliplr(MOGS_mean - MOGS_std)], [243 169 114]./255);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', 0.3);
    hold on;
    %media
    plot(MOGS_mean,'r','LineWidth',0.5);
%mongs

  %error estandar
    x_vector = [1:length(MONGS_mean), fliplr(1:length(MONGS_mean))];
    patch = fill(x_vector, [MONGS_mean + MONGS_std , fliplr(MONGS_mean - MONGS_std)],[128 193 219]./255);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', 0.3);
    hold on;
    %media
    plot(MONGS_mean,'b','LineWidth',0.5);

%Total Learning

    %error estandar
    x_vector = [1:length(TL_mean), fliplr(1:length(TL_mean))];
    patch = fill(x_vector, [TL_mean + TL_std , fliplr(TL_mean - TL_std)], [200 200 200]./255);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', 0.5);
    hold on;
    
    %media
    plot(TL_mean,'k','LineWidth',0.5);

    
if Group_Parameters.flag_norm==0 || contains(titulo,'crudo')
    ylim([-1 1])
else
    if strcmp(Group_Parameters.flag_tipo_norm,'Z')
        ylim([-10 10])
    else %norm 01
        %ylim([-3 3])
        ylim([-1 1])
    end
end 
yline(0);
xlabel('Blocks')
legend('','MOGS','', 'MONGS', '', 'Total Learning')
xlim([1 length(TL_mean)])
