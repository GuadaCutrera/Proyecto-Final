%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Display de figuras para el proyecto. A partir de los resultados
%%% grupales se generan figuras de 
%%%     - Cantidad de secuencias correctas
%%%     - Intervalo entre teclas
%%%     - Ganancias a nivel de bloques (MOGS y MONGS)
%%%     - Ganancias a nivel de ensayos (MicroMOGS y MicroMONGS)
%%% Se muestran los datos hasta el bloque 11.
%%%                     Guadalupe Cutrera 2/3/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
%%-------------------------------------------------------------------------

%SIMKEYPRESS
Group_Results_wo_ITI= load('C:\Users\guada\Desktop\LFA FMED\Codigos Guada\Guada_2023\AnalisisGuada\Behaviour\Piloto MSL time sin ITI\Group_Results.mat');
Group_Results_w_ITI_1s= load('C:\Users\guada\Desktop\LFA FMED\Codigos Guada\Guada_2023\AnalisisGuada\Behaviour\Piloto MSL time 1s\Group_Results.mat');
pathGuardar='C:\Users\guada\Desktop\LFA FMED\Codigos Guada\Guada_2023\ComparePlots\';


%% --------------------- SECUENCIAS CORRECTAS -----------------------------

Group_Results_w_ITI_1s.Group_Results.seq_group_num_mean = nanmean(Group_Results_w_ITI_1s.Group_Results.Group_Parameters.seq_group_num,1);
Group_Results_w_ITI_1s.Group_Results.seq_group_num_std = nanstd(Group_Results_w_ITI_1s.Group_Results.Group_Parameters.seq_group_num,0,1)/sqrt(16);

Group_Results_wo_ITI.Group_Results.seq_group_num_mean = nanmean(Group_Results_wo_ITI.Group_Results.Group_Parameters.seq_group_num,1);
Group_Results_wo_ITI.Group_Results.seq_group_num_std = nanstd(Group_Results_wo_ITI.Group_Results.Group_Parameters.seq_group_num,0,1)/sqrt(16);

figure(1); set(gcf,'Color','white'); box ON; hold on;
%set(gcf,'Position',get(0,'ScreenSize'));
%subplot(1,2,2);  
p1=plot(Group_Results_w_ITI_1s.Group_Results.seq_group_num_mean(1:11),'.','MarkerSize', 17,'Color','#7E2F8E','LineWidth',1,'DisplayName','CON ITI');
errorbar(1:length(Group_Results_w_ITI_1s.Group_Results.seq_group_num_mean(1:11)),Group_Results_w_ITI_1s.Group_Results.seq_group_num_mean(1:11),Group_Results_w_ITI_1s.Group_Results.seq_group_num_std(1:11),'Color','#7E2F8E','LineWidth',1)
xlabel('Bloque','FontName','Arial','FontSize',12);
ylabel('Secuencias Correctas','FontName','Arial','FontSize',12);
%title('Grupo CON ITI')
xlim([0.7 11.3]);



%subplot(1,2,1); 
p2=plot(Group_Results_wo_ITI.Group_Results.seq_group_num_mean(1:11),'.','MarkerSize', 17,'Color','#D95319','LineWidth',1,'DisplayName','SIN ITI');
errorbar(1:length(Group_Results_wo_ITI.Group_Results.seq_group_num_mean(1:11)),Group_Results_wo_ITI.Group_Results.seq_group_num_mean(1:11),Group_Results_wo_ITI.Group_Results.seq_group_num_std(1:11),'Color','#D95319','LineWidth',1)
xlabel('Bloque','FontName','Arial','FontSize',12);
ylabel('Secuencias Correctas','FontName','Arial','FontSize',12);
%title('Grupo SIN ITI')
ylim([1.5 5.5]);
xlim([0.7 11.3]);
legend([p1 p2],{'Grupo CON ITI','Grupo SIN ITI'})

%% -------------------------- IKI corr -----------------------------------

Matriz_iki(1,:)=Group_Results_wo_ITI.Group_Results.iki_visual_group_corr_mean;
Matriz_iki(2,:)=Group_Results_w_ITI_1s.Group_Results.iki_visual_group_corr_mean;
std_iki(1,:)=Group_Results_wo_ITI.Group_Results.iki_visual_group_corr_std;
std_iki(2,:)=Group_Results_w_ITI_1s.Group_Results.iki_visual_group_corr_std;

LearningCompare(Matriz_iki,std_iki,2)
ax = gca;
ax.YDir = 'reverse';
ylim([-0.3 0.5])

%%
%promedio por bloque
iki_promedio_wo_ITI=[];
iki_promedio_w_ITI=[];
cont=1;
for i=1:10:size(Group_Results_wo_ITI.Group_Results.Group_Parameters.iki_visual_group_corr,2)
    iki_promedio_wo_ITI(:,cont)=nanmean(Group_Results_wo_ITI.Group_Results.Group_Parameters.iki_visual_group_corr(:,i:i+9),2)
    iki_promedio_w_ITI(:,cont)=nanmean(Group_Results_w_ITI_1s.Group_Results.Group_Parameters.iki_visual_group_corr(:,i:i+9),2)
    cont=cont+1;
end

% Diferencia entre los cinco primeros y cinco ultimos bloques
Dif_Block_wo_ITI=nanmean(iki_promedio_wo_ITI(:,1:2),2)-nanmean(iki_promedio_wo_ITI(:,10:11),2);
Dif_Block_w_ITI=nanmean(iki_promedio_w_ITI(:,1:2),2)-nanmean(iki_promedio_w_ITI(:,10:11),2);

%% ------------------------- VIOLINES MICRO GAINS CORR -------------------

%% BLOQUE 11
Plot_violin_MicroGains(Group_Results_w_ITI_1s.Group_Results.TL_group_acumulado_corr_mean(1:11),...
    Group_Results_w_ITI_1s.Group_Results.mogs_group_acumulado_corr_mean(1:11),Group_Results_w_ITI_1s.Group_Results.mongs_group_acumulado_corr_mean(1:11),...
    Group_Results_w_ITI_1s.Group_Results.Group_Parameters.TL_acumulado_group_corr(:,1:11),Group_Results_w_ITI_1s.Group_Results.Group_Parameters.mogs_acumulado_group_corr(:,1:11),...
    Group_Results_w_ITI_1s.Group_Results.Group_Parameters.mongs_acumulado_group_corr(:,1:11),'With ITI')


Plot_violin_MicroGains(Group_Results_wo_ITI.Group_Results.TL_group_acumulado_corr_mean(1:11),...
    Group_Results_wo_ITI.Group_Results.mogs_group_acumulado_corr_mean(1:11),Group_Results_wo_ITI.Group_Results.mongs_group_acumulado_corr_mean(1:11),...
    Group_Results_wo_ITI.Group_Results.Group_Parameters.TL_acumulado_group_corr(:,1:11),Group_Results_wo_ITI.Group_Results.Group_Parameters.mogs_acumulado_group_corr(:,1:11),...
    Group_Results_wo_ITI.Group_Results.Group_Parameters.mongs_acumulado_group_corr(:,1:11),'Without ITI')


%% -------------------- Violines Micro Micro Gains ------------------------

% BLOQUE 11
Plot_violin_MicroMicroGains(Group_Results_w_ITI_1s.Group_Results.media_MicroMogs_corr_acum_mean(1:11),...
    Group_Results_w_ITI_1s.Group_Results.media_MicroMongs_corr_acum_mean(1:11),...
    Group_Results_w_ITI_1s.Group_Results.media_MicroMogs_corr_acum_median(1:11),...
    Group_Results_w_ITI_1s.Group_Results.media_MicroMongs_corr_acum_median(1:11),...
    Group_Results_w_ITI_1s.Group_Results.mediana_MicroMogs_corr_acum_mean(1:11),...
    Group_Results_w_ITI_1s.Group_Results.mediana_MicroMongs_corr_acum_mean(1:11),...
    Group_Results_w_ITI_1s.Group_Results.mediana_MicroMogs_corr_acum_median(1:11),...
    Group_Results_w_ITI_1s.Group_Results.mediana_MicroMongs_corr_acum_median(1:11),...
    Group_Results_w_ITI_1s.Group_Results.Group_Parameters.media_MicroMOGS_corr_acum_group(:,1:11),...
    Group_Results_w_ITI_1s.Group_Results.Group_Parameters.media_MicroMONGS_corr_acum_group(:,1:11),...
    Group_Results_w_ITI_1s.Group_Results.Group_Parameters.mediana_MicroMOGS_corr_acum_group(:,1:11),...
    Group_Results_w_ITI_1s.Group_Results.Group_Parameters.mediana_MicroMONGS_corr_acum_group(:,1:11),'With ITI')

Plot_violin_MicroMicroGains(Group_Results_wo_ITI.Group_Results.media_MicroMogs_corr_acum_mean(1:11),...
    Group_Results_wo_ITI.Group_Results.media_MicroMongs_corr_acum_mean(1:11),...
    Group_Results_wo_ITI.Group_Results.media_MicroMogs_corr_acum_median(1:11),...
    Group_Results_wo_ITI.Group_Results.media_MicroMongs_corr_acum_median(1:11),...
    Group_Results_wo_ITI.Group_Results.mediana_MicroMogs_corr_acum_mean(1:11),...
    Group_Results_wo_ITI.Group_Results.mediana_MicroMongs_corr_acum_mean(1:11),...
    Group_Results_wo_ITI.Group_Results.mediana_MicroMogs_corr_acum_median(1:11),...
    Group_Results_wo_ITI.Group_Results.mediana_MicroMongs_corr_acum_median(1:11),...
    Group_Results_wo_ITI.Group_Results.Group_Parameters.media_MicroMOGS_corr_acum_group(:,1:11),...
    Group_Results_wo_ITI.Group_Results.Group_Parameters.media_MicroMONGS_corr_acum_group(:,1:11),...
    Group_Results_wo_ITI.Group_Results.Group_Parameters.mediana_MicroMOGS_corr_acum_group(:,1:11),...
    Group_Results_wo_ITI.Group_Results.Group_Parameters.mediana_MicroMONGS_corr_acum_group(:,1:11),'Without ITI')


%%
Plot_violin_MicroMicroGains(Group_Results_w_ITI_1s.Group_Results.media_MicroMogs_corr_acum_mean(1:11),...
    Group_Results_w_ITI_1s.Group_Results.media_MicroMongs_corr_acum_mean(1:11),...
    Group_Results_w_ITI_1s.Group_Results.media_MicroTL_corr_acum_mean(1:11),...
    Group_Results_w_ITI_1s.Group_Results.media_MicroMogs_corr_acum_median(1:11),...
    Group_Results_w_ITI_1s.Group_Results.media_MicroMongs_corr_acum_median(1:11),...
    Group_Results_wo_ITI.Group_Results.media_MicroMogs_corr_acum_mean(1:11),...
    Group_Results_wo_ITI.Group_Results.media_MicroMongs_corr_acum_mean(1:11),...
    Group_Results_wo_ITI.Group_Results.media_MicroTL_corr_acum_mean(1:11),...
    Group_Results_wo_ITI.Group_Results.media_MicroMogs_corr_acum_median(1:11),...
    Group_Results_wo_ITI.Group_Results.media_MicroMongs_corr_acum_median(1:11),...
    Group_Results_w_ITI_1s.Group_Results.Group_Parameters.media_MicroMOGS_corr_acum_group(:,1:11),...
    Group_Results_w_ITI_1s.Group_Results.Group_Parameters.media_MicroMONGS_corr_acum_group(:,1:11),...
    Group_Results_w_ITI_1s.Group_Results.Group_Parameters.media_MicroTL_corr_acum_group(:,1:11),...
    Group_Results_wo_ITI.Group_Results.Group_Parameters.media_MicroMOGS_corr_acum_group(:,1:11),...
    Group_Results_wo_ITI.Group_Results.Group_Parameters.media_MicroMONGS_corr_acum_group(:,1:11),...
    Group_Results_wo_ITI.Group_Results.Group_Parameters.media_MicroTL_corr_acum_group(:,1:11),'Ganancias Micro Micro')


%% MICRO MICRO

for i=1:16
    MicroMOGS_corr(i,1)=Group_Results_wo_ITI.Group_Results.Group_Parameters.media_MicroMOGS_corr_acum_group(i,11);
    MicroMONGS_corr(i,1)=Group_Results_wo_ITI.Group_Results.Group_Parameters.media_MicroMONGS_corr_acum_group(i,11);
    MicroTL_corr(i,1)=Group_Results_wo_ITI.Group_Results.Group_Parameters.media_MicroTL_corr_acum_group(i,11);
    MicroMOGS_corr(i,2)=Group_Results_w_ITI_1s.Group_Results.Group_Parameters.media_MicroMOGS_corr_acum_group(i,11);
    MicroMONGS_corr(i,2)=Group_Results_w_ITI_1s.Group_Results.Group_Parameters.media_MicroMONGS_corr_acum_group(i,11);
    MicroTL_corr(i,2)=Group_Results_w_ITI_1s.Group_Results.Group_Parameters.media_MicroTL_corr_acum_group(i,11);

end


% WITHOUT ITI
Matriz_Gains_acum(:,1)=MicroTL_corr(:,1);
Matriz_Gains_acum(:,2)=MicroMOGS_corr(:,1);
Matriz_Gains_acum(:,3)=MicroMONGS_corr(:,1);

Matriz_Gains_mean(1,:)=Group_Results_wo_ITI.Group_Results.media_MicroTL_corr_acum_mean(1:11);
Matriz_Gains_mean(2,:)=Group_Results_wo_ITI.Group_Results.media_MicroMogs_corr_acum_mean(1:11);
Matriz_Gains_mean(3,:)=Group_Results_wo_ITI.Group_Results.media_MicroMongs_corr_acum_mean(1:11);

Matriz_Gains_std(1,:)=Group_Results_wo_ITI.Group_Results.media_MicroTL_corr_acum_std(1:11);
Matriz_Gains_std(2,:)=Group_Results_wo_ITI.Group_Results.media_MicroMogs_corr_acum_std(1:11);
Matriz_Gains_std(3,:)=Group_Results_wo_ITI.Group_Results.media_MicroMongs_corr_acum_std(1:11);

Plot_Gains_and_violin(Matriz_Gains_mean,Matriz_Gains_std,Matriz_Gains_acum')
sgtitle('SIN ITI')

%% WITH ITI
Matriz_Gains_acum(:,1)=MicroTL_corr(:,2);
Matriz_Gains_acum(:,2)=MicroMOGS_corr(:,2);
Matriz_Gains_acum(:,3)=MicroMONGS_corr(:,2);

Matriz_Gains_mean(1,:)=Group_Results_w_ITI_1s.Group_Results.media_MicroTL_corr_acum_mean(1:11);
Matriz_Gains_mean(2,:)=Group_Results_w_ITI_1s.Group_Results.media_MicroMogs_corr_acum_mean(1:11);
Matriz_Gains_mean(3,:)=Group_Results_w_ITI_1s.Group_Results.media_MicroMongs_corr_acum_mean(1:11);

Matriz_Gains_std(1,:)=Group_Results_w_ITI_1s.Group_Results.media_MicroTL_corr_acum_std(1:11);
Matriz_Gains_std(2,:)=Group_Results_w_ITI_1s.Group_Results.media_MicroMogs_corr_acum_std(1:11);
Matriz_Gains_std(3,:)=Group_Results_w_ITI_1s.Group_Results.media_MicroMongs_corr_acum_std(1:11);

%w ITI
Plot_Gains_and_violin(Matriz_Gains_mean,Matriz_Gains_std,Matriz_Gains_acum')
sgtitle('CON ITI')

%% MICRO GAINS


for i=1:16
    MicroMOGS_corr(i,1)=Group_Results_wo_ITI.Group_Results.Group_Parameters.mogs_acumulado_group_corr(i,11);
    MicroMONGS_corr(i,1)=Group_Results_wo_ITI.Group_Results.Group_Parameters.mongs_acumulado_group_corr(i,11);
    MicroTL_corr(i,1)=Group_Results_wo_ITI.Group_Results.Group_Parameters.TL_acumulado_group_corr(i,11);
    MicroMOGS_corr(i,2)=Group_Results_w_ITI_1s.Group_Results.Group_Parameters.mogs_acumulado_group_corr(i,11);
    MicroMONGS_corr(i,2)=Group_Results_w_ITI_1s.Group_Results.Group_Parameters.mongs_acumulado_group_corr(i,11);
    MicroTL_corr(i,2)=Group_Results_w_ITI_1s.Group_Results.Group_Parameters.TL_acumulado_group_corr(i,11);

end


% WITHOUT ITI
Matriz_Gains_acum(:,1)=MicroTL_corr(:,1);
Matriz_Gains_acum(:,2)=MicroMOGS_corr(:,1);
Matriz_Gains_acum(:,3)=MicroMONGS_corr(:,1);

Matriz_Gains_mean(1,:)=Group_Results_wo_ITI.Group_Results.TL_group_acumulado_corr_mean(1:11);
Matriz_Gains_mean(2,:)=Group_Results_wo_ITI.Group_Results.mogs_group_acumulado_corr_mean(1:11);
Matriz_Gains_mean(3,:)=Group_Results_wo_ITI.Group_Results.mongs_group_acumulado_corr_mean(1:11);

Matriz_Gains_std(1,:)=Group_Results_wo_ITI.Group_Results.TL_group_acumulado_corr_std(1:11);
Matriz_Gains_std(2,:)=Group_Results_wo_ITI.Group_Results.mogs_group_acumulado_corr_std(1:11);
Matriz_Gains_std(3,:)=Group_Results_wo_ITI.Group_Results.mongs_group_acumulado_corr_std(1:11);

Plot_Gains_and_violin(Matriz_Gains_mean,Matriz_Gains_std,Matriz_Gains_acum')
sgtitle('SIN ITI')

%% WITH ITI
Matriz_Gains_acum(:,1)=MicroTL_corr(:,2);
Matriz_Gains_acum(:,2)=MicroMOGS_corr(:,2);
Matriz_Gains_acum(:,3)=MicroMONGS_corr(:,2);

Matriz_Gains_mean(1,:)=Group_Results_w_ITI_1s.Group_Results.TL_group_acumulado_corr_mean(1:11);
Matriz_Gains_mean(2,:)=Group_Results_w_ITI_1s.Group_Results.mogs_group_acumulado_corr_mean(1:11);
Matriz_Gains_mean(3,:)=Group_Results_w_ITI_1s.Group_Results.mongs_group_acumulado_corr_mean(1:11);

Matriz_Gains_std(1,:)=Group_Results_w_ITI_1s.Group_Results.TL_group_acumulado_corr_std(1:11);
Matriz_Gains_std(2,:)=Group_Results_w_ITI_1s.Group_Results.mogs_group_acumulado_corr_std(1:11);
Matriz_Gains_std(3,:)=Group_Results_w_ITI_1s.Group_Results.mongs_group_acumulado_corr_std(1:11);

Plot_Gains_and_violin(Matriz_Gains_mean,Matriz_Gains_std,Matriz_Gains_acum')
sgtitle('CON ITI')

%% CORRELACION MICRO MICRO OFF Y MICRO OFF

%ultimo del acumulado
for i=1:16
    MOGS_corr(i,1)=Group_Results_wo_ITI.Group_Results.Group_Parameters.mogs_acumulado_group_corr(i,11);
    MOGS_corr(i,2)=Group_Results_w_ITI_1s.Group_Results.Group_Parameters.mogs_acumulado_group_corr(i,11);
    MicroMOGS_corr(i,1)=Group_Results_wo_ITI.Group_Results.Group_Parameters.media_MicroMOGS_corr_acum_group(i,11);
    MicroMOGS_corr(i,2)=Group_Results_w_ITI_1s.Group_Results.Group_Parameters.media_MicroMOGS_corr_acum_group(i,11);
end

figure; set(gcf,'Color','white'); box OFF; hold on; %sgtitle(titulo)
set(gcf,'Position',get(0,'ScreenSize'));
%SIN ITI
    subplot(1,2,1)
    plot(MOGS_corr(:,1),MicroMOGS_corr(:,1),'Marker','.','Color','#D95319','LineStyle','none','MarkerSize', 17)
    hold on;
    b = regress(MicroMOGS_corr(:,1), [ones(size(MOGS_corr,1),1) MOGS_corr(:,1)]); % Obtener los coeficientes de la recta de regresión
    [r,p] = corrcoef(MOGS_corr(:,1)', MicroMOGS_corr(:,1)'); % Calcular el coeficiente de correlación R2
    r2 = r(1,2)^2;
    plot(MOGS_corr(:,1), b(1) + b(2)*MOGS_corr(:,1),'k','LineWidth',1); % Graficar la recta de regresión
    text(0.2,0.8,{['R^2 =' num2str(r2) ', p-valor = ' num2str(p(2))],['f(x) = ' num2str(b(2)) 'x ' num2str(b(1))]});
    xlabel('MOGS','FontName','Arial','FontSize',12);
    ylabel('Micro MOGS','FontName','Arial','FontSize',12);
    title('Grupo SIN ITI')
    clear b; clear r; clear r2; clear p;
    ylim([-2 1])
    xlim([-0.5 1])
    yline(0)
    xline(0)
    
%CON ITI
    subplot(1,2,2)
    plot(MOGS_corr(:,2),MicroMOGS_corr(:,2),'Marker','.','Color','#7E2F8E','LineStyle','none','MarkerSize', 17)
    hold on;
    b = regress(MicroMOGS_corr(:,2), [ones(size(MOGS_corr,1),1) MOGS_corr(:,2)]); % Obtener los coeficientes de la recta de regresión
    [r,p] = corrcoef(MOGS_corr(:,2)', MicroMOGS_corr(:,2)'); % Calcular el coeficiente de correlación R2
    r2 = r(1,2)^2;
    plot(MOGS_corr(:,2), b(1) + b(2)*MOGS_corr(:,2),'k','LineWidth',1); % Graficar la recta de regresión
    text(0.2,0.8,{['R^2 = ' num2str(r2) ', p-valor = ' num2str(p(2))],['f(x) = ' num2str(b(2)) 'x +' num2str(b(1))]});
    xlabel('MOGS','FontName','Arial','FontSize',12);
    ylabel('Micro MOGS','FontName','Arial','FontSize',12);
    title('Grupo CON ITI')
    ylim([-2 1])
    xlim([-0.5 1])
    yline(0)
    xline(0)
    clear b; clear r; clear r2; clear p;
    
%%
figure; set(gcf,'Color','white'); box OFF; hold on; %sgtitle(titulo)
set(gcf,'Position',get(0,'ScreenSize'));
%SIN ITI
    subplot(1,2,1)
    plot(MicroMOGS_corr(:,1),MOGS_corr(:,1),'Marker','.','Color','#D95319','LineStyle','none','MarkerSize', 17)
    hold on;
    b = regress(MOGS_corr(:,1), [ones(size(MicroMOGS_corr,1),1) MicroMOGS_corr(:,1)]); % Obtener los coeficientes de la recta de regresión
    [r,p] = corrcoef(MicroMOGS_corr(:,1)', MOGS_corr(:,1)'); % Calcular el coeficiente de correlación R2
    r2 = r(1,2)^2;
    plot(MicroMOGS_corr(:,1), b(1) + b(2)*MicroMOGS_corr(:,1),'k','LineWidth',1); % Graficar la recta de regresión
    %text(0.2,0.8,{['R^2 =' num2str(r2) ', p-valor = ' num2str(p(2))],['f(x) = ' num2str(b(2)) 'x ' num2str(b(1))]});
    xlabel('MicroMOGS','FontName','Arial','FontSize',12);
    ylabel('MOGS','FontName','Arial','FontSize',12);
    title('Grupo SIN ITI')
    clear b; clear r; clear r2; clear p;
    ylim([-0.6 1])
    xlim([-2 0.7])
    yline(0)
    xline(0)
    
%CON ITI
    subplot(1,2,2)
    plot(MicroMOGS_corr(:,2),MOGS_corr(:,2),'Marker','.','Color','#7E2F8E','LineStyle','none','MarkerSize', 17)
    hold on;
    b = regress(MOGS_corr(:,2), [ones(size(MicroMOGS_corr,1),1) MicroMOGS_corr(:,2)]); % Obtener los coeficientes de la recta de regresión
    [r,p] = corrcoef(MicroMOGS_corr(:,2)', MOGS_corr(:,2)'); % Calcular el coeficiente de correlación R2
    r2 = r(1,2)^2;
    plot(MicroMOGS_corr(:,2), b(1) + b(2)*MicroMOGS_corr(:,2),'k','LineWidth',1); % Graficar la recta de regresión
    %text(0.2,0.8,{['R^2 = ' num2str(r2) ', p-valor = ' num2str(p(2))],['f(x) = ' num2str(b(2)) 'x +' num2str(b(1))]});
    xlabel('MicroMOGS','FontName','Arial','FontSize',12);
    ylabel(' MOGS','FontName','Arial','FontSize',12);
    title('Grupo CON ITI')
    ylim([-0.6 1])
    xlim([-2 0.7])
    yline(0)
    xline(0)
    clear b; clear r; clear r2; clear p;
    
