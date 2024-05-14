function MicroMicroGains_Plot_Grupal(MicroMOGS_mean,MicroMOGS_std,MicroMONGS_mean,MicroMONGS_std,MicroMOGS_subjects,MicroMONGS_subjects,titulo,Group_Parameters)
%Guadalupe Cutrera 12/3/23
figure; set(gcf,'Color','white'); box OFF; hold on; sgtitle(titulo);
set(gcf,'Position',get(0,'ScreenSize'));

%% MICRO MOGS
subplot(3,1,1)

for i=1:size(MicroMOGS_subjects,1)
    s=scatter(1:size(MicroMOGS_subjects,2),MicroMOGS_subjects(i,:),'r','filled','MarkerEdgeAlpha',0.1);
    s.MarkerFaceAlpha = 0.4;
    hold on;
end

plot(MicroMOGS_mean,'k.','MarkerSize',10);
for j=1:(Group_Parameters.cant_SeqBlock-1):length(MicroMOGS_mean)
     xline(j); hold on;
end 
xlabel('Trial','FontName','Arial','FontSize',12)
ylabel('Micro MOGs','FontName','Arial','FontSize',12);
xlim([0 length(MicroMOGS_mean)]);
if contains(titulo,'crudo')
    if contains(titulo,'no acum')
        ylim([-2 2])
    else
        ylim([-2.5 2.5])
    end
else 
    if contains(titulo,'no acum')
         if Group_Parameters.flag_norm==0
                ylim([-1 1])
            else 
                ylim([-5 5])
         end
    else
        if Group_Parameters.flag_norm==0
                ylim([-2 2])
            else 
                ylim([-10 10])
        end
    end
end
yline(0);

%% MICRO MONGS
subplot(3,1,2)

for i=1:size(MicroMONGS_subjects,1)
    s=scatter(1:size(MicroMONGS_subjects,2),MicroMONGS_subjects(i,:),'b','filled','MarkerEdgeAlpha',0.1);
    s.MarkerFaceAlpha = 0.4;
    hold on;
end

plot(MicroMONGS_mean,'k.','MarkerSize',10);
for j=1:Group_Parameters.cant_SeqBlock:length(MicroMONGS_mean)
     xline(j); hold on;
end 
xlabel('Trial','FontName','Arial','FontSize',12)
ylabel('Micro MONGs','FontName','Arial','FontSize',12);
xlim([0 length(MicroMONGS_mean)]);
if contains(titulo,'crudo')
    if contains(titulo,'no acum')
        ylim([-2 2])
    else
        ylim([-2.5 2.5])
    end
else 
    if contains(titulo,'no acum')
         if Group_Parameters.flag_norm==0
                ylim([-1 1])
            else 
                ylim([-5 5])
         end
    else
        if Group_Parameters.flag_norm==0
                ylim([-2 2])
            else 
                ylim([-10 10])
        end
    end
end
yline(0);

% GRAFICO COMPLETO 
subplot(3,1,3)

% para que sean comparables mogs y mongs debo agregar un NaN al final de 
% cada bloque de mogs asi ambos tienen la misma cantidad de elementos en cada bloque
aux_micro_mogs_mean=[];
for i=1:(Group_Parameters.cant_SeqBlock-1):(length(MicroMOGS_mean))
    aux_micro_mogs_mean=[aux_micro_mogs_mean MicroMOGS_mean(i:(i-1)+(Group_Parameters.cant_SeqBlock-1)) NaN];
end

%mogs

    nan_cont=0;
    for i=1:(Group_Parameters.cant_SeqBlock-1):length(MicroMOGS_mean)
        x_vector = [(i+nan_cont):((i+nan_cont-1)+(Group_Parameters.cant_SeqBlock-1)), fliplr((i+nan_cont):((i+nan_cont-1)+(Group_Parameters.cant_SeqBlock-1)))];
        patch = fill(x_vector, [MicroMOGS_mean(i:(i-1)+(Group_Parameters.cant_SeqBlock-1)) + MicroMOGS_std(i:(i-1)+(Group_Parameters.cant_SeqBlock-1)) , fliplr(MicroMOGS_mean(i:(i-1)+(Group_Parameters.cant_SeqBlock-1)) - MicroMOGS_std(i:(i-1)+(Group_Parameters.cant_SeqBlock-1)))], [243 169 114]./255);
        set(patch, 'edgecolor', 'none');
        set(patch, 'FaceAlpha', 0.3);
        hold on;
        nan_cont=nan_cont+1;
    end
    %media
    plot(aux_micro_mogs_mean,'r','LineWidth',0.5);
%mongs

  %error estandar
    x_vector = [1:length(MicroMONGS_mean), fliplr(1:length(MicroMONGS_mean))];
    patch = fill(x_vector, [MicroMONGS_mean + MicroMONGS_std, fliplr(MicroMONGS_mean - MicroMONGS_std)],[128 193 219]./255);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', 0.3);
    hold on;
    %media
    plot(MicroMONGS_mean,'b','LineWidth',0.5);
   
    if contains(titulo,'crudo')
        if contains(titulo,'no acum')
            ylim([-0.5 0.5])
        else
            ylim([-1 1])
        end
    else
        if contains(titulo,'no acum')
            if Group_Parameters.flag_norm==0
                ylim([-0.3 0.3])
            else %aca solo hay normalización Zscore
                ylim([-1.5 1.5])
            end
        else
            if Group_Parameters.flag_norm==0
                ylim([-0.5 0.5])
            else %aca solo hay normalización Zscore
                ylim([-3 3])
            end
        end
    end
    yline(0)