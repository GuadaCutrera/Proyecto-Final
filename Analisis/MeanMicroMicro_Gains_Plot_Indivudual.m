%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 Guadalupe Cutrera 2/8/23                            %%%
%%% Esta funcion recibe la media/o mediana (las variables dicen mean, pero%
%%% el calculo es el mismo, asi que se le puede pasar como parámetro la %%%
%%% mediana) de los Micro Micro de cada bloque de la tarea, las grafica con
%%% su respectivo desvío y hace otro gráfico para el acumulado.         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [MicroMogs_mean_acum,MicroMongs_mean_acum,MicroTL_mean_acum]=MeanMicroMicro_Gains_Plot_Indivudual(MicroMogs_mean,MicroMongs_mean,MicroTL_mean,MicroMogs_std,MicroMongs_std,MicroTL_std,titulo)
figure; set(gcf,'Color','white'); box OFF; hold on; sgtitle([titulo ' value per block'])

% -------------------------- Micro MOGS -----------------------------------
%if seq_results(1,1).flag_norm==0
subplot(2,3,1)
plot(MicroMogs_mean,'r.')
errorbar(MicroMogs_mean,MicroMogs_std,'s')
ylim([-1 1])
% else
%     subplot(1,2,1)
%     plot(MicroMogs_mean,'r.')
%     errorbar(MicroMogs_mean,MicroMogs_std,'s')
%     ylim([-5 5])
% end
yline(0)
ylabel('MicroMogs [s]')
xlabel('Block')
title('No acum')

%-------------------------- Micro MONGS---------------------------------- 
%if seq_results(1,1).flag_norm==0
subplot(2,3,2)
plot(MicroMongs_mean,'b.')
errorbar(MicroMongs_mean,MicroMongs_std,'s')
ylim([-1 1])
% else
%     subplot(1,2,2)
%     plot(MicroMongs_mean,'b.')
%     errorbar(MicroMongs_mean,MicroMongs_std,'s')
%     ylim([-5 5])
% end
yline(0)
ylabel('MicroMongs [s]')
xlabel('Block')
title('No acum')

subplot(2,3,3)
plot(MicroTL_mean,'k.')
errorbar(MicroTL_mean,MicroTL_std,'s')
ylim([-1 1])
% else
%     subplot(1,2,2)
%     plot(MicroMongs_mean,'b.')
%     errorbar(MicroMongs_mean,MicroMongs_std,'s')
%     ylim([-5 5])
% end
yline(0)
ylabel('MicroTL [s]')
xlabel('Block')
title('No acum')

%------------------------Micro MOGS Acum-----------------------------------

MicroMogs_mean_acum=cumsum(MicroMogs_mean,'omitnan');
subplot(2,3,4)
plot(MicroMogs_mean_acum,'r.')
% errorbar(MicroMogs_mean_acum,MicroMogs_std,'s')
ylim([-3 3])
yline(0)
ylabel('MicroMogs [s]')
xlabel('Block')
title('Acum')

%------------------------Micro MONGS Acum----------------------------------

MicroMongs_mean_acum=cumsum(MicroMongs_mean,'omitnan');
subplot(2,3,5)
plot(MicroMongs_mean_acum,'b.')
% errorbar(MicroMongs_mean_acum,MicroMongs_std,'s')
ylim([-3 3])
yline(0)
ylabel('MicroMongs [s]')
xlabel('Block')
title('Acum')


%----------------------- Micro TL acum-------------------------------
MicroTL_mean_acum=cumsum(MicroTL_mean,'omitnan');
subplot(2,3,6)
plot(MicroTL_mean_acum,'k.')
% errorbar(MicroMongs_mean_acum,MicroMongs_std,'s')
ylim([-3 3])
yline(0)
ylabel('MicroTL [s]')
xlabel('Block')
title('Acum')