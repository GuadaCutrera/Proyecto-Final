%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               Plots Micro Gains                                     %%%
%%%                                                                     %%%
%%% Acum MOGS: acumulative sum of Micro Offline Gains per block         %%%
%%% Acum MONGS: acumulative sum of Micro Online Gains per block         %%%
%%% Acum TL: acumulative sum of Total Learning per block                %%%
%%% Not Acum MOGS: Micro Offline Gains per block                        %%%
%%% Not Acum MONGS: Micro Online Gains per block                        %%%
%%% Not Acum TL: Total Learning per block                               %%%
%%%                      Guadalupe cutrera 5/8/23                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MicroGains_Plot_Individual(Mogs,Mongs,TL,Mogs_acum,Mongs_acum,TL_acum,titulo,seq_results)

figure; set(gcf,'Color','white'); box OFF; hold on; sgtitle(titulo)

%% MOGS acumulado
subplot(2,3,3)
plot(Mogs_acum,'r.');
xlabel('Block','FontName','Arial','FontSize',12)
ylabel('MOGs','FontName','Arial','FontSize',12);
if seq_results(1,1).flag_norm==0 || contains(titulo,'crudo')
    ylim([-1.5 1.5])
else
    if strcmp(seq_results(1,1).flag_tipo_norm,'Z')
        ylim([-12 12])
    else
        ylim([-5 5])
    end
end
xlim([0 length(Mogs_acum)]);
yline(0)


%% MONGS acumulado
subplot(2,3,2)
plot(Mongs_acum,'b.');
xlabel('Block','FontName','Arial','FontSize',12)
ylabel('MONGs','FontName','Arial','FontSize',12);
if seq_results(1,1).flag_norm==0 || contains(titulo,'crudo')
    ylim([-1.5 1.5])
else
    if strcmp(seq_results(1,1).flag_tipo_norm,'Z')
        ylim([-12 12])
    else
        ylim([-5 5])
    end
end
xlim([0 length(Mongs_acum)]);
yline(0)

%% TOTAL LEARNING acumulado
subplot(2,3,1)
plot(TL_acum,'k.');
xlabel('Block','FontName','Arial','FontSize',12)
ylabel('Total Learning','FontName','Arial','FontSize',12);
if seq_results(1,1).flag_norm==0 || contains(titulo,'crudo')
    ylim([-1.5 1.5])
else
    if strcmp(seq_results(1,1).flag_tipo_norm,'Z')
        ylim([-12 12])
    else
        ylim([-5 5])
    end
end
xlim([0 length(TL_acum)]);
title('Acumulado')
yline(0)

%% MOGS no acumulado
subplot(2,3,6)
plot(Mogs,'r.');
xlabel('Block','FontName','Arial','FontSize',12)
ylabel('MOGs','FontName','Arial','FontSize',12);
if seq_results(1,1).flag_norm==0 || contains(titulo,'crudo')
    ylim([-0.5 0.5])
else
    ylim([-5 5])
end
xlim([0 length(Mogs)]);
yline(0)

%% MONGS no acumulado
subplot(2,3,5)
plot(Mongs,'b.');
xlabel('Block','FontName','Arial','FontSize',12)
ylabel('MONGs','FontName','Arial','FontSize',12);
if seq_results(1,1).flag_norm==0 || contains(titulo,'crudo')
    ylim([-0.5 0.5])
else
    ylim([-5 5])
end
xlim([0 length(Mongs)]);
yline(0)

%% TOTAL LEARNING no acumulado
subplot(2,3,4)
plot(TL,'k.');
xlabel('Block','FontName','Arial','FontSize',12)
ylabel('Total Learning','FontName','Arial','FontSize',12);
if seq_results(1,1).flag_norm==0 || contains(titulo,'crudo')
    ylim([-0.5 0.5])
else
    ylim([-5 5])
end
xlim([0 length(TL)]);
title('No acumulado')
yline(0)