%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IKI_per_trial: Interkey Interval on each correct sequence. This type%%%
%%% of plot is also known to the lab as "the learning curve"            %%%
%%% It only plots the interval of correct sequences one after another.  %%%
%%%                                                                     %%%
%%% It's also included the filtered version of IKI_per_trial            %%%
%%%                     Guadalupe Cutrera 5/8/23                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function InterKeyInterval_Plot_Individual(seq_results,titulo)

figure; set(gcf,'Color','white'); box OFF; hold on;

%% Datos crudos
if seq_results(1,1).flag_norm==1 || seq_results(1,1).flag_filt ==1 %si se corrigieron datos -> en la figura hay dos img
    subplot(2,1,1)
end
IKI_trial=reshape(seq_results(1,1).IKI_per_trial',1, []);
%noTrials=length(IKI_trial); %Number of trials: 12 seq x 15 block = 180 trials
plot(IKI_trial);
hold on

for i=1:max(seq_results.cant_SeqBlock):length(IKI_trial)
    xline(i); %indicates a new block 
end

xlabel('Trials','FontName','Arial','FontSize',12)
ylabel('Interkeys interval','FontName','Arial','FontSize',12);
ylim([0 (max([seq_results.Interval12, seq_results.Interval23, ...
              seq_results.Interval34, seq_results.Interval45]))])
xlim([0 length(IKI_trial)]);
title('Curva de Learning - Datos crudos')

%% datos corregidos
if seq_results(1,1).flag_norm==1 || seq_results(1,1).flag_filt ==1
    subplot(2,1,2)
    IKI_trial_corr=reshape(seq_results(1,1).IKI_per_trial_corr',1, []);
    plot(IKI_trial_corr);
    hold on

    for i=1:max(seq_results.cant_SeqBlock):length(IKI_trial_corr)
        xline(i);
    end

    xlabel('Trials','FontName','Arial','FontSize',12)
    ylabel('Interkeys interval','FontName','Arial','FontSize',12);
    if seq_results(1,1).flag_norm==0
        %le dejo el mismo limite superior que en el no filtrado para conservar la
        %escala del gràfico
        ylim([0 (max([seq_results.Interval12, seq_results.Interval23, ...
                      seq_results.Interval34, seq_results.Interval45]))])
     else
        if  strcmp(seq_results(1,1).flag_tipo_norm,'Z')
            ylim([-5 5])
            yline(0)
        else
            ylim([0 1]) 
        end
        
    end
    xlim([0 length(IKI_trial_corr)]);
    title(['Curva de Learning - ' titulo])
end %end if datos corregidos