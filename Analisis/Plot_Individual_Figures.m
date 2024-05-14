function seq_results=Plot_Individual_Figures(seq_results,path,fname)

%En caso de que la data sea corregida por filtrado o normalización
%  20/5/2023
if seq_results(1,1).flag_norm==0
    titulo_norm='sin norm';
else
    titulo_norm=seq_results(1,1).flag_tipo_norm;
end
if seq_results(1,1).flag_filt ==1
    titulo_filt='filt';
else
    titulo_filt='sin filt';
end
titulo=[titulo_filt ' - ' titulo_norm];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                             FIGURE 1:                               %%%
%%% Block Duraion: Time in between start of green cross and start of the%%%
%%% red one.                                                            %%%
%%%                                                                     %%%
%%% Correct Sequence Duration: mean value of the difference of last     %%% 
%%% and first key presses of correct sequences, per block               %%%
%%%                                                                     %%%
%%% Number of correct sequences: sum of correct sequences per block     %%%
%%%                                                                     %%%
%%% Interkey Interval: mean value of interkey intervals of correct      %%%
%%% sequences per block.                                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
figure; set(gcf,'Color','white'); box OFF; hold on;

%% Plot Block duration
subplot(2,2,1); 
% errorbar(1:1:length(seq_results.GOduration),seq_results.GOduration, seq_results.standard);
plot(1:1:length(seq_results.GOduration),seq_results.GOduration,'+');
hold on;
xlabel('Blocks','FontName','Arial','FontSize',12);
ylabel('Block duration','FontName','Arial','FontSize',12); 
ylim([0 (max(seq_results.GOduration))+min(seq_results.GOduration)]);
xlim([0 length(seq_results.GOduration)+1]);

%% Plot correct sequences duration
subplot(2,2,2); 
errorbar(1:1:length(seq_results.SEQduration),seq_results.SEQduration,seq_results.SEQstandard); 
hold on;
xlabel('Blocks','FontName','Arial','FontSize',12); 
ylabel('Correct Sequences duration','FontName','Arial','FontSize',12); 
ylim([0 (max(seq_results.SEQduration)+min(seq_results.SEQduration))]);
xlim([0 length(seq_results.SEQduration)+1]);

%% Number of correct sequences
subplot(2,2,3);  plot(seq_results.correct,'bo','MarkerSize', 6);
xlabel('Blocks','FontName','Arial','FontSize',12);
ylabel('Correct Sequences','FontName','Arial','FontSize',12);
ylim([0 (max(seq_results.correct))+3]);
xlim([0 length(seq_results.correct)+1]);
%% Plot Interkeys interval per block --- HABRIA QUE PONER UN GRÀFICO DE ESTE FILTRADO
subplot(2,2,4)
plot(seq_results.Intervalmean,'k.','MarkerSize', 12);
xlabel('Blocks','FontName','Arial','FontSize',12)
ylabel('Interkeys interval','FontName','Arial','FontSize',12);
ylim([0 (max([seq_results.Interval12, seq_results.Interval23, ...
              seq_results.Interval34, seq_results.Interval45]))])
xlim([0 length(seq_results.Interval12)+1]);

%% save
saveas(gcf,[path fname(1:end-4) '_CorrectSeq_Duration.' 'fig']);
saveas(gcf,[path fname(1:end-4) '_CorrectSeq_Duration.' 'png']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             FIGURE 2:                               %%%
%%% IKI_per_trial: Interkey Interval on each correct sequence. This type%%%
%%% of plot is also known to the lab as "the learning curve"            %%%
%%% It only plots the interval of correct sequences one after another.  %%%
%%%                                                                     %%%
%%% It's also included the filtered version of IKI_per_trial            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InterKeyInterval_Plot_Individual(seq_results,titulo)
% save
saveas(gcf,[path fname(1:end-4) '_CurvaLearning.' 'fig']);
saveas(gcf,[path fname(1:end-4) '_CurvaLearning.' 'png']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                        FIGURE 3:                                    %%%
%%%                                                                     %%%
%%%               Plots Micro Gains as raw data.                        %%%
%%%                                                                     %%%
%%% Acum MOGS: acumulative sum of Micro Offline Gains per block         %%%
%%% Acum MONGS: acumulative sum of Micro Online Gains per block         %%%
%%% Acum TL: acumulative sum of Total Learning per block                %%%
%%% Not Acum MOGS: Micro Offline Gains per block                        %%%
%%% Not Acum MONGS: Micro Online Gains per block                        %%%
%%% Not Acum TL: Total Learning per block                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MOGS acumulado
seq_results.MOGS_acumulativo=cumsum(seq_results.MOGS,'omitnan');
% MONGS acumulado
seq_results.MONGS_acumulativo=cumsum(seq_results.MONGS,'omitnan');
% TOTAL LEARNING acumulado
seq_results.Total_Learning_acumulativo=cumsum(seq_results.Total_Learning,'omitnan');


MicroGains_Plot_Individual(seq_results.MOGS,seq_results.MONGS,seq_results.Total_Learning,...
    seq_results.MOGS_acumulativo,seq_results.MONGS_acumulativo,seq_results.Total_Learning_acumulativo,['Micro Gains - crudo'],seq_results)

% Save
saveas(gcf,[path fname(1:end-4) '_MicroGains.' 'fig']);
saveas(gcf,[path fname(1:end-4) '_MicroGains.' 'png']);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                        FIGURE 4:                                    %%%
%%%                                                                     %%%
%%%               Plots Micro Gains as filtered data.                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if seq_results(1,1).flag_norm==1 || seq_results(1,1).flag_filt ==1

    % MOGS acumulado
    seq_results.MOGS_acumulativo_corr=cumsum(seq_results.MOGS_corr,'omitnan');
    % MONGS acumulado
    seq_results.MONGS_acumulativo_corr=cumsum(seq_results.MONGS_corr,'omitnan');
    % TOTAL LEARNING
    seq_results.Total_Learning_acumulativo_corr=cumsum(seq_results.Total_Learning_corr,'omitnan');
   
    MicroGains_Plot_Individual(seq_results.MOGS_corr,seq_results.MONGS_corr,seq_results.Total_Learning_corr,...
        seq_results.MOGS_acumulativo_corr,seq_results.MONGS_acumulativo_corr,seq_results.Total_Learning_acumulativo_corr,['Micro Gains - ' titulo],seq_results)

    % Save figure: .fig and .png
    saveas(gcf,[path fname(1:end-4) '_MicroGains_corr.' 'fig']);
    saveas(gcf,[path fname(1:end-4) '_MicroGains_corr.' 'png']);
end % end if data_corr
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                            FIGURE 5                                 %%%
%%%                  Micro Micro Gains as Raw data                      %%%
%%%                                                                     %%%
%%% MicroMogs Not Acum: Micro offline gains per trial                   %%%
%%% MicroMongs Not Acum: Micro online gains per trial                   %%%
%%% MicroMogs Acum: sum of Micro Offline gains per trial                %%%
%%% MicroMongs Acum: sum of Micro Online gains per trial                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[seq_results.MicroMogs_acum,seq_results.MicroMongs_acum,seq_results.MicroTL_acum]= MicroMicroGains_Plot_Individual(seq_results.MicroMOGS,seq_results.MicroMONGS,seq_results.MicroTL,['Micro Micro Gains - crudo'],seq_results);
saveas(gcf,[path fname(1:end-4) '_MicroMicroGains_crudo.' 'fig']);
saveas(gcf,[path fname(1:end-4) '_MicroMicroGains_crudo.' 'png']);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                            FIGURE 6                                 %%%
%%%             Micro Micro Gains as Filtered and/or Normalized data    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if seq_results(1,1).flag_norm==1 || seq_results(1,1).flag_filt ==1
    [seq_results.MicroMogs_corr_acum,seq_results.MicroMongs_corr_acum,seq_results.MicroTL_corr_acum]= MicroMicroGains_Plot_Individual(seq_results.MicroMOGS_corr,seq_results.MicroMONGS_corr,seq_results.MicroTL_corr,['Micro Micro Gains - ' titulo ],seq_results);

    saveas(gcf,[path fname(1:end-4) '_MicroMicroGains_corr.' 'fig']);
    saveas(gcf,[path fname(1:end-4) '_MicroMicroGains_corr.' 'png']);
end %end if data_corr


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             FIGURE 8                                %%%
%%% Mean of micro micro gains per block: Based on figures 5 and 6,      %%% 
%%% mean value is calculated for each block. This meassure allows to    %%%
%%% see if the micro micro gain of the whole block is positive or not.  %%%
%%% The mean value is calculated for the acumulative and non acumulative
%%% gains if the data is not normalized. Otherwise is only calculated   %%%
%%% for the non acumulative gains.                                      %%%
%%%                                                                     %%%
%%%                             29/12/22                                %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MOGS NO ACUM
seq_results(1,1).media_MicroMogs=NaN(1,size(seq_results.MicroMOGS,1));
seq_results(1,1).mediana_MicroMogs=NaN(1,size(seq_results.MicroMOGS,1));
seq_results(1,1).std_MicroMogs=NaN(1,size(seq_results.MicroMOGS,1));

for i=1:size(seq_results.MicroMOGS,1)
    seq_results(1,1).media_MicroMogs(i)=nanmean(seq_results.MicroMOGS(i,:));
    seq_results(1,1).mediana_MicroMogs(i)=nanmedian(seq_results.MicroMOGS(i,:));
    seq_results(1,1).std_MicroMogs(i)=nanstd(seq_results.MicroMOGS(i,:));
end


%MONGS NO ACUM
seq_results(1,1).media_MicroMongs=NaN(1,size(seq_results.MicroMONGS,1));
seq_results(1,1).mediana_MicroMongs=NaN(1,size(seq_results.MicroMONGS,1));
seq_results(1,1).std_MicroMongs=NaN(1,size(seq_results.MicroMONGS,1));

for i=1:size(seq_results.MicroMONGS,1)
    seq_results(1,1).media_MicroMongs(i)=nanmean(seq_results.MicroMONGS(i,:));
    seq_results(1,1).mediana_MicroMongs(i)=nanmedian(seq_results.MicroMONGS(i,:));
    seq_results(1,1).std_MicroMongs(i)=nanstd(seq_results.MicroMONGS(i,:));
end

%TL NO ACUM
seq_results(1,1).media_MicroTL=NaN(1,size(seq_results.MicroTL,1));
seq_results(1,1).mediana_MicroTL=NaN(1,size(seq_results.MicroTL,1));
seq_results(1,1).std_MicroTL=NaN(1,size(seq_results.MicroTL,1));

for i=1:size(seq_results.MicroTL,1)
    seq_results(1,1).media_MicroTL(i)=nanmean(seq_results.MicroTL(i,:));
    seq_results(1,1).mediana_MicroTL(i)=nanmedian(seq_results.MicroTL(i,:));
    seq_results(1,1).std_MicroTL(i)=nanstd(seq_results.MicroTL(i,:));
end

[seq_results(1,1).media_MicroMogs_acum,seq_results(1,1).media_MicroMongs_acum,seq_results(1,1).media_MicroTL_acum]=MeanMicroMicro_Gains_Plot_Indivudual(seq_results(1,1).media_MicroMogs,seq_results(1,1).media_MicroMongs,seq_results(1,1).media_MicroTL,seq_results(1,1).std_MicroMogs,seq_results(1,1).std_MicroMongs,seq_results(1,1).std_MicroTL,'Mean');

saveas(gcf,[path fname(1:end-4) ['_MicroMicroGains_crudo_Mean' '.'] 'fig']);
saveas(gcf,[path fname(1:end-4) ['_MicroMicroGains_curdo_Mean' '.'] 'png']);

% 6/9/23
[seq_results(1,1).mediana_MicroMogs_acum,seq_results(1,1).mediana_MicroMongs_acum,seq_results(1,1).mediana_MicroTL_acum]=MeanMicroMicro_Gains_Plot_Indivudual(seq_results(1,1).mediana_MicroMogs,seq_results(1,1).mediana_MicroMongs,seq_results(1,1).mediana_MicroTL,seq_results(1,1).std_MicroMogs,seq_results(1,1).std_MicroMongs,seq_results(1,1).std_MicroTL,'Median');

saveas(gcf,[path fname(1:end-4) ['_MicroMicroGains_crudo_Median' '.'] 'fig']);
saveas(gcf,[path fname(1:end-4) ['_MicroMicroGains_curdo_Median' '.'] 'png']);

if seq_results(1,1).flag_norm==1 || seq_results(1,1).flag_filt ==1
    %MOGS NO ACUM
    seq_results(1,1).media_MicroMogs_corr=NaN(1,size(seq_results.MicroMOGS_corr,1));
    seq_results(1,1).mediana_MicroMogs_corr=NaN(1,size(seq_results.MicroMOGS_corr,1));
    seq_results(1,1).std_MicroMogs_corr=NaN(1,size(seq_results.MicroMOGS_corr,1));

    for i=1:size(seq_results.MicroMOGS,1)
        seq_results(1,1).media_MicroMogs_corr(i)=nanmean(seq_results.MicroMOGS_corr(i,:));
        seq_results(1,1).mediana_MicroMogs_corr(i)=nanmedian(seq_results.MicroMOGS_corr(i,:));
        seq_results(1,1).std_MicroMogs_corr(i)=nanstd(seq_results.MicroMOGS_corr(i,:));
    end


    %MONGS NO ACUM
    seq_results(1,1).media_MicroMongs_corr=NaN(1,size(seq_results.MicroMONGS_corr,1));
    seq_results(1,1).mediana_MicroMongs_corr=NaN(1,size(seq_results.MicroMONGS_corr,1));
    seq_results(1,1).std_MicroMongs_corr=NaN(1,size(seq_results.MicroMONGS_corr,1));

    for i=1:size(seq_results.MicroMONGS_corr,1)
        seq_results(1,1).media_MicroMongs_corr(i)=nanmean(seq_results.MicroMONGS_corr(i,:));
        seq_results(1,1).mediana_MicroMongs_corr(i)=nanmedian(seq_results.MicroMONGS_corr(i,:));
        seq_results(1,1).std_MicroMongs_corr(i)=nanstd(seq_results.MicroMONGS_corr(i,:));
    end

    %TL NO ACUM
    seq_results(1,1).media_MicroTL_corr=NaN(1,size(seq_results.MicroTL_corr,1));
    seq_results(1,1).mediana_MicroTL_corr=NaN(1,size(seq_results.MicroTL_corr,1));
    seq_results(1,1).std_MicroTL_corr=NaN(1,size(seq_results.MicroTL_corr,1));

    for i=1:size(seq_results.MicroTL_corr,1)
        seq_results(1,1).media_MicroTL_corr(i)=nanmean(seq_results.MicroTL_corr(i,:));
        seq_results(1,1).mediana_MicroTL_corr(i)=nanmedian(seq_results.MicroTL_corr(i,:));
        seq_results(1,1).std_MicroTL_corr(i)=nanstd(seq_results.MicroTL_corr(i,:));
    end
    
    
    [seq_results(1,1).media_MicroMogs_corr_acum,seq_results(1,1).media_MicroMongs_corr_acum,seq_results(1,1).media_MicroTL_corr_acum]=MeanMicroMicro_Gains_Plot_Indivudual(seq_results(1,1).media_MicroMogs_corr,seq_results(1,1).media_MicroMongs_corr,seq_results(1,1).media_MicroTL_corr,seq_results(1,1).std_MicroMogs_corr,seq_results(1,1).std_MicroMongs_corr,seq_results(1,1).std_MicroTL_corr,'Mean');

    saveas(gcf,[path fname(1:end-4) ['_MicroMicroGains_corr_Mean' '.'] 'fig']);
    saveas(gcf,[path fname(1:end-4) ['_MicroMicroGains_corr_Mean' '.'] 'png']);

    % 6/9/23
    [seq_results(1,1).mediana_MicroMogs_corr_acum,seq_results(1,1).mediana_MicroMongs_corr_acum,seq_results(1,1).mediana_MicroTL_corr_acum]=MeanMicroMicro_Gains_Plot_Indivudual(seq_results(1,1).mediana_MicroMogs_corr,seq_results(1,1).mediana_MicroMongs_corr,seq_results(1,1).mediana_MicroTL_corr,seq_results(1,1).std_MicroMogs,seq_results(1,1).std_MicroMongs_corr,seq_results(1,1).std_MicroTL_corr,'Median');

    saveas(gcf,[path fname(1:end-4) ['_MicroMicroGains_corr_Median' '.'] 'fig']);
    saveas(gcf,[path fname(1:end-4) ['_MicroMicroGains_corr_Median' '.'] 'png']);
    
end

% MOGS ACUM
% 
% plot(seq_results(1,1).mediana_MicroMogs_acum,'r.')
% errorbar(seq_results(1,1).mediana_MicroMogs_acum,seq_results(1,1).std_MicroMogs_acum,'r')
% yline(0)
% ylabel('MicroMogs acum [s]')
% xlabel('Block')
% ylim([-3 3])
% title('acum')
% 
% % MONGS ACUM
% subplot(2,2,4)
% plot(seq_results(1,1).mediana_MicroMongs_acum,'b.')
% errorbar(seq_results(1,1).mediana_MicroMongs_acum,seq_results(1,1).std_MicroMongs_acum,'b')
% yline(0)
% ylabel('MicroMongs acum [s]')
% xlabel('Block')
% ylim([-3 3])
% title('acum')


% save
% saveas(gcf,[path fname(1:end-4) ['_MicroMicroGains_crudo_Mean' '.'] 'fig']);
% saveas(gcf,[path fname(1:end-4) ['_MicroMicroGains_curdo_Mean' '.'] 'png']);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     Learning Curve visualization                      %%%
%%%                          23/6/23                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seq_results.IKI_visual=Visualization_LearningCurve(seq_results(1,1).IKI_per_trial,seq_results(1,1).seqduration,seq_results);
saveas(gcf,[path fname(1:end-4) '_CurvaLearningVisual.' 'fig']);
saveas(gcf,[path fname(1:end-4) '_CurvaLearningVisual.' 'png']);

seq_results.IKI_visual_corr=Visualization_LearningCurve(seq_results(1,1).IKI_per_trial_corr,seq_results(1,1).seqduration,seq_results);
saveas(gcf,[path fname(1:end-4) '_CurvaLearningVisual_corr.' 'fig']);
saveas(gcf,[path fname(1:end-4) '_CurvaLearningVisual_corr.' 'png']);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%              MicroMicro Gains visualization                         %%%
%%%                              6/8/23                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[seq_results.MicroMOGS_visual,seq_results.MicroMONGS_visual] = Visualization_MicroMicroGains_Individual(seq_results.MicroMOGS,seq_results.MicroMONGS);
%saveas(gcf,[path fname(1:end-4) '_MicroMicroGains_Visual.' 'fig']);
%saveas(gcf,[path fname(1:end-4) '_MicroMicroGains_Visual.' 'png']);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             FIGURE 9                                %%%
%%% Sequence Duration: This is an other representation of the learning  %%%
%%% curve in which instead of ploting the mean value of the transitions %%%
%%% on each sequence, we plot the sum of them.                          %%%
%%% Subplot 1: Per trial                                                %%%
%%% Subplot 2: Per block. Mean value of the correct sequences           %%%
%%%                                                                     %%%
%%%                            3/1/23                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(9); set(gcf,'Color','white'); box OFF; hold on;
% 
% %subplot(1,2,1)
% plot(reshape(seq_results(1,1).seqduration_trial',1,[]))
% for i=1:12:noTrials
%     xline(i); %indicates a new block 
% end

% %% save
% saveas(gcf,[path fname(1:end-4) '_CorrectSeq_Duration.' 'fig']);
% saveas(gcf,[path fname(1:end-4) '_CorrectSeq_Duration.' 'png']);
end