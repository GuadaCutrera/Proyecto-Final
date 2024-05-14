%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function is called from ID_RunAnalysisGroup2.m or              %%%
%%%                                                                     %%%
%%% Parameters:                                                         %%%
%%%    - Group_Parameters: has all the parameters calculated on each    %%%
%%%    individual (from "Analisis_Seq_Individual.m") from all subjects  %%%
%%%    being analized as a group.                                       %%%
%%%    - Path: where to store the figures                               %%%
%%%    - Titulo: how the group is called. example: "With ITI"           %%%
%%%    - Paradigm_flag: wich paradigm is being analyzed "time" or "key" %%%
%%%    - N: subjects                                                    %%%
%%%                                                                     %%%
%%% Group_Results: stores the results as mean value and standard error  %%% 
%%%                                                                     %%%
%%%                   Guadalupe Cutrera 19/12/22                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Group_Results = Plot_Grupal_Figures(Group_Parameters,path,titulo,paradigm_flag,N)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             FIGURE 1:                               %%%
%%% Block Duraion: Time in between start of green cross and start of the%%%
%%% red one.                                                            %%%
%%%                                                                     %%%
%%% Correct Sequence Duration: mean value of the difference of last     %%% 
%%% and first intertap interval of correct sequences, per block         %%%
%%%                                                                     %%%
%%% Number of correct sequences: sum of correct sequences per block     %%%
%%%                                                                     %%%
%%% Interkey Interval: mean value of interkey intervals of correct      %%%
%%% sequences per block.                                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); set(gcf,'Color','white'); box OFF; hold on; sgtitle(titulo)
set(gcf,'Position',get(0,'ScreenSize'));

subplot(2,2,1); 
Group_Results.duration_group_mean = nanmean(Group_Parameters.duration_group,1);
Group_Results.duration_group_std = nanstd(Group_Parameters.duration_group,0,1)/sqrt(N);

plot(1:1:length(Group_Results.duration_group_mean),Group_Results.duration_group_mean,'+');
errorbar(1:1:length(Group_Results.duration_group_mean),Group_Results.duration_group_mean, Group_Results.duration_group_std);
hold on;
xlabel('Blocks','FontName','Arial','FontSize',12);
ylabel('Block duration','FontName','Arial','FontSize',12); 
ylim([0 (max(Group_Results.duration_group_mean)+min(Group_Results.duration_group_mean))]);
xlim([0 length(Group_Results.duration_group_mean)+1]);

%% Plot correct sequences duration
Group_Results.seq_group_mean = nanmean(Group_Parameters.seq_group,1);
Group_Results.seq_group_std = nanstd(Group_Parameters.seq_group,0,1)/sqrt(N);

subplot(2,2,2); 
errorbar(1:1:length(Group_Results.seq_group_mean),Group_Results.seq_group_mean,Group_Results.seq_group_std); 
hold on;
xlabel('Blocks','FontName','Arial','FontSize',12); 
ylabel('Correct Sequences duration','FontName','Arial','FontSize',12); 
ylim([0 (max(Group_Results.seq_group_mean)+min(Group_Results.seq_group_mean))]);
xlim([0 length(Group_Results.seq_group_mean)+1]);

%% Number of correct sequences
Group_Results.seq_group_num_mean = nanmean(Group_Parameters.seq_group_num,1);
Group_Results.seq_group_num_std = nanstd(Group_Parameters.seq_group_num,0,1)/sqrt(N);

subplot(2,2,3);  plot(Group_Results.seq_group_num_mean,'+','MarkerSize', 6);
errorbar(1:length(Group_Results.seq_group_num_mean),Group_Results.seq_group_num_mean,Group_Results.seq_group_num_std)
xlabel('Blocks','FontName','Arial','FontSize',12);
ylabel('Correct Sequences','FontName','Arial','FontSize',12);
ylim([0 (max(Group_Results.seq_group_num_mean))+3]);
xlim([0 length(Group_Results.seq_group_num_mean)+1]);

%% Plot Interkeys interval
Group_Results.key_group_mean = nanmean(Group_Parameters.key_group,1);
Group_Results.key_group_std = nanstd(Group_Parameters.key_group,0,1)/sqrt(N);

subplot(2,2,4)
plot(Group_Results.key_group_mean,'k.','MarkerSize', 10); hold on
errorbar(1:1:length(Group_Results.key_group_mean),Group_Results.key_group_mean,Group_Results.key_group_std); 
xlabel('Blocks','FontName','Arial','FontSize',12)
ylabel('Interkeys interval','FontName','Arial','FontSize',12);
ylim([0 max(Group_Results.key_group_mean+0.1)])
xlim([0 length(Group_Results.key_group_mean)+1]);

%% save
saveas(gcf,[path 'CorrectSeq_Duration_Group.' 'fig']);
saveas(gcf,[path 'CorrectSeq_Duration_Group.' 'png']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             FIGURE 2:                               %%%
%%% IKI_per_trial: Interkey Interval on each correct sequence. This type%%%
%%% of plot is also known to the lab as "the learning curve"            %%%
%%% It only plots the mean interval of correct sequences one after another.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Group_Results.iki_group_mean = nanmean(Group_Parameters.iki_group,1);              
Group_Results.iki_group_std = nanstd(Group_Parameters.iki_group,0,1)/sqrt(N);               

% InterKeyInterval_Plot_Grupal(Group_Results.iki_group_mean,Group_Results.iki_group_std,paradigm_flag,Group_Parameters,[titulo ' - crudo']);
% 
% saveas(gcf,[path 'CurvaLearning_Group.' 'fig']);
% saveas(gcf,[path 'CurvaLearning_Group.' 'png']);
% 
% clear x_vector; clear x_vector2; clear barra; clear i; clear patch; clear cont; 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIGURE 3: It's also included the corrected version of IKI_per_trial  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Group_Parameters.flag_norm==1 || Group_Parameters.flag_filt==1
    
    Group_Results.iki_group_corr_mean = nanmean(Group_Parameters.iki_group_corr,1);              
    Group_Results.iki_group_corr_std = nanstd(Group_Parameters.iki_group_corr,0,1)/sqrt(N);               

%      InterKeyInterval_Plot_Grupal(Group_Results.iki_group_corr_mean,Group_Results.iki_group_corr_std,paradigm_flag,Group_Parameters,[titulo ' - ' Group_Parameters.titulo_analisis]);
% 
%     saveas(gcf,[path 'CurvaLearning_corr_Group.' 'fig']);
%     saveas(gcf,[path 'CurvaLearning_corr_Group.' 'png']);
    
end % end si esta corregido
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                        FIGURE 4:                                    %%%
%%%                                                                     %%%
%%%               Plots Micro Gains as raw data.                        %%%
%%%                                                                     %%%
%%% Acum MOGS: acumulative sum of Micro Offline Gains per block         %%%
%%% Acum MONGS: acumulative sum of Micro Online Gains per block         %%%
%%% Acum TL: acumulative sum of Total Learning per block                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MOGS acumulado
Group_Results.mogs_group_acumulado_mean = nanmean(Group_Parameters.mogs_acumulado_group,1);                                       
Group_Results.mogs_group_acumulado_std = nanstd(Group_Parameters.mogs_acumulado_group,0,1)/sqrt(N);    

% MONGS acumulado
Group_Results.mongs_group_acumulado_mean = nanmean(Group_Parameters.mongs_acumulado_group,1);                                       
Group_Results.mongs_group_acumulado_std = nanstd(Group_Parameters.mongs_acumulado_group,0,1)/sqrt(N);    

% TOTAL LEARNING
Group_Results.TL_group_acumulado_mean = nanmean(Group_Parameters.TL_acumulado_group,1);                                       
Group_Results.TL_group_acumulado_std = nanstd(Group_Parameters.TL_acumulado_group,0,1)/sqrt(N);  


% MicroGains_Plot_Grupal(Group_Results.mogs_group_acumulado_mean,Group_Results.mogs_group_acumulado_std, ... 
%     Group_Results.mongs_group_acumulado_mean,Group_Results.mongs_group_acumulado_std,Group_Results.TL_group_acumulado_mean, ...
%     Group_Results.TL_group_acumulado_std,Group_Parameters.mogs_acumulado_group,Group_Parameters.mongs_acumulado_group,Group_Parameters.TL_acumulado_group,[titulo ' - crudo'],Group_Parameters)
% 
% 
% saveas(gcf,[path '_MicroGains_acum_crudo_Group.' 'fig']);
% saveas(gcf,[path '_MicroGains_acum_crudo_Group.' 'png']);
% saveas(gcf,[path '_MicroGains_acum_crudo_Group.' 'svg']);
    
 %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                          FIGURE 5                                   %%%
%%%             Plot MOGs no acumulativo crudo                          %%%
%%%                                                                     %%%
%%% Not Acum MOGS: Micro Offline Gains per block                        %%%
%%% Not Acum MONGS: Micro Online Gains per block                        %%%
%%% Not Acum TL: Total Learning per block                               %%%
%%%                                                                     %%%
%%%                          7/3/23                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MOGS no acumulado
Group_Results.mogs_group_mean = nanmean(Group_Parameters.mogs_group,1);                                       
Group_Results.mogs_group_std = nanstd(Group_Parameters.mogs_group,0,1)/sqrt(N);               

% MONGS no acumulado
Group_Results.mongs_group_mean = nanmean(Group_Parameters.mongs_group,1);                                       
Group_Results.mongs_group_std = nanstd(Group_Parameters.mongs_group,0,1)/sqrt(N);    

% TOTAL LEARNING no acumulado
Group_Results.TL_group_mean = nanmean(Group_Parameters.TL_group,1);                                       
Group_Results.TL_group_std = nanstd(Group_Parameters.TL_group,0,1)/sqrt(N);  


% MicroGains_Plot_Grupal(Group_Results.mogs_group_mean,Group_Results.mogs_group_std, ... 
%     Group_Results.mongs_group_mean,Group_Results.mongs_group_std,Group_Results.TL_group_mean, ...
%     Group_Results.TL_group_std,Group_Parameters.mogs_group,Group_Parameters.mongs_group,Group_Parameters.TL_group,[titulo ' - crudo - no acum'],Group_Parameters)
% 
% 
% saveas(gcf,[path '_MicroGains_crudo_Group.' 'fig']);
% saveas(gcf,[path '_MicroGains_crudo_Group.' 'png']);

%% Violin plots                  26/12/23

Matriz_Ganancias(:,1)=Group_Parameters.TL_acumulado_group(:,end);
Matriz_Ganancias(:,2)=Group_Parameters.mogs_acumulado_group(:,end);
Matriz_Ganancias(:,3)=Group_Parameters.mongs_acumulado_group(:,end-1);
% 
% figure;set(gcf,'Color','white'); box OFF; hold on; sgtitle([titulo])
% set(gcf,'Position',get(0,'ScreenSize'));
% violin(Matriz_Ganancias,'xlabel',{'TL','MOGS','MONGS'},'facecolor',[.3 .3 .3; 1 0 0; 0 0 1],'edgecolor','k','mc','k','ylabel','Sum of deltas (s)','plotlegend',0,'medc',[]);
% plot(Matriz_Ganancias','k.','Markersize',15)
% yline(0)
% ylim([-2 2])
% 
% saveas(gcf,[path '_MicroGains_violin.' 'fig']);
% saveas(gcf,[path '_MicroGains_violin.' 'png']);
% %saveas(gcf,[path '_MicroGains_violin.' 'svg']);
% 
% clear Matriz_Ganancias;
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                        FIGURE 6:                                    %%%
%%%                                                                     %%%
%%%               Plots Micro Gains as filtered data.                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Group_Parameters.flag_norm==1 || Group_Parameters.flag_filt==1
  
    % MOGS acumulado
    Group_Results.mogs_group_acumulado_corr_mean = nanmean(Group_Parameters.mogs_acumulado_group_corr,1);                                       
    Group_Results.mogs_group_acumulado_corr_std = nanstd(Group_Parameters.mogs_acumulado_group_corr,0,1)/sqrt(N);    

    % MONGS acumulado
    Group_Results.mongs_group_acumulado_corr_mean = nanmean(Group_Parameters.mongs_acumulado_group_corr,1);                                       
    Group_Results.mongs_group_acumulado_corr_std = nanstd(Group_Parameters.mongs_acumulado_group_corr,0,1)/sqrt(N);    

    % TOTAL LEARNING
    Group_Results.TL_group_acumulado_corr_mean = nanmean(Group_Parameters.TL_acumulado_group_corr,1);                                       
    Group_Results.TL_group_acumulado_corr_std = nanstd(Group_Parameters.TL_acumulado_group_corr,0,1)/sqrt(N);  

    
    MicroGains_Plot_Grupal(Group_Results.mogs_group_acumulado_corr_mean,Group_Results.mogs_group_acumulado_corr_std, ... 
        Group_Results.mongs_group_acumulado_corr_mean,Group_Results.mongs_group_acumulado_corr_std,Group_Results.TL_group_acumulado_corr_mean, ...
        Group_Results.TL_group_acumulado_corr_std,Group_Parameters.mogs_acumulado_group_corr,Group_Parameters.mongs_acumulado_group_corr,...
        Group_Parameters.TL_acumulado_group_corr,[titulo ' - ' Group_Parameters.titulo_analisis '- acum'],Group_Parameters)


    saveas(gcf,[path '_MicroGains_corr_acum_Group.' 'fig']);
    saveas(gcf,[path '_MicroGains_corr_acum_Group.' 'png']);
    saveas(gcf,[path '_MicroGains_corr_acum_Group.' 'svg']);
    
    %% Violin plots                  26/12/23

    Matriz_Ganancias(:,1)=Group_Parameters.TL_acumulado_group_corr(:,end);
    Matriz_Ganancias(:,2)=Group_Parameters.mogs_acumulado_group_corr(:,end);
    Matriz_Ganancias(:,3)=Group_Parameters.mongs_acumulado_group_corr(:,end-1);

    figure;set(gcf,'Color','white'); box OFF; hold on; sgtitle([titulo])
    set(gcf,'Position',get(0,'ScreenSize'));
    violin(Matriz_Ganancias,'xlabel',{'TL','MOGS','MONGS'},'facecolor',[.3 .3 .3; 1 0 0; 0 0 1],'edgecolor','k','mc','k','ylabel','Sum of deltas (s)','plotlegend',0,'medc',[]);
    plot(Matriz_Ganancias','k.','Markersize',15)
    yline(0)
    ylim([-2 2])

    saveas(gcf,[path '_MicroGains_corr_violin.' 'fig']);
    saveas(gcf,[path '_MicroGains_corr_violin.' 'png']);
    saveas(gcf,[path '_MicroGains_corr_violin.' 'svg']);

end %if si esta corregido

 %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                          FIGURE 7                                   %%%
%%%             Plot MOGs no acumulativo filtrado                       %%%
%%%                                                                     %%%
%%% Not Acum MOGS: Micro Offline Gains per block                        %%%
%%% Not Acum MONGS: Micro Online Gains per block                        %%%
%%% Not Acum TL: Total Learning per block                               %%%
%%%                                                                     %%%
%%%                          7/3/23                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Group_Parameters.flag_norm==1 || Group_Parameters.flag_filt==1

    %MOGS
    Group_Results.mogs_group_corr_mean = nanmean(Group_Parameters.mogs_group_corr,1);                                       
    Group_Results.mogs_group_corr_std = nanstd(Group_Parameters.mogs_group_corr,0,1)/sqrt(N);               

    % MONGS no acumulado
    Group_Results.mongs_group_corr_mean = nanmean(Group_Parameters.mongs_group_corr,1);                                       
    Group_Results.mongs_group_corr_std = nanstd(Group_Parameters.mongs_group_corr,0,1)/sqrt(N);    

    % TOTAL LEARNING
    Group_Results.TL_group_corr_mean = nanmean(Group_Parameters.TL_group_corr,1);                                       
    Group_Results.TL_group_corr_std = nanstd(Group_Parameters.TL_group_corr,0,1)/sqrt(N);  

    
    MicroGains_Plot_Grupal(Group_Results.mogs_group_corr_mean,Group_Results.mogs_group_corr_std, ... 
        Group_Results.mongs_group_corr_mean,Group_Results.mongs_group_corr_std,Group_Results.TL_group_corr_mean, ...
        Group_Results.TL_group_corr_std,Group_Parameters.mogs_group_corr,Group_Parameters.mongs_group_corr,...
        Group_Parameters.TL_group_corr,[titulo ' - ' Group_Parameters.titulo_analisis '- no acum'],Group_Parameters)


    saveas(gcf,[path '_MicroGains_corr_Group.' 'fig']);
    saveas(gcf,[path '_MicroGains_corr_Group.' 'png']);
    
end % if si esta corregido

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                            FIGURE 8                                 %%%
%%%                  Micro Micro Gains as Raw data                      %%%
%%%                                                                     %%%
%%% MicroMogs Not Acum: Micro offline gains per trial                   %%%
%%% MicroMongs Not Acum: Micro online gains per trial                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MICRO MOGS no acumulado
Group_Results.micro_mogs_group_mean=nanmean(Group_Parameters.micro_mogs_group);
Group_Results.micro_mogs_group_std=nanstd(Group_Parameters.micro_mogs_group)/sqrt(N); %error estandar

% MICRO MONGS no acumulado
Group_Results.micro_mongs_group_mean=nanmean(Group_Parameters.micro_mongs_group);
Group_Results.micro_mongs_group_std=nanstd(Group_Parameters.micro_mongs_group)/sqrt(N); %error estandar

% MICRO TL no acumulado
Group_Results.micro_TL_group_mean=nanmean(Group_Parameters.micro_TL_group);
Group_Results.micro_TL_group_std=nanstd(Group_Parameters.micro_TL_group)/sqrt(N); %error estandar

% MicroMicroGains_Plot_Grupal(Group_Results.micro_mogs_group_mean,Group_Results.micro_mogs_group_std,...
%     Group_Results.micro_mongs_group_mean,Group_Results.micro_mongs_group_std, ...
%     Group_Parameters.micro_mogs_group,Group_Parameters.micro_mongs_group,[titulo ' - crudo - no acum'],Group_Parameters)
% 
% saveas(gcf,[path '_MicroMicroGains_crudo_Group.' 'fig']);
% saveas(gcf,[path '_MicroMicroGains_crudo_Group.' 'png']);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             FIGURE 9 -Raw                           %%%
%%% MicroMogs Acum: sum of Micro Offline gains per trial                %%%
%%% MicroMongs Acum: sum of Micro Online gains per trial                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MICRO MOGS acumulado
Group_Results.micro_mogs_acumulado_group_mean=nanmean(Group_Parameters.micro_mogs_acumulado_group);
Group_Results.micro_mogs_acumulado_group_std=nanstd(Group_Parameters.micro_mogs_acumulado_group)/sqrt(N); %error estandar


% MICRO MONGS acumulado
Group_Results.micro_mongs_acumulado_group_mean=nanmean(Group_Parameters.micro_mongs_acumulado_group);
Group_Results.micro_mongs_acumulado_group_std=nanstd(Group_Parameters.micro_mongs_acumulado_group)/size(Group_Parameters.micro_mongs_acumulado_group,1); %error estandar


%MICRO TL acumulado
Group_Results.micro_TL_acumulado_group_mean=nanmean(Group_Parameters.micro_TL_acumulado_group);
Group_Results.micro_TL_acumulado_group_std=nanstd(Group_Parameters.micro_TL_acumulado_group)/sqrt(N); %error estandar

% MicroMicroGains_Plot_Grupal(Group_Results.micro_mogs_acumulado_group_mean,Group_Results.micro_mogs_acumulado_group_std,...
%     Group_Results.micro_mongs_acumulado_group_mean,Group_Results.micro_mongs_acumulado_group_std, ...
%     Group_Parameters.micro_mogs_acumulado_group,Group_Parameters.micro_mongs_acumulado_group,[titulo ' - crudo - acum'],Group_Parameters)
% 
% 
% saveas(gcf,[path '_MicroMicroGains_acum_crudo_Group.' 'fig']);
% saveas(gcf,[path '_MicroMicroGains_acum_crudo_Group.' 'png']);


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                            FIGURE 10                                %%%
%%%             Micro Micro Gains as Filtered data                      %%%
%%%             non acumulative                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if Group_Parameters.flag_norm==1 || Group_Parameters.flag_filt==1
% % 
% %     
%     %MICRO MOGS no acumulado
%     Group_Results.micro_mogs_group_mean_corr=nanmean(Group_Parameters.micro_mogs_group_corr);
%     Group_Results.micro_mogs_group_std_corr=nanstd(Group_Parameters.micro_mogs_group_corr)/sqrt(N); %error estandar

%     % MICRO MONGS no acumulado
%     Group_Results.micro_mongs_group_mean_corr=nanmean(Group_Parameters.micro_mongs_group_corr);
%     Group_Results.micro_mongs_group_std_corr=nanstd(Group_Parameters.micro_mongs_group_corr)/sqrt(N); %error estandar

    % MicroMicroGains_Plot_Grupal(Group_Results.micro_mogs_group_mean_corr,Group_Results.micro_mogs_group_mean_corr,...
    %     Group_Results.micro_mongs_group_mean_corr,Group_Results.micro_mongs_group_std_corr, ...
    %     Group_Parameters.micro_mogs_group_corr,Group_Parameters.micro_mongs_group_corr,[titulo ' - ' Group_Parameters.titulo_analisis],Group_Parameters)
    % 
    
%     saveas(gcf,[path '_MicroMicroGains_corr_Group.' 'fig']);
%     saveas(gcf,[path '_MicroMicroGains_corr_Group.' 'png']);
% end %if si esta corregido
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             FIGURE 11                               %%%
%%%                 Micro Micro Gains as filtered data                  %%%
%%%                 Acumulative                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  if Group_Parameters.flag_norm==1 || Group_Parameters.flag_filt==1

%     %MICRO MOGS acumulado
%     Group_Results.micro_mogs_acumulado_group_mean_corr=nanmean(Group_Parameters.micro_mogs_acumulado_group_corr);
%     Group_Results.micro_mogs_acumulado_group_std_corr=nanstd(Group_Parameters.micro_mogs_acumulado_group_corr)/sqrt(N); %error estandar

%      % MICRO MONGS acumulado
%     Group_Results.micro_mongs_acumulado_group_mean_corr=nanmean(Group_Parameters.micro_mongs_acumulado_group_corr);
%     Group_Results.micro_mongs_acumulado_group_std_corr=nanstd(Group_Parameters.micro_mongs_acumulado_group_corr)/sqrt(N); %error estandar

    % MicroMicroGains_Plot_Grupal(Group_Results.micro_mogs_acumulado_group_mean_corr,Group_Results.micro_mogs_acumulado_group_std_corr,...
    %     Group_Results.micro_mongs_acumulado_group_mean_corr,Group_Results.micro_mongs_acumulado_group_std_corr, ...
    %     Group_Parameters.micro_mogs_acumulado_group_corr,Group_Parameters.micro_mongs_acumulado_group_corr,[titulo ' - ' Group_Parameters.titulo_analisis],Group_Parameters)
  
%     saveas(gcf,[path '_MicroMicroGains_acum_corr_Group.' 'fig']);
%     saveas(gcf,[path '_MicroMicroGains_acum_corr_Group.' 'png']);

%end %si esta corregido

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             Figure 13                               %%%
%%%  Micro Micro Gains mean value per block. Shown data is raw.       %%%
%%%  Non acumulative values are shown if data is normalized but acumulative
%%%  values are not.                                                    %%%
%%%                              6/1/23                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MOGS NO ACUM -------------------------------------------------------------
Group_Results.media_MicroMogs_mean=nanmean(Group_Parameters.media_MicroMOGS_group);
Group_Results.media_MicroMogs_median=nanmedian(Group_Parameters.media_MicroMOGS_group); %media por bloque - mediana across subj
Group_Results.media_MicroMogs_std=nanstd(Group_Parameters.media_MicroMOGS_group)/sqrt(N); %error estandar
for i=1:8
    Group_Results.media_MicroMogs_sem(i)=bootstrap_SEM_median(Group_Parameters.media_MicroMOGS_group(:,i));
end

%MONGS NO ACUM ------------------------------------------------------------
Group_Results.media_MicroMongs_mean=nanmean(Group_Parameters.media_MicroMONGS_group);
Group_Results.media_MicroMongs_median=nanmedian(Group_Parameters.media_MicroMONGS_group);
Group_Results.media_MicroMongs_std=nanstd(Group_Parameters.media_MicroMONGS_group)/sqrt(N); %error estandar
for i=1:8
    Group_Results.media_MicroMongs_sem(i)=bootstrap_SEM_median(Group_Parameters.media_MicroMONGS_group(:,i));
end

%TL NO ACUM ------------------------------------------------------------
Group_Results.media_MicroTL_mean=nanmean(Group_Parameters.media_MicroTL_group);
Group_Results.media_MicroTL_median=nanmedian(Group_Parameters.media_MicroTL_group);
Group_Results.media_MicroTL_std=nanstd(Group_Parameters.media_MicroTL_group)/sqrt(N); %error estandar
for i=1:8
    Group_Results.media_MicroTL_sem(i)=bootstrap_SEM_median(Group_Parameters.media_MicroTL_group(:,i));
end


%MOGS ACUM ----------------------------------------------------------------
Group_Results.media_MicroMogs_acum_mean=nanmean(Group_Parameters.media_MicroMOGS_acum_group);
Group_Results.media_MicroMogs_acum_median=nanmedian(Group_Parameters.media_MicroMOGS_acum_group);
Group_Results.media_MicroMogs_acum_std=nanstd(Group_Parameters.media_MicroMOGS_acum_group)/sqrt(N); %error estandar
for i=1:8
    Group_Results.media_MicroMogs_acum_sem(i)=bootstrap_SEM_median(Group_Parameters.media_MicroMOGS_acum_group(:,i));
end
%MONGS ACUM ---------------------------------------------------------------
Group_Results.media_MicroMongs_acum_mean=nanmean(Group_Parameters.media_MicroMONGS_acum_group);
Group_Results.media_MicroMongs_acum_median=nanmedian(Group_Parameters.media_MicroMONGS_acum_group);
Group_Results.media_MicroMongs_acum_std=nanstd(Group_Parameters.media_MicroMONGS_acum_group)/sqrt(N); %error estandar
for i=1:8
    Group_Results.media_MicroMongs_acum_sem(i)=bootstrap_SEM_median(Group_Parameters.media_MicroMONGS_acum_group(:,i));
end

%TL ACUM ---------------------------------------------------------------
Group_Results.media_MicroTL_acum_mean=nanmean(Group_Parameters.media_MicroTL_acum_group);
Group_Results.media_MicroTL_acum_median=nanmedian(Group_Parameters.media_MicroTL_acum_group);
Group_Results.media_MicroTL_acum_std=nanstd(Group_Parameters.media_MicroTL_acum_group)/sqrt(N); %error estandar
for i=1:8
    Group_Results.media_MicroTL_acum_sem(i)=bootstrap_SEM_median(Group_Parameters.media_MicroTL_acum_group(:,i));
end


%media across subj
% MeanMicroMicroGains_Plot_Grupal(Group_Results.media_MicroMongs_mean,Group_Results.media_MicroMogs_mean,...
% Group_Results.media_MicroMongs_std,Group_Results.media_MicroMogs_std,Group_Results.media_MicroMogs_acum_mean,...
% Group_Results.media_MicroMongs_acum_mean,Group_Parameters.media_MicroMOGS_acum_group,Group_Parameters.media_MicroMONGS_acum_group,titulo,Group_Parameters,'mean')

% saveas(gcf,[path '_Mean_MicroMicroGains_media_crudo.' 'fig']);
% saveas(gcf,[path '_Mean_MicroMicroGains_media_crudo.' 'png']);
% 
% %mediana across subj
% MeanMicroMicroGains_Plot_Grupal(Group_Results.media_MicroMongs_median,Group_Results.media_MicroMogs_median,...
% Group_Results.media_MicroMongs_sem,Group_Results.media_MicroMogs_sem,Group_Results.media_MicroMogs_acum_median,...
% Group_Results.media_MicroMongs_acum_median,Group_Parameters.media_MicroMOGS_acum_group,Group_Parameters.media_MicroMONGS_acum_group,titulo,Group_Parameters,'mean')
% 
% saveas(gcf,[path '_Mean_MicroMicroGains_mediana_crudo.' 'fig']);
% saveas(gcf,[path '_Mean_MicroMicroGains_mediana_crudo.' 'png']);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             Figure 13                               %%%
%%%  Same as previous figure but with median value per block            %%%
%%%                                7/9/23                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MOGS NO ACUM -------------------------------------------------------------
Group_Results.mediana_MicroMogs_mean=nanmean(Group_Parameters.mediana_MicroMOGS_group);
Group_Results.mediana_MicroMogs_median=nanmedian(Group_Parameters.mediana_MicroMOGS_group); %mediana por bloque - mediana across subj
Group_Results.mediana_MicroMogs_std=nanstd(Group_Parameters.mediana_MicroMOGS_group)/sqrt(N); %error estandar
for i=1:length(Group_Results.mediana_MicroMogs_mean)
    Group_Results.mediana_MicroMogs_sem(i)=bootstrap_SEM_median(Group_Parameters.mediana_MicroMOGS_group(:,i)); %error estandar
end
%MONGS NO ACUM ------------------------------------------------------------
Group_Results.mediana_MicroMongs_mean=nanmean(Group_Parameters.mediana_MicroMONGS_group);
Group_Results.mediana_MicroMongs_median=nanmedian(Group_Parameters.mediana_MicroMONGS_group);
Group_Results.mediana_MicroMongs_std=nanstd(Group_Parameters.mediana_MicroMONGS_group)/sqrt(N); %error estandar
for i=1:length(Group_Results.mediana_MicroMongs_mean)
    Group_Results.mediana_MicroMongs_sem(i)=bootstrap_SEM_median(Group_Parameters.mediana_MicroMONGS_group(:,i)); %error estandar
end

%TL NO ACUM ------------------------------------------------------------
Group_Results.mediana_MicroTL_mean=nanmean(Group_Parameters.mediana_MicroTL_group);
Group_Results.mediana_MicroTL_median=nanmedian(Group_Parameters.mediana_MicroTL_group);
Group_Results.mediana_MicroTL_std=nanstd(Group_Parameters.mediana_MicroTL_group)/sqrt(N); %error estandar
for i=1:length(Group_Results.mediana_MicroTL_mean)
    Group_Results.mediana_MicroTL_sem(i)=bootstrap_SEM_median(Group_Parameters.mediana_MicroTL_group(:,i)); %error estandar
end

%MOGS ACUM ----------------------------------------------------------------
Group_Results.mediana_MicroMogs_acum_mean=nanmean(Group_Parameters.mediana_MicroMOGS_acum_group);
Group_Results.mediana_MicroMogs_acum_median=nanmedian(Group_Parameters.mediana_MicroMOGS_acum_group);
Group_Results.mediana_MicroMogs_acum_std=nanstd(Group_Parameters.mediana_MicroMOGS_acum_group)/sqrt(N); %error estandar
for i=1:length(Group_Results.mediana_MicroMogs_acum_mean)
    Group_Results.mediana_MicroMogs_acum_sem(i)=bootstrap_SEM_median(Group_Parameters.mediana_MicroMOGS_acum_group(:,i)); %error estandar
end
%MONGS ACUM ---------------------------------------------------------------
Group_Results.mediana_MicroMongs_acum_mean=nanmean(Group_Parameters.mediana_MicroMONGS_acum_group);
Group_Results.mediana_MicroMongs_acum_median=nanmedian(Group_Parameters.mediana_MicroMONGS_acum_group);
Group_Results.mediana_MicroMongs_acum_std=nanstd(Group_Parameters.mediana_MicroMONGS_acum_group)/sqrt(N); %error estandar
for i=1:length(Group_Results.mediana_MicroMongs_acum_mean)
    Group_Results.mediana_MicroMoNgs_acum_sem(i)=bootstrap_SEM_median(Group_Parameters.mediana_MicroMONGS_acum_group(:,i)); %error estandar
end

%TL ACUM ---------------------------------------------------------------
Group_Results.mediana_MicroTL_acum_mean=nanmean(Group_Parameters.mediana_MicroTL_acum_group);
Group_Results.mediana_MicroTL_acum_median=nanmedian(Group_Parameters.mediana_MicroTL_acum_group);
Group_Results.mediana_MicroTL_acum_std=nanstd(Group_Parameters.mediana_MicroTL_acum_group)/sqrt(N); %error estandar
for i=1:length(Group_Results.mediana_MicroTL_acum_mean)
    Group_Results.mediana_MicroTL_acum_sem(i)=bootstrap_SEM_median(Group_Parameters.mediana_MicroTL_acum_group(:,i)); %error estandar
end

% %media across subj
% MeanMicroMicroGains_Plot_Grupal(Group_Results.mediana_MicroMongs_mean,Group_Results.mediana_MicroMogs_mean,...
% Group_Results.mediana_MicroMongs_std,Group_Results.mediana_MicroMogs_std,Group_Results.mediana_MicroMogs_acum_mean,...
% Group_Results.mediana_MicroMongs_acum_mean,Group_Parameters.mediana_MicroMOGS_acum_group,Group_Parameters.mediana_MicroMONGS_acum_group,titulo,Group_Parameters,'median')
% 
% saveas(gcf,[path '_Median_MicroMicroGains_media_crudo.' 'fig']);
% saveas(gcf,[path '_Median_MicroMicroGains_media_crudo.' 'png']);
% 
% %mediana across subj
% MeanMicroMicroGains_Plot_Grupal(Group_Results.mediana_MicroMongs_median,Group_Results.mediana_MicroMogs_median,...
% Group_Results.mediana_MicroMongs_sem,Group_Results.mediana_MicroMogs_sem,Group_Results.mediana_MicroMogs_acum_median,...
% Group_Results.mediana_MicroMongs_acum_median,Group_Parameters.mediana_MicroMOGS_acum_group,Group_Parameters.mediana_MicroMONGS_acum_group,titulo,Group_Parameters,'median')
% 
% saveas(gcf,[path '_Median_MicroMicroGains_mediana_crudo.' 'fig']);
% saveas(gcf,[path '_Median_MicroMicroGains_mediana_crudo.' 'png']);


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                            FIGURE 14 y 15                           %%%
%%%           Correlation Block Duration Vs MOGS and MONGS              %%%
%%% Because some of the subjects didn't have the expected behaviour, we %%%
%%% think that the Micro Gains may be retalted to the time each subject %%%
%%% takes to complete each block.                                       %%%
%%% The mean values are shown in the figure, but a matrix for each subject%
%%%is created.                                                          %%%
%%% In order to get the R^2 value of the correlation, we make a linear  %%%
%%% regression.
%%%                          12/5/2023                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CorrelacionPlot(Group_Results,Group_Results.mogs_group_mean,Group_Results.mongs_group_mean,Group_Results.mogs_group_acumulado_mean,Group_Results.mongs_group_acumulado_mean,[titulo ' - crudo'],paradigm_flag);
% 
% saveas(gcf,[path '_Correlacion_crudo_Group.' 'fig']);
% saveas(gcf,[path '_Correlacion_crudo_Group.' 'png']);

if Group_Parameters.flag_norm==1 || Group_Parameters.flag_filt==1
    
%    CorrelacionPlot(Group_Results,Group_Results.mogs_group_corr_mean,Group_Results.mongs_group_corr_mean,Group_Results.mogs_group_acumulado_corr_mean,...
%        Group_Results.mongs_group_acumulado_corr_mean,[titulo ' - ' Group_Parameters.titulo_analisis],paradigm_flag);
% 
%    saveas(gcf,[path '_Correlacion_corr_Group.' 'fig']);
%    saveas(gcf,[path '_Correlacion_corr_Group.' 'png']);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             FIGURE 16 y 17                          %%%
%%% Figure with 3 axis to show how Micro Gains change with the block    %%%
%%% duration in each block.                                             %%%
%%% This figure shows the first aproach with non acumulative and acumulative
%%% raw data.                                                           %%%
%%%                         12/5/2023                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BlockDur_vs_MicroGains(Group_Results,Group_Results.mogs_group_mean,Group_Results.mongs_group_mean,Group_Results.mogs_group_acumulado_mean,...
%     Group_Results.mongs_group_acumulado_mean,[titulo ' - crudo'],paradigm_flag);
% 
% saveas(gcf,[path '_BlockDur_Vs_MicroGains_crudo_Group.' 'fig']);
% saveas(gcf,[path '_BlockDur_Vs_MicroGains_crudo_Group.' 'png']);

if Group_Parameters.flag_norm==1 || Group_Parameters.flag_filt==1
    
%    BlockDur_vs_MicroGains(Group_Results,Group_Results.mogs_group_corr_mean,Group_Results.mongs_group_corr_mean,Group_Results.mogs_group_acumulado_corr_mean,...
%        Group_Results.mongs_group_acumulado_corr_mean,[titulo ' - ' Group_Parameters.titulo_analisis],paradigm_flag);
% 
%     saveas(gcf,[path '_BlockDur_Vs_MicroGains_corr_Group.' 'fig']);
%     saveas(gcf,[path '_BlockDur_Vs_MicroGains_corr_Group.' 'png']);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         FIGURE                                      %%%
%%% Tapping Speed
%%%                      11/6/2023                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Group_Results.speed_group_mean = nanmean(1./Group_Parameters.iki_group,1);              
Group_Results.speed_group_std = nanstd(1./Group_Parameters.iki_group,0,1)/sqrt(N);  

% tapping_speed(Group_Results.speed_group_mean,Group_Results.speed_group_std,Group_Parameters,paradigm_flag,[titulo ' - crudo']);
% 
% saveas(gcf,[path '_tapping_speed_crudo_Group.' 'fig']);
% saveas(gcf,[path '_tapping_speed_crudo_Group.' 'png']);

if Group_Parameters.flag_norm==1 || Group_Parameters.flag_filt==1
    
    Group_Results.speed_group_corr_mean = nanmean(1./Group_Parameters.iki_group_corr,1);              
    Group_Results.speed_group_corr_std = nanstd(1./Group_Parameters.iki_group_corr,0,1)/sqrt(N);
 %   tapping_speed(Group_Results.speed_group_corr_mean,Group_Results.speed_group_corr_std,Group_Parameters,paradigm_flag,[titulo ' - ' Group_Parameters.titulo_analisis]);
    
  %  saveas(gcf,[path '_tapping_speed_corr_Group.' 'fig']);
  %  saveas(gcf,[path '_tapping_speed_corr_Group.' 'png']);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         FIGURE                                      %%%
%%% Learning Curve and Tapping Speed, visualization                     %%%
%%%                                                                     %%%
%%%                      23/6/2023                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------ IKI ---------------------------------------
Group_Results.iki_visual_group_mean = nanmean(Group_Parameters.iki_visual_group,1);    
Group_Results.iki_visual_group_median = nanmedian(Group_Parameters.iki_visual_group,1); 
Group_Results.iki_visual_group_std = nanstd(Group_Parameters.iki_visual_group,0,1)/sqrt(N);  

% Visualization_Learning_Grupal(Group_Results.iki_visual_group_mean,Group_Results.iki_visual_group_std,'iki',[titulo ' - crudo']);
% saveas(gcf,[path '_Learning_visual_crudo_Group.' 'fig']);
% saveas(gcf,[path '_Learning_visual_crudo_Group.' 'png']);

%Learning_print(Group_Results.iki_visual_group_mean,Group_Parameters.iki_visual_group,[path 'Sujetos\'],'crudo');

 if Group_Parameters.flag_norm==1 || Group_Parameters.flag_filt==1
    Group_Results.iki_visual_group_corr_mean = nanmean(Group_Parameters.iki_visual_group_corr,1);              
    Group_Results.iki_visual_group_corr_std = nanstd(Group_Parameters.iki_visual_group_corr,0,1)/sqrt(N);  

   Visualization_Learning_Grupal(Group_Results.iki_visual_group_corr_mean,Group_Results.iki_visual_group_corr_std,'iki',[titulo ' - ' Group_Parameters.titulo_analisis]);
    saveas(gcf,[path '_Learning_visual_corr_Group.' 'fig']);
    saveas(gcf,[path '_Learning_visual_corr_Group.' 'png']);
    
    %yyaxis left;
    sgtitle('learning inverted')
    ax = gca;
    ax.YDir = 'reverse';
    saveas(gcf,[path ['_Learning_visual_inverted_corr_Group_' titulo '.'] 'fig']);
    saveas(gcf,[path ['_Learning_visual_inverted_corr_Group_' titulo '.'] 'png']);
    saveas(gcf,[path ['_Learning_visual_inverted_corr_Group_' titulo '.'] 'svg']);
   % Learning_print(Group_Results.iki_visual_group_corr_mean,Group_Parameters.iki_visual_group_corr,[path 'Sujetos\'],'corr');
     
 end

%------------------------------ TAPPING SPEED ------------------------------
Group_Results.speed_visual_group_mean = nanmean(1./Group_Parameters.iki_visual_group,1);              
Group_Results.speed_visual_group_std = nanstd(1./Group_Parameters.iki_visual_group,0,1)/sqrt(N);  

% Visualization_Learning_Grupal(Group_Results.speed_visual_group_mean,Group_Results.speed_visual_group_std,'speed',[titulo ' - crudo']);
% 
% saveas(gcf,[path '_tapping_speed_visual_crudo_Group.' 'fig']);
% saveas(gcf,[path '_tapping_speed_visual_crudo_Group.' 'png']);

if Group_Parameters.flag_norm==1 || Group_Parameters.flag_filt==1
    
    Group_Results.speed_visual_group_corr_mean = nanmean(1./Group_Parameters.iki_visual_group_corr,1);              
    Group_Results.speed_visual_group_corr_std = nanstd(1./Group_Parameters.iki_visual_group_corr,0,1)/sqrt(N);
    
%     Visualization_Learning_Grupal(Group_Results.speed_visual_group_corr_mean,Group_Results.speed_visual_group_corr_std,'speed',[titulo ' - ' Group_Parameters.titulo_analisis]);
%     
%     saveas(gcf,[path '_tapping_speed_visual_corr_Group.' 'fig']);
%     saveas(gcf,[path '_tapping_speed_visual_corr_Group.' 'png']);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                            as Figure 13                             %%%
%%%  Micro Micro Gains median value per block. Shown data is CORRECTED  %%%
%%%  Non acumulative values are shown if data is normalized but acumulative
%%%  values are not.                                                    %%%
%%%                              7/12/23                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MOGS NO ACUM -------------------------------------------------------------
Group_Results.media_MicroMogs_corr_mean=nanmean(Group_Parameters.media_MicroMOGS_corr_group);
Group_Results.media_MicroMogs_corr_median=nanmedian(Group_Parameters.media_MicroMOGS_corr_group); %media por bloque - mediana across subj
Group_Results.media_MicroMogs_corr_std=nanstd(Group_Parameters.media_MicroMOGS_corr_group)/sqrt(N); %error estandar
for i=1:length(Group_Results.media_MicroMogs_corr_mean)
    Group_Results.media_MicroMogs_corr_sem(i)=bootstrap_SEM_median(Group_Parameters.media_MicroMOGS_corr_group(:,i));
end
%MONGS NO ACUM ------------------------------------------------------------
Group_Results.media_MicroMongs_corr_mean=nanmean(Group_Parameters.media_MicroMONGS_corr_group);
Group_Results.media_MicroMongs_corr_median=nanmedian(Group_Parameters.media_MicroMONGS_corr_group);
Group_Results.media_MicroMongs_corr_std=nanstd(Group_Parameters.media_MicroMONGS_corr_group)/sqrt(N); %error estandar
for i=1:length(Group_Results.media_MicroMongs_corr_mean)
    Group_Results.media_MicroMongs_corr_sem(i)=bootstrap_SEM_median(Group_Parameters.media_MicroMONGS_corr_group(:,i));
end
%TL NO ACUM ------------------------------------------------------------
Group_Results.media_MicroTL_corr_mean=nanmean(Group_Parameters.media_MicroTL_corr_group);
Group_Results.media_MicroTL_corr_median=nanmedian(Group_Parameters.media_MicroTL_corr_group);
Group_Results.media_MicroTL_corr_std=nanstd(Group_Parameters.media_MicroTL_corr_group)/sqrt(N); %error estandar
for i=1:length(Group_Results.media_MicroTL_corr_mean)
    Group_Results.media_MicroTL_corr_sem(i)=bootstrap_SEM_median(Group_Parameters.media_MicroTL_corr_group(:,i));
end

%MOGS ACUM ----------------------------------------------------------------
Group_Results.media_MicroMogs_corr_acum_mean=nanmean(Group_Parameters.media_MicroMOGS_corr_acum_group);
Group_Results.media_MicroMogs_corr_acum_median=nanmedian(Group_Parameters.media_MicroMOGS_corr_acum_group);
Group_Results.media_MicroMogs_corr_acum_std=nanstd(Group_Parameters.media_MicroMOGS_corr_acum_group)/sqrt(N); %error estandar
for i=1:length(Group_Results.media_MicroMogs_corr_acum_mean)
    Group_Results.media_MicroMogs_corr_acum_sem(i)=bootstrap_SEM_median(Group_Parameters.media_MicroMOGS_corr_acum_group(:,i));
end
%MONGS ACUM ---------------------------------------------------------------
Group_Results.media_MicroMongs_corr_acum_mean=nanmean(Group_Parameters.media_MicroMONGS_corr_acum_group);
Group_Results.media_MicroMongs_corr_acum_median=nanmedian(Group_Parameters.media_MicroMONGS_corr_acum_group);
Group_Results.media_MicroMongs_corr_acum_std=nanstd(Group_Parameters.media_MicroMONGS_corr_acum_group)/sqrt(N); %error estandar
for i=1:length(Group_Results.media_MicroMongs_corr_acum_mean)
    Group_Results.media_MicroMongs_corr_acum_sem(:,i)=bootstrap_SEM_median(Group_Parameters.media_MicroMONGS_corr_acum_group(:,i));
end
%TL ACUM ---------------------------------------------------------------
Group_Results.media_MicroTL_corr_acum_mean=nanmean(Group_Parameters.media_MicroTL_corr_acum_group);
Group_Results.media_MicroTL_corr_acum_median=nanmedian(Group_Parameters.media_MicroTL_corr_acum_group);
Group_Results.media_MicroTL_corr_acum_std=nanstd(Group_Parameters.media_MicroTL_corr_acum_group)/sqrt(N); %error estandar
for i=1:length(Group_Results.media_MicroTL_corr_acum_mean)
    Group_Results.media_MicroTL_corr_acum_sem(:,i)=bootstrap_SEM_median(Group_Parameters.media_MicroTL_corr_acum_group(:,i));
end

%media across subj
MeanMicroMicroGains_Plot_Grupal(Group_Results.media_MicroMongs_corr_mean,Group_Results.media_MicroMogs_corr_mean,Group_Results.media_MicroTL_corr_mean,...
Group_Results.media_MicroMongs_corr_std,Group_Results.media_MicroMogs_corr_std,Group_Results.media_MicroTL_corr_std,Group_Results.media_MicroMogs_corr_acum_mean,...
Group_Results.media_MicroMongs_corr_acum_mean,Group_Results.media_MicroTL_corr_acum_mean,Group_Parameters.media_MicroMOGS_corr_acum_group,Group_Parameters.media_MicroMONGS_corr_acum_group,Group_Parameters.media_MicroTL_corr_acum_group,titulo,Group_Parameters,'mean')

saveas(gcf,[path '_Mean_MicroMicroGains_media_corr.' 'fig']);
saveas(gcf,[path '_Mean_MicroMicroGains_media_corr.' 'png']);

% %mediana across subj
% MeanMicroMicroGains_Plot_Grupal(Group_Results.media_MicroMongs_corr_median,Group_Results.media_MicroMogs_corr_median,...
% Group_Results.media_MicroMongs_corr_sem,Group_Results.media_MicroMogs_corr_sem,Group_Results.media_MicroMogs_corr_acum_median,...
% Group_Results.media_MicroMongs_corr_acum_median,Group_Parameters.media_MicroMOGS_corr_acum_group,Group_Parameters.media_MicroMONGS_corr_acum_group,titulo,Group_Parameters,'mean')
% 
% saveas(gcf,[path '_Mean_MicroMicroGains_mediana_corr.' 'fig']);
% saveas(gcf,[path '_Mean_MicroMicroGains_mediana_corr.' 'png']);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                           as  Figure 13 but corrected               %%%
%%%  Same as previous figure but with median value per block            %%%
%%%                              7/9/23                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MOGS NO ACUM -------------------------------------------------------------
Group_Results.mediana_MicroMogs_corr_mean=nanmean(Group_Parameters.mediana_MicroMOGS_corr_group);
Group_Results.mediana_MicroMogs_corr_median=nanmedian(Group_Parameters.mediana_MicroMOGS_corr_group); %mediana por bloque - mediana across subj
Group_Results.mediana_MicroMogs_corr_std=nanstd(Group_Parameters.mediana_MicroMOGS_corr_group)/sqrt(N); %error estandar
for i=1:length(Group_Results.mediana_MicroMogs_corr_mean)
    Group_Results.mediana_MicroMogs_corr_sem(i)=bootstrap_SEM_median(Group_Parameters.mediana_MicroMOGS_corr_group(:,i)); %error estandar
end
%MONGS NO ACUM ------------------------------------------------------------
Group_Results.mediana_MicroMongs_corr_mean=nanmean(Group_Parameters.mediana_MicroMONGS_corr_group);
Group_Results.mediana_MicroMongs_corr_median=nanmedian(Group_Parameters.mediana_MicroMONGS_corr_group);
Group_Results.mediana_MicroMongs_corr_std=nanstd(Group_Parameters.mediana_MicroMONGS_corr_group)/sqrt(N); %error estandar
for i=1:length(Group_Results.mediana_MicroMongs_corr_mean)
    Group_Results.mediana_MicroMongs_corr_sem(i)=bootstrap_SEM_median(Group_Parameters.mediana_MicroMONGS_corr_group(:,i)); %error estandar
end
%MOGS ACUM ----------------------------------------------------------------
Group_Results.mediana_MicroMogs_corr_acum_mean=nanmean(Group_Parameters.mediana_MicroMOGS_corr_acum_group);
Group_Results.mediana_MicroMogs_corr_acum_median=nanmedian(Group_Parameters.mediana_MicroMOGS_corr_acum_group);
Group_Results.mediana_MicroMogs_corr_acum_std=nanstd(Group_Parameters.mediana_MicroMOGS_corr_acum_group)/sqrt(N); %error estandar
for i=1:length(Group_Results.mediana_MicroMogs_corr_acum_mean)
    Group_Results.mediana_MicroMogs_corr_acum_sem(i)=bootstrap_SEM_median(Group_Parameters.mediana_MicroMOGS_corr_acum_group(:,i)); %error estandar
end
%MONGS ACUM ---------------------------------------------------------------
Group_Results.mediana_MicroMongs_corr_acum_mean=nanmean(Group_Parameters.mediana_MicroMONGS_corr_acum_group);
Group_Results.mediana_MicroMongs_corr_acum_median=nanmedian(Group_Parameters.mediana_MicroMONGS_corr_acum_group);
Group_Results.mediana_MicroMongs_corr_acum_std=nanstd(Group_Parameters.mediana_MicroMONGS_corr_acum_group)/sqrt(N); %error estandar
for i=1:length(Group_Results.mediana_MicroMongs_corr_acum_mean)
    Group_Results.mediana_MicroMongs_corr_acum_sem(i)=bootstrap_SEM_median(Group_Parameters.mediana_MicroMONGS_corr_acum_group(:,i)); %error estandar
end
% % %media across subj
% MeanMicroMicroGains_Plot_Grupal(Group_Results.mediana_MicroMongs_corr_mean,Group_Results.mediana_MicroMogs_corr_mean,...
% Group_Results.mediana_MicroMongs_corr_std,Group_Results.mediana_MicroMogs_corr_std,Group_Results.mediana_MicroMogs_corr_acum_mean,...
% Group_Results.mediana_MicroMongs_corr_acum_mean,Group_Parameters.mediana_MicroMOGS_corr_acum_group,Group_Parameters.mediana_MicroMONGS_corr_acum_group,titulo,Group_Parameters,'median')
% 
% saveas(gcf,[path '_Median_MicroMicroGains_media_corr.' 'fig']);
% saveas(gcf,[path '_Median_MicroMicroGains_media_corr.' 'png']);
% 
% %mediana across subj
% MeanMicroMicroGains_Plot_Grupal(Group_Results.mediana_MicroMongs_corr_median,Group_Results.mediana_MicroMogs_corr_median,...
% Group_Results.mediana_MicroMongs_corr_sem,Group_Results.mediana_MicroMogs_corr_sem,Group_Results.mediana_MicroMogs_corr_acum_median,...
% Group_Results.mediana_MicroMongs_corr_acum_median,Group_Parameters.mediana_MicroMOGS_corr_acum_group,Group_Parameters.mediana_MicroMONGS_corr_acum_group,titulo,Group_Parameters,'median')
% 
% saveas(gcf,[path '_Median_MicroMicroGains_mediana_corr.' 'fig']);
% saveas(gcf,[path '_Median_MicroMicroGains_mediana_corr.' 'png']);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIGURE: PLOT INTERVALOS. Se plotean la media/mediana por bloque para
%%% cada uno de los cuatro intervalos.                                  %%%
%%%                      16/12/23                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mediana original
Group_Results.Int12_MedianPerBlock_group_mean=nanmean(Group_Parameters.Int12_MedianPerBlock_group,1);
Group_Results.Int12_MedianPerBlock_group_std=nanstd(Group_Parameters.Int12_MedianPerBlock_group)/sqrt(N);
Group_Results.Int23_MedianPerBlock_group_mean=nanmean(Group_Parameters.Int23_MedianPerBlock_group,1);
Group_Results.Int23_MedianPerBlock_group_std=nanstd(Group_Parameters.Int23_MedianPerBlock_group)/sqrt(N);
Group_Results.Int34_MedianPerBlock_group_mean=nanmean(Group_Parameters.Int34_MedianPerBlock_group,1);
Group_Results.Int34_MedianPerBlock_group_std=nanstd(Group_Parameters.Int34_MedianPerBlock_group)/sqrt(N);
Group_Results.Int45_MedianPerBlock_group_mean=nanmean(Group_Parameters.Int45_MedianPerBlock_group,1);
Group_Results.Int45_MedianPerBlock_group_std=nanstd(Group_Parameters.Int45_MedianPerBlock_group)/sqrt(N);

% plot_intervalos(Group_Parameters.Int12_MedianPerBlock_group,Group_Parameters.Int23_MedianPerBlock_group,Group_Parameters.Int34_MedianPerBlock_group,...
%     Group_Parameters.Int45_MedianPerBlock_group,Group_Results.Int12_MedianPerBlock_group_mean,Group_Results.Int23_MedianPerBlock_group_mean,...
%     Group_Results.Int34_MedianPerBlock_group_mean,Group_Results.Int45_MedianPerBlock_group_mean,Group_Results.Int12_MedianPerBlock_group_std,...
%     Group_Results.Int23_MedianPerBlock_group_std,Group_Results.Int34_MedianPerBlock_group_std,Group_Results.Int45_MedianPerBlock_group_std,[titulo '- mediana por bloque - original'])
% 
% gcf22=figure(22);
% saveas(gcf22,[path '_Intervals_MedianPerBlock.' 'fig']);
% saveas(gcf22,[path '_Intervals_MedianPerBlock.' 'png']);
% gcf23=figure(23);
% saveas(gcf23,[path '_Intervals_MedianPerBlock_completo.' 'fig']);
% saveas(gcf23,[path '_Intervals_MedianPerBlock_completo.' 'png']);

%mediana corregida
Group_Results.Int12_MedianPerBlock_group_corr_mean=nanmean(Group_Parameters.Int12_MedianPerBlock_corr_group,1);
Group_Results.Int12_MedianPerBlock_group_corr_std=nanstd(Group_Parameters.Int12_MedianPerBlock_corr_group)/sqrt(N);
Group_Results.Int23_MedianPerBlock_group_corr_mean=nanmean(Group_Parameters.Int23_MedianPerBlock_corr_group,1);
Group_Results.Int23_MedianPerBlock_group_corr_std=nanstd(Group_Parameters.Int23_MedianPerBlock_corr_group)/sqrt(N);
Group_Results.Int34_MedianPerBlock_group_corr_mean=nanmean(Group_Parameters.Int34_MedianPerBlock_corr_group,1);
Group_Results.Int34_MedianPerBlock_group_corr_std=nanstd(Group_Parameters.Int34_MedianPerBlock_corr_group)/sqrt(N);
Group_Results.Int45_MedianPerBlock_group_corr_mean=nanmean(Group_Parameters.Int45_MedianPerBlock_corr_group,1);
Group_Results.Int45_MedianPerBlock_group_corr_std=nanstd(Group_Parameters.Int45_MedianPerBlock_corr_group)/sqrt(N);

% plot_intervalos(Group_Parameters.Int12_MedianPerBlock_corr_group,Group_Parameters.Int23_MedianPerBlock_corr_group,Group_Parameters.Int34_MedianPerBlock_corr_group,...
%     Group_Parameters.Int45_MedianPerBlock_corr_group,Group_Results.Int12_MedianPerBlock_group_corr_mean,Group_Results.Int23_MedianPerBlock_group_corr_mean,...
%     Group_Results.Int34_MedianPerBlock_group_corr_mean,Group_Results.Int45_MedianPerBlock_group_corr_mean,Group_Results.Int12_MedianPerBlock_group_corr_std,...
%     Group_Results.Int23_MedianPerBlock_group_corr_std,Group_Results.Int34_MedianPerBlock_group_corr_std,Group_Results.Int45_MedianPerBlock_group_corr_std,[titulo '- mediana por bloque - corregido'])
% 
% gcf24=figure(24);
% saveas(gcf24,[path '_Intervals_MedianPerBlock_corr.' 'fig']);
% saveas(gcf24,[path '_Intervals_MedianPerBlock_corr.' 'png']);
% gcf25=figure(25);
% saveas(gcf25,[path '_Intervals_MedianPerBlock_completo_corr.' 'fig']);
% saveas(gcf25,[path '_Intervals_MedianPerBlock_completo_corr.' 'png']);

%media original
Group_Results.Int12_MeanPerBlock_group_mean=nanmean(Group_Parameters.Int12_MeanPerBlock_group,1);
Group_Results.Int12_MeanPerBlock_group_std=nanstd(Group_Parameters.Int12_MeanPerBlock_group)/sqrt(N);
Group_Results.Int23_MeanPerBlock_group_mean=nanmean(Group_Parameters.Int23_MeanPerBlock_group,1);
Group_Results.Int23_MeanPerBlock_group_std=nanstd(Group_Parameters.Int23_MeanPerBlock_group)/sqrt(N);
Group_Results.Int34_MeanPerBlock_group_mean=nanmean(Group_Parameters.Int34_MeanPerBlock_group,1);
Group_Results.Int34_MeanPerBlock_group_std=nanstd(Group_Parameters.Int34_MeanPerBlock_group)/sqrt(N);
Group_Results.Int45_MeanPerBlock_group_mean=nanmean(Group_Parameters.Int45_MeanPerBlock_group,1);
Group_Results.Int45_MeanPerBlock_group_std=nanstd(Group_Parameters.Int45_MeanPerBlock_group)/sqrt(N);

% plot_intervalos(Group_Parameters.Int12_MeanPerBlock_group,Group_Parameters.Int23_MeanPerBlock_group,Group_Parameters.Int34_MeanPerBlock_group,...
%     Group_Parameters.Int45_MeanPerBlock_group,Group_Results.Int12_MeanPerBlock_group_mean,Group_Results.Int23_MeanPerBlock_group_mean,...
%     Group_Results.Int34_MeanPerBlock_group_mean,Group_Results.Int45_MeanPerBlock_group_mean,Group_Results.Int12_MeanPerBlock_group_std,...
%     Group_Results.Int23_MeanPerBlock_group_std,Group_Results.Int34_MeanPerBlock_group_std,Group_Results.Int45_MeanPerBlock_group_std,[titulo '- media por bloque - original'])
% 
% gcf26=figure(26);
% saveas(gcf26,[path '_Intervals_MeanPerBlock.' 'fig']);
% saveas(gcf26,[path '_Intervals_MeanPerBlock.' 'png']);
% gcf27=figure(27);
% saveas(gcf27,[path '_Intervals_MeanPerBlock_completo.' 'fig']);
% saveas(gcf27,[path '_Intervals_MeanPerBlock_completo.' 'png']);

%media corregida
Group_Results.Int12_MeanPerBlock_group_corr_mean=nanmean(Group_Parameters.Int12_MeanPerBlock_corr_group,1);
Group_Results.Int12_MeanPerBlock_group_corr_std=nanstd(Group_Parameters.Int12_MeanPerBlock_corr_group)/sqrt(N);
Group_Results.Int23_MeanPerBlock_group_corr_mean=nanmean(Group_Parameters.Int23_MeanPerBlock_corr_group,1);
Group_Results.Int23_MeanPerBlock_group_corr_std=nanstd(Group_Parameters.Int23_MeanPerBlock_corr_group)/sqrt(N);
Group_Results.Int34_MeanPerBlock_group_corr_mean=nanmean(Group_Parameters.Int34_MeanPerBlock_corr_group,1);
Group_Results.Int34_MeanPerBlock_group_corr_std=nanstd(Group_Parameters.Int34_MeanPerBlock_corr_group)/sqrt(N);
Group_Results.Int45_MeanPerBlock_group_corr_mean=nanmean(Group_Parameters.Int45_MeanPerBlock_corr_group,1);
Group_Results.Int45_MeanPerBlock_group_corr_std=nanstd(Group_Parameters.Int45_MeanPerBlock_corr_group)/sqrt(N);

% plot_intervalos(Group_Parameters.Int12_MeanPerBlock_corr_group,Group_Parameters.Int23_MeanPerBlock_corr_group,Group_Parameters.Int34_MeanPerBlock_corr_group,...
%     Group_Parameters.Int45_MeanPerBlock_corr_group,Group_Results.Int12_MeanPerBlock_group_corr_mean,Group_Results.Int23_MeanPerBlock_group_corr_mean,...
%     Group_Results.Int34_MeanPerBlock_group_corr_mean,Group_Results.Int45_MeanPerBlock_group_corr_mean,Group_Results.Int12_MeanPerBlock_group_corr_std,...
%     Group_Results.Int23_MeanPerBlock_group_corr_std,Group_Results.Int34_MeanPerBlock_group_corr_std,Group_Results.Int45_MeanPerBlock_group_corr_std,[titulo '- media por bloque - corregido'])
% gcf28=figure(28);
% saveas(gcf28,[path '_Intervals_MeanPerBlock_corr.' 'fig']);
% saveas(gcf28,[path '_Intervals_MeanPerBlock_corr.' 'png']);
% gcf29=figure(29);
% saveas(gcf29,[path '_Intervals_MeanPerBlock_completo_corr.' 'fig']);
% saveas(gcf29,[path '_Intervals_MeanPerBlock_completo_corr.' 'png']);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIGURE VISUALIZATION DE CADA TRANSICION POR SEPARDO.                %%%
%%% Uso las mismas funciones que iki visual 
%%%                  20/12/23                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Group_Results.Int12_visual_group_mean=nanmean(Group_Parameters.Int12_visual_group);
Group_Results.Int12_visual_group_std=nanstd(Group_Parameters.Int12_visual_group)/sqrt(N);
mad_12=mad(Group_Parameters.Int12_visual_group,1);
Group_Results.Int23_visual_group_mean=nanmean(Group_Parameters.Int23_visual_group);
Group_Results.Int23_visual_group_std=nanstd(Group_Parameters.Int23_visual_group)/sqrt(N);
mad_23=mad(Group_Parameters.Int23_visual_group,1);
Group_Results.Int34_visual_group_mean=nanmean(Group_Parameters.Int34_visual_group);
Group_Results.Int34_visual_group_std=nanstd(Group_Parameters.Int34_visual_group)/sqrt(N);
mad_34=mad(Group_Parameters.Int34_visual_group,1);
Group_Results.Int45_visual_group_mean=nanmean(Group_Parameters.Int45_visual_group);
Group_Results.Int45_visual_group_std=nanstd(Group_Parameters.Int45_visual_group)/sqrt(N);
mad_45=mad(Group_Parameters.Int45_visual_group,1);
% 
Matriz_media=[Group_Results.Int12_visual_group_mean;Group_Results.Int23_visual_group_mean;Group_Results.Int34_visual_group_mean;Group_Results.Int45_visual_group_mean];
Matriz_std=[Group_Results.Int12_visual_group_std;Group_Results.Int23_visual_group_std;Group_Results.Int34_visual_group_std;Group_Results.Int45_visual_group_std];
% 
% Matriz_media=[Group_Results.Int12_visual_group_mean;NaN(size(Group_Results.Int12_visual_group_mean));NaN(size(Group_Results.Int12_visual_group_mean));Group_Results.Int45_visual_group_mean];
% Matriz_std=[Group_Results.Int12_visual_group_std;NaN(size(Group_Results.Int12_visual_group_mean));NaN(size(Group_Results.Int12_visual_group_mean));Group_Results.Int45_visual_group_std];


% LearningCompare(Matriz_media,Matriz_std,4)
% sgtitle([titulo ' - Learning Curve - Transiciones'])
% saveas(gcf,[path '_Intervalos_Visual_media.' 'fig']);
% saveas(gcf,[path '_Intervalos_Visual_media.' 'png']);
% saveas(gcf,[path '_Intervalos_Visual_media_12_45.' 'fig']);
% saveas(gcf,[path '_Intervalos_Visual_media_12_45.' 'png']);


% PlotInteractivo(Group_Parameters.Int12_visual_group,Group_Results.Int12_visual_group_mean,3*mad_12,'')
% sgtitle([titulo ' - Transicion 12'])
% saveas(gcf,[path '_Plot_completo_Int12.' 'fig']);
% saveas(gcf,[path '_Plot_completo_Int12.' 'png']);
% 
% PlotInteractivo(Group_Parameters.Int23_visual_group,Group_Results.Int23_visual_group_mean,3*mad_23,'')
% sgtitle([titulo ' - Transicion 23'])
% saveas(gcf,[path '_Plot_completo_Int23.' 'fig']);
% saveas(gcf,[path '_Plot_completo_Int23.' 'png']);

% PlotInteractivo(Group_Parameters.Int34_visual_group,Group_Results.Int34_visual_group_mean,3*mad_34,'')
% sgtitle([titulo ' - Transicion 34'])
% saveas(gcf,[path '_Plot_completo_Int34.' 'fig']);
% saveas(gcf,[path '_Plot_completo_Int34.' 'png']);
% 
% PlotInteractivo(Group_Parameters.Int45_visual_group,Group_Results.Int45_visual_group_mean,3*mad_45,'')
% sgtitle([titulo ' - Transicion 45'])
% saveas(gcf,[path '_Plot_completo_Int45.' 'fig']);
% saveas(gcf,[path '_Plot_completo_Int45.' 'png']);

if Group_Parameters.flag_norm==1 || Group_Parameters.flag_filt==1
    Group_Results.Int12_visual_corr_group_mean=nanmean(Group_Parameters.Int12_visual_corr_group);
    Group_Results.Int12_visual_corr_group_std=nanstd(Group_Parameters.Int12_visual_corr_group)/sqrt(N);
    Group_Results.Int23_visual_corr_group_mean=nanmean(Group_Parameters.Int23_visual_corr_group);
    Group_Results.Int23_visual_corr_group_std=nanstd(Group_Parameters.Int23_visual_corr_group)/sqrt(N);
    Group_Results.Int34_visual_corr_group_mean=nanmean(Group_Parameters.Int34_visual_corr_group);
    Group_Results.Int34_visual_corr_group_std=nanstd(Group_Parameters.Int34_visual_corr_group)/sqrt(N);
    Group_Results.Int45_visual_corr_group_mean=nanmean(Group_Parameters.Int45_visual_corr_group);
    Group_Results.Int45_visual_corr_group_std=nanstd(Group_Parameters.Int45_visual_corr_group)/sqrt(N);

    Matriz_media=[Group_Results.Int12_visual_corr_group_mean;Group_Results.Int23_visual_corr_group_mean;Group_Results.Int34_visual_corr_group_mean;Group_Results.Int45_visual_corr_group_mean];
    Matriz_std=[Group_Results.Int12_visual_corr_group_std;Group_Results.Int23_visual_corr_group_std;Group_Results.Int34_visual_corr_group_std;Group_Results.Int45_visual_corr_group_std];

%     LearningCompare(Matriz_media,Matriz_std,4)
%     sgtitle([titulo ' - Learning Curve - Transiciones corregido'])
%     saveas(gcf,[path '_Intervalos_Visual_media_corr.' 'fig']);
%     saveas(gcf,[path '_Intervalos_Visual_media_corr.' 'png']);
% 
%     PlotInteractivo(Group_Parameters.Int12_visual_corr_group,Group_Results.Int12_visual_corr_group_mean,'')
%     sgtitle([titulo ' - Transicion 12_corr'])
%     saveas(gcf,[path '_Plot_completo_Int12_corr.' 'fig']);
%     saveas(gcf,[path '_Plot_completo_Int12_corr.' 'png']);
% 
%     PlotInteractivo(Group_Parameters.Int23_visual_corr_group,Group_Results.Int23_visual_corr_group_mean,'')
%     sgtitle([titulo ' - Transicion 23 corregido'])
%     saveas(gcf,[path '_Plot_completo_Int23_corr.' 'fig']);
%     saveas(gcf,[path '_Plot_completo_Int23_corr.' 'png']);
% 
%     PlotInteractivo(Group_Parameters.Int34_visual_corr_group,Group_Results.Int34_visual_corr_group_mean,'')
%     sgtitle([titulo ' - Transicion 34 corregido'])
%     saveas(gcf,[path '_Plot_completo_Int34_corr.' 'fig']);
%     saveas(gcf,[path '_Plot_completo_Int34_corr.' 'png']);
% 
%     PlotInteractivo(Group_Parameters.Int45_visual_corr_group,Group_Results.Int45_visual_corr_group_mean,'')
%     sgtitle([titulo ' - Transicion 45 corregido'])
%     saveas(gcf,[path '_Plot_completo_Int45_corr.' 'fig']);
%     saveas(gcf,[path '_Plot_completo_Int45_corr.' 'png']);
end

%%
% 
% Group_Results.Int12_visual_group_mean=nanmean(Group_Parameters.Int12_visual_group);
% Group_Results.Int12_visual_group_std=nanstd(Group_Parameters.Int12_visual_group)/sqrt(N);
% mad_12=mad(Group_Parameters.Int12_visual_group,1);
% Group_Results.Int23_visual_group_mean=nanmean(Group_Parameters.Int23_visual_group);
% Group_Results.Int23_visual_group_std=nanstd(Group_Parameters.Int23_visual_group)/sqrt(N);
% mad_23=mad(Group_Parameters.Int23_visual_group,1);
% Group_Results.Int34_visual_group_mean=nanmean(Group_Parameters.Int34_visual_group);
% Group_Results.Int34_visual_group_std=nanstd(Group_Parameters.Int34_visual_group)/sqrt(N);
% mad_34=mad(Group_Parameters.Int34_visual_group,1);
% Group_Results.Int45_visual_group_mean=nanmean(Group_Parameters.Int45_visual_group);
% Group_Results.Int45_visual_group_std=nanstd(Group_Parameters.Int45_visual_group)/sqrt(N);
% mad_45=mad(Group_Parameters.Int45_visual_group,1);
% riz_std=[Group_Results.Int12_visual_group_std;Group_Results.Int23_visual_group_std;Group_Results.Int34_visual_group_std;Group_Results.Int45_visual_group_std];
% % 
% Matriz_media=[Group_Results.Int12_visual_group_mean;NaN(size(Group_Results.Int12_visual_group_mean));NaN(size(Group_Results.Int12_visual_group_mean));Group_Results.Int45_visual_group_mean];
% Matriz_std=[Group_Results.Int12_visual_group_std;NaN(size(Group_Results.Int12_visual_group_mean));NaN(size(Group_Results.Int12_visual_group_mean));Group_Results.Int45_visual_group_std];
% 
% 
% LearningCompare(Matriz_media,Matriz_std,4)
% sgtitle([titulo ' - Learning Curve - Transiciones'])
% saveas(gcf,[path '_Intervalos_Visual_media_12_45.' 'fig']);
% saveas(gcf,[path '_Intervalos_Visual_media_12_45.' 'png']);
% 
%     Group_Results.Int12_visual_corr_group_mean=nanmean(Group_Parameters.Int12_visual_corr_group);
%     Group_Results.Int12_visual_corr_group_std=nanstd(Group_Parameters.Int12_visual_corr_group)/sqrt(N);
%     Group_Results.Int23_visual_corr_group_mean=nanmean(Group_Parameters.Int23_visual_corr_group);
%     Group_Results.Int23_visual_corr_group_std=nanstd(Group_Parameters.Int23_visual_corr_group)/sqrt(N);
%     Group_Results.Int34_visual_corr_group_mean=nanmean(Group_Parameters.Int34_visual_corr_group);
%     Group_Results.Int34_visual_corr_group_std=nanstd(Group_Parameters.Int34_visual_corr_group)/sqrt(N);
%     Group_Results.Int45_visual_corr_group_mean=nanmean(Group_Parameters.Int45_visual_corr_group);
%     Group_Results.Int45_visual_corr_group_std=nanstd(Group_Parameters.Int45_visual_corr_group)/sqrt(N);
% 
% Matriz_media=[Group_Results.Int12_visual_corr_group_mean;NaN(size(Group_Results.Int12_visual_group_mean));NaN(size(Group_Results.Int12_visual_group_mean));Group_Results.Int45_visual_corr_group_mean];
% Matriz_std=[Group_Results.Int12_visual_corr_group_std;NaN(size(Group_Results.Int12_visual_group_mean));NaN(size(Group_Results.Int12_visual_group_mean));Group_Results.Int45_visual_corr_group_std];
% 
% 
% LearningCompare(Matriz_media,Matriz_std,4)
% sgtitle([titulo ' - Learning Curve - Transiciones - corregido'])
% saveas(gcf,[path '_Intervalos_Visual_media_12_45_corr.' 'fig']);
% saveas(gcf,[path '_Intervalos_Visual_media_12_45_corr.' 'png']);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%              FIGURE: visualization micro micro                    %%%
%%%                              7/8/23                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Group_Results.MicroMogs_visual_group_mean = nanmean(Group_Parameters.MicroMOGS_visual_group,1);              
%Group_Results.MicroMogs_visual_group_std = nanstd(Group_Parameters.MicroMOGS_visual_group,0,1)/sqrt(N); 

%Group_Results.MicroMongs_visual_group_mean = nanmean(Group_Parameters.MicroMONGS_visual_group,1);              
%Group_Results.MicroMongs_visual_group_std = nanstd(Group_Parameters.MicroMONGS_visual_group,0,1)/sqrt(N); 

% Visualization_MicroMicroGains_Grupal(Group_Results.MicroMogs_visual_group_mean,Group_Results.MicroMongs_visual_group_mean,...
%   Group_Results.MicroMongs_visual_group_std,Group_Results.MicroMongs_visual_group_std,['Micro Micro Gains - ' titulo ' - Crudo']);
% 
% saveas(gcf,[path '_MicroMicroGains_visual_crudo_Group.' 'fig']);
% saveas(gcf,[path '_MicroMicroGains_visual_crudo_Group.' 'png']);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Micro Micro one point: me quedo con un solo punto por cada micro micro
%%% para poder hacer comparaciones. 
%%% Data original
%%%                  24/12/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %MEDIA por bloque de todos los micro micro 

%media de todas las medias por bloque (mediaTask) - media/mediana across sujeto (mediaS)
% [Group_Results.mediaB_mediaT_mediaS_mogs,Group_Results.mediaB_mediaT_mediaS_mongs,...
%     Group_Results.mediaB_mediaT_medianaS_mogs, Group_Results.mediaB_mediaT_medianaS_mongs]=...
%     MicroMicro_one_point_plot(Group_Parameters.media_MicroMOGS_group,Group_Parameters.media_MicroMONGS_group, 'media',[titulo '- media por bloque - original']);
% 
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_mediaT_crudo_Group.' 'fig']);
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_mediaT_crudo_Group.' 'png']);
% 
% %mediana de todas las medias por bloque (mediaTask) - media/mediana across sujeto (mediaS)
% [Group_Results.mediaB_medianaT_mediaS_mogs,Group_Results.mediaB_medianaT_mediaS_mongs,...
%     Group_Results.mediaB_medianaT_medianaS_mogs, Group_Results.mediaB_medianaT_medianaS_mongs]=...
%     MicroMicro_one_point_plot(Group_Parameters.media_MicroMOGS_group,Group_Parameters.media_MicroMONGS_group, 'mediana', [titulo '- media por bloque - original']);
% 
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_medianaT_crudo_Group.' 'fig']);
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_medianaT_crudo_Group.' 'png']);
% 
% %ultimo punto del acum (acumTask) - media/mediana across sujeto (mediaS)
% [Group_Results.mediaB_acumT_mediaS_mogs,Group_Results.mediaB_acumT_mediaS_mongs,...
%     Group_Results.mediaB_acumT_medianaS_mogs, Group_Results.mediaB_acumT_medianaS_mongs]=...
%     MicroMicro_one_point_plot(Group_Parameters.media_MicroMOGS_acum_group,Group_Parameters.media_MicroMONGS_acum_group, 'acum', [titulo '- media por bloque - original']);
% 
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_acumT_crudo_Group.' 'fig']);
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_acumT_crudo_Group.' 'png']);
% 
%                 %MEDIANA por bloque de todos los micro micro 
% 
% %media de todas las medias por bloque (mediaTask) - media/mediana across sujeto (mediaS)
% [Group_Results.mediaB_mediaT_mediaS_mogs,Group_Results.mediaB_mediaT_mediaS_mongs,...
%     Group_Results.mediaB_mediaT_medianaS_mogs, Group_Results.mediaB_mediaT_medianaS_mongs]=...
%     MicroMicro_one_point_plot(Group_Parameters.media_MicroMOGS_group,Group_Parameters.media_MicroMONGS_group, 'media', [titulo '- mediana por bloque - original']);
% 
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_mediaT_crudo_Group.' 'fig']);
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_mediaT_crudo_Group.' 'png']);
% 
% %mediana de todas las medias por bloque (mediaTask) - media/mediana across sujeto (mediaS)
% [Group_Results.mediaB_medianaT_mediaS_mogs,Group_Results.mediaB_medianaT_mediaS_mongs,...
%     Group_Results.mediaB_medianaT_medianaS_mogs, Group_Results.mediaB_medianaT_medianaS_mongs]=...
%     MicroMicro_one_point_plot(Group_Parameters.media_MicroMOGS_group,Group_Parameters.media_MicroMONGS_group, 'mediana', [titulo '- mediana por bloque - original']);
% 
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_medianaT_crudo_Group.' 'fig']);
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_medianaT_crudo_Group.' 'png']);
% 
% %ultimo punto del acum (acumTask) - media/mediana across sujeto (mediaS)
% [Group_Results.mediaB_acumT_mediaS_mogs,Group_Results.mediaB_acumT_mediaS_mongs,...
%     Group_Results.mediaB_acumT_medianaS_mogs, Group_Results.mediaB_acumT_medianaS_mongs]=...
%     MicroMicro_one_point_plot(Group_Parameters.media_MicroMOGS_acum_group,Group_Parameters.media_MicroMONGS_acum_group, 'acum', [titulo '- mediana por bloque - original']);
% 
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_acumT_crudo_Group.' 'fig']);
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_acumT_crudo_Group.' 'png']);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Micro Micro one point: me quedo con un solo punto por cada micro micro
%%% para poder hacer comparaciones. 
%%% Data CORREGIDA
%%%                  24/12/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %MEDIA por bloque de todos los micro micro 

% %media de todas las medias por bloque (mediaTask) - media/mediana across sujeto (mediaS)
% [Group_Results.mediaB_mediaT_mediaS_mogs,Group_Results.mediaB_mediaT_mediaS_mongs,...
%     Group_Results.mediaB_mediaT_medianaS_mogs, Group_Results.mediaB_mediaT_medianaS_mongs]=...
%     MicroMicro_one_point_plot(Group_Parameters.media_MicroMOGS_corr_group,Group_Parameters.media_MicroMONGS_corr_group, 'media', [titulo '- media por bloque - corr']);
% 
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_mediaT_corr_Group.' 'fig']);
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_mediaT_corr_Group.' 'png']);
% 
% %mediana de todas las medias por bloque (mediaTask) - media/mediana across sujeto (mediaS)
% [Group_Results.mediaB_medianaT_mediaS_mogs,Group_Results.mediaB_medianaT_mediaS_mongs,...
%     Group_Results.mediaB_medianaT_medianaS_mogs, Group_Results.mediaB_medianaT_medianaS_mongs]=...
%     MicroMicro_one_point_plot(Group_Parameters.media_MicroMOGS_corr_group,Group_Parameters.media_MicroMONGS_corr_group, 'mediana',  [titulo '- media por bloque - corr']);
% 
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_medianaT_corr_Group.' 'fig']);
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_medianaT_corr_Group.' 'png']);
% 
% %ultimo punto del acum (acumTask) - media/mediana across sujeto (mediaS)
% [Group_Results.mediaB_acumT_mediaS_mogs,Group_Results.mediaB_acumT_mediaS_mongs,...
%     Group_Results.mediaB_acumT_medianaS_mogs, Group_Results.mediaB_acumT_medianaS_mongs]=...
%     MicroMicro_one_point_plot(Group_Parameters.media_MicroMOGS_corr_acum_group,Group_Parameters.media_MicroMONGS_corr_acum_group, 'acum',  [titulo '- media por bloque - corr']);
% 
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_acumT_corr_Group.' 'fig']);
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_acumT_corr_Group.' 'png']);


%% VIOLINES MICRO MICRO
% 31/12/23  
% for i=1:size(Group_Parameters.media_MicroMOGS_corr_group,1)
%                 %media por bloque (mediaB)
%     %media task
%     Group_Parameters.mediaB_mediaT_mogs(i)=nanmean(Group_Parameters.media_MicroMOGS_corr_group(i,:));
%     Group_Parameters.mediaB_mediaT_mongs(i)=nanmean(Group_Parameters.media_MicroMONGS_corr_group(i,:));
%     %mediana task
%     Group_Parameters.mediaB_medianaT_mogs(i)=nanmedian(Group_Parameters.media_MicroMOGS_corr_group(i,:));
%     Group_Parameters.mediaB_medianaT_mongs(i)=nanmedian(Group_Parameters.media_MicroMONGS_corr_group(i,:));
%     %acum task
%     Group_Parameters.mediaB_acumT_mogs=Group_Parameters.media_MicroMOGS_corr_acum_group(:,end);
%     Group_Parameters.mediaB_acumT_mongs=Group_Parameters.media_MicroMONGS_corr_acum_group(:,end);
%     
%             %mediana por bloque (MedianaB)
%     %media task
%     Group_Parameters.medianaB_mediaT_mogs(i)=nanmean(Group_Parameters.mediana_MicroMOGS_corr_group(i,:));
%     Group_Parameters.medianaB_mediaT_mongs(i)=nanmean(Group_Parameters.mediana_MicroMONGS_corr_group(i,:));
%     %mediana task
%     Group_Parameters.medianaB_medianaT_mogs(i)=nanmedian(Group_Parameters.mediana_MicroMOGS_corr_group(i,:));
%     Group_Parameters.medianaB_medianaT_mongs(i)=nanmedian(Group_Parameters.mediana_MicroMONGS_corr_group(i,:));
%     %acum task
%     Group_Parameters.medianaB_acumT_mogs=Group_Parameters.mediana_MicroMOGS_corr_acum_group(:,end);
%     Group_Parameters.medianaB_acumT_mongs=Group_Parameters.mediana_MicroMONGS_corr_acum_group(:,end);
% end
% 
% figure; set(gcf,'Color','white'); box OFF; hold on; sgtitle([titulo])
% set(gcf,'Position',get(0,'ScreenSize'));
% subplot(2,3,1)
% [h,L,aux_media,aux_mediana,bw]=violin([Group_Parameters.mediaB_mediaT_mogs' Group_Parameters.mediaB_mediaT_mongs'],'xlabel',{'Micro MOGS','Micro MONGS'},'facecolor',[1 0 0; 0 0 1],'edgecolor','k','bw',0.3,'mc','k','ylabel','Mean Value Task (s)','plotlegend',0)
% plot([Group_Parameters.mediaB_mediaT_mogs' Group_Parameters.mediaB_mediaT_mongs']','k.','Markersize',15)
% yline(0)
% ylim([-1.5 1.5])
% title('Media Bloque - Media Task')
% Group_Results.mediaB_mediaT_mediaS_mogs=aux_media(1);
% Group_Results.mediaB_mediaT_mediaS_mongs=aux_media(2);
% Group_Results.mediaB_mediaT_medianaS_mogs=aux_mediana(1);
% Group_Results.mediaB_mediaT_medianaS_mongs=aux_mediana(2);
% clear aux_media;clear aux_mediana;
% 
% 
% subplot(2,3,2)
% [h,L,aux_media,aux_mediana,bw]= violin([Group_Parameters.mediaB_medianaT_mogs' Group_Parameters.mediaB_medianaT_mongs'],'xlabel',{'Micro MOGS','Micro MONGS'},'facecolor',[1 0 0; 0 0 1],'edgecolor','k','bw',0.3,'mc','k','ylabel','Median Value Task (s)','plotlegend',0)
% plot([Group_Parameters.mediaB_mediaT_mogs' Group_Parameters.mediaB_mediaT_mongs']','k.','Markersize',15)
% yline(0)
% ylim([-1.5 1.5])
% title('Media Bloque - Mediana Task')
% Group_Results.mediaB_medianaT_mediaS_mogs= aux_media(1); 
% Group_Results.mediaB_medianaT_mediaS_mongs=aux_media(2);
% Group_Results.mediaB_medianaT_medianaS_mogs=aux_mediana(1);
% Group_Results.mediaB_medianaT_medianaS_mongs=aux_mediana(2);
% clear aux_media;clear aux_mediana;
% 
% subplot(2,3,3)
% [h,L,aux_media,aux_mediana,bw]=violin([Group_Parameters.mediaB_acumT_mogs Group_Parameters.mediaB_acumT_mongs],'xlabel',{'Micro MOGS','Micro MONGS'},'facecolor',[1 0 0; 0 0 1],'edgecolor','k','bw',0.3,'mc','k','ylabel','Last Acum Value Task (s)','plotlegend',0)
% plot([Group_Parameters.mediaB_acumT_mogs Group_Parameters.mediaB_acumT_mongs]','k.','Markersize',15)
% yline(0)
% ylim([-1.5 1.5])
% title('Media Bloque - Acum Task')
% legend('','Mean','Median')
% Group_Results.mediaB_acumT_mediaS_mogs= aux_media(1); 
% Group_Results.mediaB_acumT_mediaS_mongs =aux_media(2);
% Group_Results.mediaB_acumT_medianaS_mogs =aux_mediana(1);
% Group_Results.mediaB_acumT_medianaS_mongs=aux_mediana(2);
% clear aux_media;clear aux_mediana;
% 
% subplot(2,3,4)
% [h,L,aux_media,aux_mediana,bw]=violin([Group_Parameters.medianaB_mediaT_mogs' Group_Parameters.medianaB_mediaT_mongs'],'xlabel',{'Micro MOGS','Micro MONGS'},'facecolor',[1 0 0; 0 0 1],'edgecolor','k','bw',0.3,'mc','k','ylabel','Mean Value Task (s)','plotlegend',0)
% plot([Group_Parameters.medianaB_mediaT_mogs' Group_Parameters.medianaB_mediaT_mongs']','k.','Markersize',15)
% yline(0)
% ylim([-1.5 1.5])
% title('Mediana Bloque - Media Task')
% Group_Results.medianaB_mediaT_mediaS_mogs=aux_media(1);
% Group_Results.medianaB_mediaT_mediaS_mongs=aux_media(2);
% Group_Results.medianaB_mediaT_medianaS_mogs=aux_mediana(1);
% Group_Results.medianaB_mediaT_medianaS_mongs=aux_mediana(2);
% clear aux_media;clear aux_mediana;
% 
% subplot(2,3,5)
% [h,L,aux_media,aux_mediana,bw] =violin([Group_Parameters.medianaB_medianaT_mogs' Group_Parameters.medianaB_medianaT_mongs'],'xlabel',{'Micro MOGS','Micro MONGS'},'facecolor',[1 0 0; 0 0 1],'edgecolor','k','bw',0.3,'mc','k','ylabel','Median Value Task (s)','plotlegend',0)
% plot([Group_Parameters.medianaB_mediaT_mogs' Group_Parameters.medianaB_mediaT_mongs']','k.','Markersize',15)
% yline(0)
% ylim([-1.5 1.5])
% title('Mediana Bloque - Mediana Task')
% Group_Results.medianaB_medianaT_mediaS_mogs=aux_media(1);
% Group_Results.medianaB_medianaT_mediaS_mongs=aux_media(2);
% Group_Results.medianaB_medianaT_medianaS_mogs=aux_mediana(1);
% Group_Results.medianaB_medianaT_medianaS_mongs=aux_mediana(2);
% clear aux_media;clear aux_mediana;
% 
% 
% subplot(2,3,6)
% [h,L,aux_media,aux_mediana,bw]=violin([Group_Parameters.medianaB_acumT_mogs Group_Parameters.medianaB_acumT_mongs],'xlabel',{'Micro MOGS','Micro MONGS'},'facecolor',[1 0 0; 0 0 1],'edgecolor','k','bw',0.3,'mc','k','ylabel','Last Acum Value Task (s)','plotlegend',0)
% plot([Group_Parameters.medianaB_acumT_mogs Group_Parameters.medianaB_acumT_mongs]','k.','Markersize',15)
% yline(0)
% ylim([-1.5 1.5])
% title('Mediana Bloque - Acum Task')
% legend('','Mean','Median')
% Group_Results.medianaB_acumT_mediaS_mogs=aux_media(1);
% Group_Results.medianaB_acumT_mediaS_mongs=aux_media(2);
% Group_Results.medianaB_acumT_medianaS_mogs=aux_mediana(1);
% Group_Results.medianaB_acumT_medianaS_mongs=aux_mediana(2);
% clear aux_media;clear aux_mediana;
% 
% 
% saveas(gcf,[path '_MicroMicroGains_corr_violin.' 'fig']);
% saveas(gcf,[path '_MicroMicroGains_corr_violin.' 'png']);

%%
% 
%                 %MEDIANA por bloque de todos los micro micro 
% 
% %media de todas las medias por bloque (mediaTask) - media/mediana across sujeto (mediaS)
% [Group_Results.mediaB_mediaT_mediaS_mogs,Group_Results.mediaB_mediaT_mediaS_mongs,...
%     Group_Results.mediaB_mediaT_medianaS_mogs, Group_Results.mediaB_mediaT_medianaS_mongs]=...
%     MicroMicro_one_point_plot(Group_Parameters.media_MicroMOGS_corr_group,Group_Parameters.media_MicroMONGS_corr_group, 'media',  [titulo '- mediana por bloque - corr']);
% 
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_mediaT_corr_Group.' 'fig']);
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_mediaT_corr_Group.' 'png']);
% 
% %mediana de todas las medias por bloque (mediaTask) - media/mediana across sujeto (mediaS)
% [Group_Results.mediaB_medianaT_mediaS_mogs,Group_Results.mediaB_medianaT_mediaS_mongs,...
%     Group_Results.mediaB_medianaT_medianaS_mogs, Group_Results.mediaB_medianaT_medianaS_mongs]=...
%     MicroMicro_one_point_plot(Group_Parameters.media_MicroMOGS_corr_group,Group_Parameters.media_MicroMONGS_corr_group, 'mediana',  [titulo '- mediana por bloque - corr']);
% 
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_medianaT_corr_Group.' 'fig']);
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_medianaT_corr_Group.' 'png']);
% 
% %ultimo punto del acum (acumTask) - media/mediana across sujeto (mediaS)
% [Group_Results.mediaB_acumT_mediaS_mogs,Group_Results.mediaB_acumT_mediaS_mongs,...
%     Group_Results.mediaB_acumT_medianaS_mogs, Group_Results.mediaB_acumT_medianaS_mongs]=...
%     MicroMicro_one_point_plot(Group_Parameters.media_MicroMOGS_corr_acum_group,Group_Parameters.media_MicroMONGS_corr_acum_group, 'acum',  [titulo '- mediana por bloque - corr']);
% 
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_acumT_corr_Group.' 'fig']);
% saveas(gcf,[path '_MicroMicroGains_onePoint_mediaB_acumT_corr_Group.' 'png']);

