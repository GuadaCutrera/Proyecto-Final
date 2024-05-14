%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function allows a within subjects analysis. The dependent      %%%
%%% variables of each subject are concatenated in a matrix which is     %%%
%%% the used to calculate mean value and standar error.                 %%%
%%%                                                                     %%%    
%%% Some of the variables require a reshape because they are in a matrix%%%
%%% form and we need them in a vector form.                             %%%
%%%                                                                     %%%
%%% In this scipt is also considered the "time paradigm" in wich the    %%%
%%% subjects not always excute the same amount of sequences in each block%%
%%% In this case, the variables need completion with NaN values, so that%%%
%%% the blocks have the same number of columns and/or rows.             %%%
%%%                                                                     %%%
%%% Parameters:                                                         %%%
%%%   - SUJETOS: array of structs. Each struct has loaded the archive   %%%
%%%   results.mat that comes from Analisis_Seq_Individual 
%%%
%%%   - Path: where to store the figures and results                    %%%
%%%                                                                     %%%
%%%   - paradigm_flag: indicates if the analysis is on time or key paradigm
%%%                                                                     %%%
%%%   - titulo: indicates if the paradigm has  an ITI or not            %%%
%%%                                                                     %%%
%%%                       Guadalupe Cutrera 16/12/2022                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Group_Results=ID_RunAnalysisGroup2(SUJETOS,path,paradigm_flag,titulo)
%% Inicialización de variables

% CorrectSeq_Dur
duration_group=[]; seq_group=[]; 
seq_group_num=[];key_group=[];

%Interkey Interval raw and corrected
iki_group=[]; iki_group_corr=[];
%Interkey Interval visualization raw and corrected
iki_visual_group=[]; iki_visual_group_corr=[];

%Micro Gains Raw
mogs_group=[]; mogs_acumulado_group=[];
mongs_group=[]; mongs_acumulado_group=[];
TL_acumulado_group=[]; TL_group=[];

% Micro Gains corrected
mogs_group_corr=[]; mogs_acumulado_group_corr=[];
mongs_group_corr=[]; mongs_acumulado_group_corr=[];
TL_acumulado_group_corr=[]; TL_group_corr=[];

% Micro Micro Gains Raw
micro_mogs_group=[]; micro_mogs_acumulado_group=[];
micro_mongs_group=[]; micro_mongs_acumulado_group=[];
micro_TL_group=[]; micro_TL_acumulado_group=[];

% Micro Micro Gains corrected
micro_mogs_group_corr=[]; micro_mogs_acumulado_group_corr=[];
micro_mongs_group_corr=[]; micro_mongs_acumulado_group_corr=[];
micro_TL_group_corr=[]; micro_TL_acumulado_group_corr=[];


% Interkey_Matrix - Learning curve per trial
interkey_mat=[]; interkey_mat_corr=[];

%6/1/23
%Micro Gains Mean y median 
media_MicroMOGS_group=[]; media_MicroMONGS_group=[];media_MicroTL_group=[];
media_MicroMOGS_acum_group=[]; media_MicroMONGS_acum_group=[];media_MicroTL_acum_group=[];
mediana_MicroMOGS_group=[]; mediana_MicroMONGS_group=[]; mediana_MicroTL_group=[];
mediana_MicroMOGS_acum_group=[]; mediana_MicroMONGS_acum_group=[]; mediana_MicroTL_acum_group=[];

media_MicroMOGS_corr_group=[]; media_MicroMONGS_corr_group=[]; media_MicroTL_corr_group=[];
media_MicroMOGS_corr_acum_group=[]; media_MicroMONGS_corr_acum_group=[]; media_MicroTL_corr_acum_group=[];
mediana_MicroMOGS_corr_group=[]; mediana_MicroMONGS_corr_group=[];  mediana_MicroTL_corr_group=[];
mediana_MicroMOGS_corr_acum_group=[]; mediana_MicroMONGS_corr_acum_group=[];  mediana_MicroTL_corr_acum_group=[];
% 7/8/23
%MicroMOGS_visual_group=[]; MicroMONGS_visual_group=[]; 

% 15/12/23
%Intervalos
Int12_MeanPerBlock_group=[];Int23_MeanPerBlock_group=[];Int34_MeanPerBlock_group=[];Int45_MeanPerBlock_group=[];
Int12_MeanPerBlock_corr_group=[];Int23_MeanPerBlock_corr_group=[];Int34_MeanPerBlock_corr_group=[];Int45_MeanPerBlock_corr_group=[];

Int12_MedianPerBlock_group=[];Int23_MedianPerBlock_group=[];Int34_MedianPerBlock_group=[];Int45_MedianPerBlock_group=[];
Int12_MedianPerBlock_corr_group=[];Int23_MedianPerBlock_corr_group=[];Int34_MedianPerBlock_corr_group=[];Int45_MedianPerBlock_corr_group=[];

%20/12/23
%intervalos pero con visualizacion de bonstrup
Int12_visual_group=[]; Int23_visual_group=[]; Int34_visual_group=[]; Int45_visual_group=[];
Int12_visual_corr_group=[]; Int23_visual_corr_group=[]; Int34_visual_corr_group=[]; Int45_visual_corr_group=[];

%% visualizacion de los intervalos
% los calculo aca para no tener que pasar a todos los sujetos de nuevo por
% el analisis pero deberia ser parte del analisis individual
for i=1:length(SUJETOS)
    Int12_visual_group(i,:)=Visualization_LearningCurve(SUJETOS(i).S.seq_results.intervalo_12,SUJETOS(i).S.seq_results.seqduration,SUJETOS(i).S.seq_results);
    Int23_visual_group(i,:)=Visualization_LearningCurve(SUJETOS(i).S.seq_results.intervalo_23,SUJETOS(i).S.seq_results.seqduration,SUJETOS(i).S.seq_results);
    Int34_visual_group(i,:)=Visualization_LearningCurve(SUJETOS(i).S.seq_results.intervalo_34,SUJETOS(i).S.seq_results.seqduration,SUJETOS(i).S.seq_results);
    Int45_visual_group(i,:)=Visualization_LearningCurve(SUJETOS(i).S.seq_results.intervalo_45,SUJETOS(i).S.seq_results.seqduration,SUJETOS(i).S.seq_results);
    
    Int12_visual_corr_group(i,:)=Visualization_LearningCurve(SUJETOS(i).S.seq_results.intervalo_12_corr,SUJETOS(i).S.seq_results.seqduration,SUJETOS(i).S.seq_results);
    Int23_visual_corr_group(i,:)=Visualization_LearningCurve(SUJETOS(i).S.seq_results.intervalo_23_corr,SUJETOS(i).S.seq_results.seqduration,SUJETOS(i).S.seq_results);
    Int34_visual_corr_group(i,:)=Visualization_LearningCurve(SUJETOS(i).S.seq_results.intervalo_34_corr,SUJETOS(i).S.seq_results.seqduration,SUJETOS(i).S.seq_results);
    Int45_visual_corr_group(i,:)=Visualization_LearningCurve(SUJETOS(i).S.seq_results.intervalo_45_corr,SUJETOS(i).S.seq_results.seqduration,SUJETOS(i).S.seq_results);
end

%% reorganización de matrices

% PARADIGMA POR TIEMPO
if strcmp(paradigm_flag,'tiempo')==1
    
    %11/6/2023
    %Si el paradigma corta por tiempo debo obtener el sujeto con más cantidad de secuencias realizadas para poder analizar todos los
    %bloques de los sujetos por igual. Esto es porque cada bloque de cada sujeto tiene una cantidad diferente de secuencias. 
    max_seq=0;
    for i=1:length(SUJETOS)
        if max(SUJETOS(i).S.seq_results.correct)>max_seq
            max_seq=max(SUJETOS(i).S.seq_results.correct);
        end
    end
    
    % Sabiendo que el bloque "general" termina en max_seq, completo los bloques con menos secuencias con NaN. Esto permite que todos los
    % sujetos (en todos sus bloques) tengan la "misma cantidad de secuencias". 
    
    for i=1:length(SUJETOS)
        %filas= cantidad de bloques
        %columnas= la diferencia entre el màximo y el actual
        columnas_nan=nan(size(SUJETOS(i).S.seq_results.IKI_per_trial,1),max_seq-size(SUJETOS(i).S.seq_results.IKI_per_trial,2));
        
        % IKI_trial
        SUJETOS(i).S.seq_results.IKI_per_trial=[SUJETOS(i).S.seq_results.IKI_per_trial columnas_nan];
        
        %MICRO MICRO OFFLINE
        SUJETOS(i).S.seq_results.MicroMOGS=[SUJETOS(i).S.seq_results.MicroMOGS columnas_nan];
        SUJETOS(i).S.seq_results.MicroMogs_acum=[SUJETOS(i).S.seq_results.MicroMogs_acum columnas_nan];
        
        
        %MICRO MICRO ONLINE
        SUJETOS(i).S.seq_results.MicroMONGS=[SUJETOS(i).S.seq_results.MicroMONGS columnas_nan]; 
        SUJETOS(i).S.seq_results.MicroMongs_acum=[SUJETOS(i).S.seq_results.MicroMongs_acum columnas_nan];
         
        % MICRO TL
        SUJETOS(i).S.seq_results.MicroTL=[SUJETOS(i).S.seq_results.MicroTL columnas_nan];
        SUJETOS(i).S.seq_results.MicroTL_acum=[SUJETOS(i).S.seq_results.MicroTL_acum columnas_nan];
        
        
        if SUJETOS(i).S.seq_results.flag_norm==1 || SUJETOS(i).S.seq_results.flag_filt==1
            
            SUJETOS(i).S.seq_results.IKI_per_trial_corr=[SUJETOS(i).S.seq_results.IKI_per_trial_corr columnas_nan];
            
            SUJETOS(i).S.seq_results.MicroMOGS_corr=[SUJETOS(i).S.seq_results.MicroMOGS_corr columnas_nan];
           SUJETOS(i).S.seq_results.MicroMogs_corr_acum=[SUJETOS(i).S.seq_results.MicroMogs_corr_acum columnas_nan];
            
            SUJETOS(i).S.seq_results.MicroMONGS_corr=[SUJETOS(i).S.seq_results.MicroMONGS_corr columnas_nan];
           SUJETOS(i).S.seq_results.MicroMongs_corr_acum=[SUJETOS(i).S.seq_results.MicroMongs_corr_acum columnas_nan];  
            
           SUJETOS(i).S.seq_results.MicroTL_corr=[SUJETOS(i).S.seq_results.MicroTL_corr columnas_nan];
           SUJETOS(i).S.seq_results.MicroTL_corr_acum=[SUJETOS(i).S.seq_results.MicroTL_corr_acum columnas_nan];
            
        end
        
    end
   
end

% A los intervalos le calculo la media / mediana por bloque 

for i=1:length(SUJETOS)
    SUJETOS(i).S.seq_results.media_por_bloque_int12=nanmean(SUJETOS(i).S.seq_results.intervalo_12,2);
    SUJETOS(i).S.seq_results.media_por_bloque_int23=nanmean(SUJETOS(i).S.seq_results.intervalo_23,2);
    SUJETOS(i).S.seq_results.media_por_bloque_int34=nanmean(SUJETOS(i).S.seq_results.intervalo_34,2);
    SUJETOS(i).S.seq_results.media_por_bloque_int45=nanmean(SUJETOS(i).S.seq_results.intervalo_45,2);
    
    SUJETOS(i).S.seq_results.mediana_por_bloque_int12=nanmedian(SUJETOS(i).S.seq_results.intervalo_12,2);
    SUJETOS(i).S.seq_results.mediana_por_bloque_int23=nanmedian(SUJETOS(i).S.seq_results.intervalo_23,2);
    SUJETOS(i).S.seq_results.mediana_por_bloque_int34=nanmedian(SUJETOS(i).S.seq_results.intervalo_34,2);
    SUJETOS(i).S.seq_results.mediana_por_bloque_int45=nanmedian(SUJETOS(i).S.seq_results.intervalo_45,2);
    
     if SUJETOS(i).S.seq_results.flag_norm==1 || SUJETOS(i).S.seq_results.flag_filt==1
        
         SUJETOS(i).S.seq_results.media_por_bloque_int12_corr=nanmean(SUJETOS(i).S.seq_results.intervalo_12_corr,2);
         SUJETOS(i).S.seq_results.media_por_bloque_int23_corr=nanmean(SUJETOS(i).S.seq_results.intervalo_23_corr,2);
         SUJETOS(i).S.seq_results.media_por_bloque_int34_corr=nanmean(SUJETOS(i).S.seq_results.intervalo_34_corr,2);
         SUJETOS(i).S.seq_results.media_por_bloque_int45_corr=nanmean(SUJETOS(i).S.seq_results.intervalo_45_corr,2);
         
         SUJETOS(i).S.seq_results.mediana_por_bloque_int12_corr=nanmedian(SUJETOS(i).S.seq_results.intervalo_12_corr,2);
         SUJETOS(i).S.seq_results.mediana_por_bloque_int23_corr=nanmedian(SUJETOS(i).S.seq_results.intervalo_23_corr,2);
         SUJETOS(i).S.seq_results.mediana_por_bloque_int34_corr=nanmedian(SUJETOS(i).S.seq_results.intervalo_34_corr,2);
         SUJETOS(i).S.seq_results.mediana_por_bloque_int45_corr=nanmedian(SUJETOS(i).S.seq_results.intervalo_45_corr,2);
     end
end

% Matriz IKI_per_trial y a los micro micro - la pongo en una fila para poder graficarla. Los vector es de noBlock x noSeq
for i=1:length(SUJETOS)
    SUJETOS(i).S.seq_results.IKI_per_trial=reshape(SUJETOS(i).S.seq_results.IKI_per_trial',1,[]);
    
    SUJETOS(i).S.seq_results.MicroMOGS=reshape(SUJETOS(i).S.seq_results.MicroMOGS',1,[]);
    SUJETOS(i).S.seq_results.MicroMONGS=reshape(SUJETOS(i).S.seq_results.MicroMONGS',1,[]);
    SUJETOS(i).S.seq_results.MicroTL=reshape(SUJETOS(i).S.seq_results.MicroTL',1,[]);
    
    SUJETOS(i).S.seq_results.MicroMogs_acum=reshape(SUJETOS(i).S.seq_results.MicroMogs_acum',1,[]);
    SUJETOS(i).S.seq_results.MicroMongs_acum=reshape(SUJETOS(i).S.seq_results.MicroMongs_acum',1,[]);
    SUJETOS(i).S.seq_results.MicroTL_acum=reshape(SUJETOS(i).S.seq_results.MicroTL_acum',1,[]);
    
    SUJETOS(i).S.seq_results.media_por_bloque_int12=SUJETOS(i).S.seq_results.media_por_bloque_int12';
    SUJETOS(i).S.seq_results.media_por_bloque_int23=SUJETOS(i).S.seq_results.media_por_bloque_int23';
    SUJETOS(i).S.seq_results.media_por_bloque_int34=SUJETOS(i).S.seq_results.media_por_bloque_int34';
    SUJETOS(i).S.seq_results.media_por_bloque_int45=SUJETOS(i).S.seq_results.media_por_bloque_int45';
    
    SUJETOS(i).S.seq_results.mediana_por_bloque_int12=SUJETOS(i).S.seq_results.mediana_por_bloque_int12';
    SUJETOS(i).S.seq_results.mediana_por_bloque_int23=SUJETOS(i).S.seq_results.mediana_por_bloque_int23';
    SUJETOS(i).S.seq_results.mediana_por_bloque_int34=SUJETOS(i).S.seq_results.mediana_por_bloque_int34';
    SUJETOS(i).S.seq_results.mediana_por_bloque_int45=SUJETOS(i).S.seq_results.mediana_por_bloque_int45';

    
    %11/6/2203
    if SUJETOS(i).S.seq_results.flag_norm==1 || SUJETOS(i).S.seq_results.flag_filt==1
       
        SUJETOS(i).S.seq_results.IKI_per_trial_corr=reshape(SUJETOS(i).S.seq_results.IKI_per_trial_corr',1,[]);
        
        SUJETOS(i).S.seq_results.MicroMOGS_corr=reshape(SUJETOS(i).S.seq_results.MicroMOGS_corr',1,[]);
        SUJETOS(i).S.seq_results.MicroMONGS_corr=reshape(SUJETOS(i).S.seq_results.MicroMONGS_corr',1,[]);
        SUJETOS(i).S.seq_results.MicroTL_corr=reshape(SUJETOS(i).S.seq_results.MicroTL_corr',1,[]);

        SUJETOS(i).S.seq_results.MicroMogs_corr_acum=reshape(SUJETOS(i).S.seq_results.MicroMogs_corr_acum',1,[]);
        SUJETOS(i).S.seq_results.MicroMongs_corr_acum=reshape(SUJETOS(i).S.seq_results.MicroMongs_corr_acum',1,[]);
        SUJETOS(i).S.seq_results.MicroTL_corr_acum=reshape(SUJETOS(i).S.seq_results.MicroTL_corr_acum',1,[]);
        
        SUJETOS(i).S.seq_results.media_por_bloque_int12_corr=SUJETOS(i).S.seq_results.media_por_bloque_int12_corr';
        SUJETOS(i).S.seq_results.media_por_bloque_int23_corr=SUJETOS(i).S.seq_results.media_por_bloque_int23_corr';
        SUJETOS(i).S.seq_results.media_por_bloque_int34_corr=SUJETOS(i).S.seq_results.media_por_bloque_int34_corr';
        SUJETOS(i).S.seq_results.media_por_bloque_int45_corr=SUJETOS(i).S.seq_results.media_por_bloque_int45_corr';

        SUJETOS(i).S.seq_results.mediana_por_bloque_int12_corr=SUJETOS(i).S.seq_results.mediana_por_bloque_int12_corr';
        SUJETOS(i).S.seq_results.mediana_por_bloque_int23_corr=SUJETOS(i).S.seq_results.mediana_por_bloque_int23_corr';
        SUJETOS(i).S.seq_results.mediana_por_bloque_int34_corr=SUJETOS(i).S.seq_results.mediana_por_bloque_int34_corr';
        SUJETOS(i).S.seq_results.mediana_por_bloque_int45_corr=SUJETOS(i).S.seq_results.mediana_por_bloque_int45_corr';
    end

    
end

%% todos los sujetos en forma matricial
 for i=1:length(SUJETOS)
    duration_group=[duration_group;SUJETOS(i).S.seq_results.GOduration];
    seq_group=[seq_group;SUJETOS(i).S.seq_results.SEQduration];
    seq_group_num=[seq_group_num;SUJETOS(i).S.seq_results.correct];
    key_group=[key_group;SUJETOS(i).S.seq_results.Intervalmean];
    
    iki_group=[iki_group;SUJETOS(i).S.seq_results.IKI_per_trial];
    iki_visual_group=[iki_visual_group;SUJETOS(i).S.seq_results.IKI_visual];
    
    mogs_group=[mogs_group;SUJETOS(i).S.seq_results.MOGS'];
    mogs_acumulado_group=[mogs_acumulado_group;SUJETOS(i).S.seq_results.MOGS_acumulativo'];
    mongs_group=[mongs_group;SUJETOS(i).S.seq_results.MONGS];
    mongs_acumulado_group=[mongs_acumulado_group;SUJETOS(i).S.seq_results.MONGS_acumulativo];
    TL_group=[TL_group;SUJETOS(i).S.seq_results.Total_Learning];
    TL_acumulado_group=[TL_acumulado_group;SUJETOS(i).S.seq_results.Total_Learning_acumulativo];
    
    
    micro_mogs_group=[micro_mogs_group;SUJETOS(i).S.seq_results.MicroMOGS];
    micro_mongs_group=[micro_mongs_group;SUJETOS(i).S.seq_results.MicroMONGS];
    micro_TL_group=[micro_TL_group;SUJETOS(i).S.seq_results.MicroTL];
    micro_mogs_acumulado_group=[micro_mogs_acumulado_group;SUJETOS(i).S.seq_results.MicroMogs_acum];
    micro_mongs_acumulado_group=[micro_mongs_acumulado_group;SUJETOS(i).S.seq_results.MicroMongs_acum];
    micro_TL_acumulado_group=[micro_TL_acumulado_group;SUJETOS(i).S.seq_results.MicroTL_acum];

%     interkey_mat=[interkey_mat;SUJETOS(i).S.seq_results.interkey_matrix];
% 
%     if SUJETOS(i).S.seq_results.flag_filt==1
%         interkey_mat_corr=[interkey_mat_corr;SUJETOS(i).S.seq_results.interkey_matrix_corr];
%     end
    
    %6/1/23
    media_MicroMOGS_group=[media_MicroMOGS_group;SUJETOS(i).S.seq_results.media_MicroMogs];
    media_MicroMONGS_group=[media_MicroMONGS_group;SUJETOS(i).S.seq_results.media_MicroMongs];
    media_MicroTL_group=[media_MicroTL_group;SUJETOS(i).S.seq_results.media_MicroTL];
    media_MicroMOGS_acum_group=[media_MicroMOGS_acum_group;SUJETOS(i).S.seq_results.media_MicroMogs_acum];
    media_MicroMONGS_acum_group=[media_MicroMONGS_acum_group;SUJETOS(i).S.seq_results.media_MicroMongs_acum];
    media_MicroTL_acum_group=[media_MicroTL_acum_group;SUJETOS(i).S.seq_results.media_MicroTL_acum];
    
    mediana_MicroMOGS_group=[mediana_MicroMOGS_group;SUJETOS(i).S.seq_results.mediana_MicroMogs];
    mediana_MicroMONGS_group=[mediana_MicroMONGS_group;SUJETOS(i).S.seq_results.mediana_MicroMongs];
    mediana_MicroTL_group=[mediana_MicroTL_group;SUJETOS(i).S.seq_results.mediana_MicroTL];
    mediana_MicroMOGS_acum_group=[mediana_MicroMOGS_acum_group;SUJETOS(i).S.seq_results.mediana_MicroMogs_acum];
    mediana_MicroMONGS_acum_group=[mediana_MicroMONGS_acum_group;SUJETOS(i).S.seq_results.mediana_MicroMongs_acum];
    mediana_MicroTL_acum_group=[mediana_MicroTL_acum_group;SUJETOS(i).S.seq_results.mediana_MicroTL_acum];
 
    % data corregida
    media_MicroMOGS_corr_group=[media_MicroMOGS_corr_group;SUJETOS(i).S.seq_results.media_MicroMogs_corr];
    media_MicroMONGS_corr_group=[media_MicroMONGS_corr_group;SUJETOS(i).S.seq_results.media_MicroMongs_corr];
    media_MicroTL_corr_group=[media_MicroTL_corr_group;SUJETOS(i).S.seq_results.media_MicroTL_corr];
    media_MicroMOGS_corr_acum_group=[media_MicroMOGS_corr_acum_group;SUJETOS(i).S.seq_results.media_MicroMogs_corr_acum];
    media_MicroMONGS_corr_acum_group=[media_MicroMONGS_corr_acum_group;SUJETOS(i).S.seq_results.media_MicroMongs_corr_acum];
    media_MicroTL_corr_acum_group=[media_MicroTL_corr_acum_group;SUJETOS(i).S.seq_results.media_MicroTL_corr_acum];
    
    mediana_MicroMOGS_corr_group=[mediana_MicroMOGS_corr_group;SUJETOS(i).S.seq_results.mediana_MicroMogs_corr];
    mediana_MicroMONGS_corr_group=[mediana_MicroMONGS_corr_group;SUJETOS(i).S.seq_results.mediana_MicroMongs_corr];
    mediana_MicroTL_corr_group=[mediana_MicroTL_corr_group;SUJETOS(i).S.seq_results.mediana_MicroTL_corr];
    mediana_MicroMOGS_corr_acum_group=[mediana_MicroMOGS_corr_acum_group;SUJETOS(i).S.seq_results.mediana_MicroMogs_corr_acum];
    mediana_MicroMONGS_corr_acum_group=[mediana_MicroMONGS_corr_acum_group;SUJETOS(i).S.seq_results.mediana_MicroMongs_corr_acum];
    mediana_MicroTL_corr_acum_group=[mediana_MicroTL_corr_acum_group;SUJETOS(i).S.seq_results.mediana_MicroTL_corr_acum];

    % 7/8/23
%      MicroMOGS_visual_group=[MicroMOGS_visual_group;SUJETOS(i).S.seq_results.MicroMOGS_visual];
%      MicroMONGS_visual_group=[MicroMONGS_visual_group;SUJETOS(i).S.seq_results.MicroMONGS_visual];
    
    %INTERVALOS  15/12/23
    Int12_MeanPerBlock_group=[Int12_MeanPerBlock_group;SUJETOS(i).S.seq_results.media_por_bloque_int12];
    Int23_MeanPerBlock_group=[Int23_MeanPerBlock_group;SUJETOS(i).S.seq_results.media_por_bloque_int23];
    Int34_MeanPerBlock_group=[Int34_MeanPerBlock_group;SUJETOS(i).S.seq_results.media_por_bloque_int34];
    Int45_MeanPerBlock_group=[Int45_MeanPerBlock_group;SUJETOS(i).S.seq_results.media_por_bloque_int45];

    Int12_MedianPerBlock_group=[Int12_MedianPerBlock_group;SUJETOS(i).S.seq_results.mediana_por_bloque_int12];
    Int23_MedianPerBlock_group=[Int23_MedianPerBlock_group;SUJETOS(i).S.seq_results.mediana_por_bloque_int23];
    Int34_MedianPerBlock_group=[Int34_MedianPerBlock_group;SUJETOS(i).S.seq_results.mediana_por_bloque_int34];
    Int45_MedianPerBlock_group=[Int45_MedianPerBlock_group;SUJETOS(i).S.seq_results.mediana_por_bloque_int45];
    
    
    if SUJETOS(i).S.seq_results.flag_norm==1 || SUJETOS(i).S.seq_results.flag_filt==1
        % 11/6/2023
        iki_group_corr=[iki_group_corr;SUJETOS(i).S.seq_results.IKI_per_trial_corr];
        iki_visual_group_corr=[iki_visual_group_corr;SUJETOS(i).S.seq_results.IKI_visual_corr];
            
        mogs_group_corr=[mogs_group_corr;SUJETOS(i).S.seq_results.MOGS_corr'];
        mogs_acumulado_group_corr=[mogs_acumulado_group_corr;SUJETOS(i).S.seq_results.MOGS_acumulativo_corr'];
        mongs_group_corr=[mongs_group_corr;SUJETOS(i).S.seq_results.MONGS_corr];
        mongs_acumulado_group_corr=[mongs_acumulado_group_corr;SUJETOS(i).S.seq_results.MONGS_acumulativo_corr];
        TL_group_corr=[TL_group_corr;SUJETOS(i).S.seq_results.Total_Learning_corr];
        TL_acumulado_group_corr=[TL_acumulado_group_corr;SUJETOS(i).S.seq_results.Total_Learning_acumulativo_corr];
    
        micro_mogs_group_corr=[micro_mogs_group_corr;SUJETOS(i).S.seq_results.MicroMOGS_corr];
        micro_mongs_group_corr=[micro_mongs_group_corr;SUJETOS(i).S.seq_results.MicroMONGS_corr];
        micro_TL_group_corr=[micro_TL_group_corr;SUJETOS(i).S.seq_results.MicroTL_corr];
        micro_mogs_acumulado_group_corr=[micro_mogs_acumulado_group_corr;SUJETOS(i).S.seq_results.MicroMogs_corr_acum];
        micro_mongs_acumulado_group_corr=[micro_mongs_acumulado_group_corr;SUJETOS(i).S.seq_results.MicroMongs_corr_acum];
        micro_TL_acumulado_group_corr=[micro_TL_acumulado_group_corr;SUJETOS(i).S.seq_results.MicroTL_corr_acum];
        
        %INTERVALOS  15/12/23
        Int12_MeanPerBlock_corr_group=[Int12_MeanPerBlock_corr_group;SUJETOS(i).S.seq_results.media_por_bloque_int12_corr];
        Int23_MeanPerBlock_corr_group=[Int23_MeanPerBlock_corr_group;SUJETOS(i).S.seq_results.media_por_bloque_int23_corr];
        Int34_MeanPerBlock_corr_group=[Int34_MeanPerBlock_corr_group;SUJETOS(i).S.seq_results.media_por_bloque_int34_corr];
        Int45_MeanPerBlock_corr_group=[Int45_MeanPerBlock_corr_group;SUJETOS(i).S.seq_results.media_por_bloque_int45_corr];
        
        Int12_MedianPerBlock_corr_group=[Int12_MedianPerBlock_corr_group;SUJETOS(i).S.seq_results.mediana_por_bloque_int12_corr];
        Int23_MedianPerBlock_corr_group=[Int23_MedianPerBlock_corr_group;SUJETOS(i).S.seq_results.mediana_por_bloque_int23_corr];
        Int34_MedianPerBlock_corr_group=[Int34_MedianPerBlock_corr_group;SUJETOS(i).S.seq_results.mediana_por_bloque_int34_corr];
        Int45_MedianPerBlock_corr_group=[Int45_MedianPerBlock_corr_group;SUJETOS(i).S.seq_results.mediana_por_bloque_int45_corr];
    end


 end

 %% Pongo los parámetros en un struct
Group_Parameters.duration_group=duration_group;
Group_Parameters.seq_group=seq_group;
Group_Parameters.seq_group_num=seq_group_num;
Group_Parameters.key_group=key_group;

Group_Parameters.iki_group=iki_group;
Group_Parameters.iki_visual_group=iki_visual_group;

mogs_acumulado_group(mogs_acumulado_group==0)=NaN;
mongs_acumulado_group(mongs_acumulado_group==0)=NaN;
TL_acumulado_group(TL_acumulado_group==0)=NaN;

Group_Parameters.mogs_group=mogs_group;
Group_Parameters.mogs_acumulado_group=mogs_acumulado_group;
Group_Parameters.mongs_group=mongs_group;
Group_Parameters.mongs_acumulado_group=mongs_acumulado_group;
Group_Parameters.TL_group=TL_group;
Group_Parameters.TL_acumulado_group=TL_acumulado_group;


Group_Parameters.micro_mogs_group=micro_mogs_group;
Group_Parameters.micro_mogs_acumulado_group=micro_mogs_acumulado_group;
Group_Parameters.micro_mongs_group=micro_mongs_group;
Group_Parameters.micro_mongs_acumulado_group=micro_mongs_acumulado_group;
Group_Parameters.micro_TL_group=micro_TL_group;
Group_Parameters.micro_TL_acumulado_group=micro_TL_acumulado_group;

% Group_Parameters.MicroMOGS_visual_group=MicroMOGS_visual_group;
% Group_Parameters.MicroMONGS_visual_group=MicroMONGS_visual_group;

%16/12/23
Group_Parameters.Int12_MeanPerBlock_group=Int12_MeanPerBlock_group;
Group_Parameters.Int23_MeanPerBlock_group=Int23_MeanPerBlock_group;
Group_Parameters.Int34_MeanPerBlock_group=Int34_MeanPerBlock_group;
Group_Parameters.Int45_MeanPerBlock_group=Int45_MeanPerBlock_group;

Group_Parameters.Int12_MedianPerBlock_group=Int12_MedianPerBlock_group;
Group_Parameters.Int23_MedianPerBlock_group=Int23_MedianPerBlock_group;
Group_Parameters.Int34_MedianPerBlock_group=Int34_MedianPerBlock_group;
Group_Parameters.Int45_MedianPerBlock_group=Int45_MedianPerBlock_group;

Group_Parameters.Int12_visual_group=Int12_visual_group;
Group_Parameters.Int23_visual_group=Int23_visual_group;
Group_Parameters.Int34_visual_group=Int34_visual_group;
Group_Parameters.Int45_visual_group=Int45_visual_group;


% 11/6/2023
if SUJETOS(i).S.seq_results.flag_norm==1 || SUJETOS(i).S.seq_results.flag_filt==1
    
    Group_Parameters.iki_group_corr=iki_group_corr;
    Group_Parameters.iki_visual_group_corr=iki_visual_group_corr;

    mogs_acumulado_group_corr(mogs_acumulado_group_corr==0)=NaN;
    mongs_acumulado_group_corr(mongs_acumulado_group_corr==0)=NaN;
    TL_acumulado_group_corr(TL_acumulado_group_corr==0)=NaN;
    
    Group_Parameters.mogs_group_corr=mogs_group_corr;
    Group_Parameters.mogs_acumulado_group_corr=mogs_acumulado_group_corr;
    Group_Parameters.mongs_group_corr=mongs_group_corr;
    Group_Parameters.mongs_acumulado_group_corr=mongs_acumulado_group_corr;
    Group_Parameters.TL_group_corr=TL_group_corr;
    Group_Parameters.TL_acumulado_group_corr=TL_acumulado_group_corr;

    Group_Parameters.micro_mogs_group_corr=micro_mogs_group_corr;
    Group_Parameters.micro_mogs_acumulado_group_corr=micro_mogs_acumulado_group_corr;
    Group_Parameters.micro_mongs_group_corr=micro_mongs_group_corr;
    Group_Parameters.micro_mongs_acumulado_group_corr=micro_mongs_acumulado_group_corr;
    Group_Parameters.micro_TL_group_corr=micro_TL_group_corr;
    Group_Parameters.micro_TL_acumulado_group_corr=micro_TL_acumulado_group_corr;

    % 7/12/23
    Group_Parameters.media_MicroMOGS_corr_group=media_MicroMOGS_corr_group;
    Group_Parameters.media_MicroMONGS_corr_group=media_MicroMONGS_corr_group;
    Group_Parameters.media_MicroTL_corr_group=media_MicroTL_corr_group;

    media_MicroMOGS_corr_acum_group(media_MicroMOGS_corr_acum_group==0)=NaN;
    media_MicroMONGS_corr_acum_group(media_MicroMONGS_corr_acum_group==0)=NaN;
    media_MicroTL_corr_acum_group(media_MicroTL_corr_acum_group==0)=NaN;
    Group_Parameters.media_MicroMOGS_corr_acum_group=media_MicroMOGS_corr_acum_group;
    Group_Parameters.media_MicroMONGS_corr_acum_group=media_MicroMONGS_corr_acum_group;
    Group_Parameters.media_MicroTL_corr_acum_group=media_MicroTL_corr_acum_group;

    Group_Parameters.mediana_MicroMOGS_corr_group=mediana_MicroMOGS_corr_group;
    Group_Parameters.mediana_MicroMONGS_corr_group=mediana_MicroMONGS_corr_group;
    Group_Parameters.mediana_MicroTL_corr_group=mediana_MicroTL_corr_group;

    mediana_MicroMOGS_corr_acum_group(mediana_MicroMOGS_corr_acum_group==0)=NaN;
    mediana_MicroMONGS_corr_acum_group(mediana_MicroMONGS_corr_acum_group==0)=NaN;
    mediana_MicroTL_corr_acum_group(mediana_MicroTL_corr_acum_group==0)=NaN;
    Group_Parameters.mediana_MicroMOGS_corr_acum_group=mediana_MicroMOGS_corr_acum_group;
    Group_Parameters.mediana_MicroMONGS_corr_acum_group=mediana_MicroMONGS_corr_acum_group;
    Group_Parameters.mediana_MicroTL_corr_acum_group=mediana_MicroTL_corr_acum_group;
    
    %16/12/23
    Group_Parameters.Int12_MeanPerBlock_corr_group=Int12_MeanPerBlock_corr_group;
    Group_Parameters.Int23_MeanPerBlock_corr_group=Int23_MeanPerBlock_corr_group;
    Group_Parameters.Int34_MeanPerBlock_corr_group=Int34_MeanPerBlock_corr_group;
    Group_Parameters.Int45_MeanPerBlock_corr_group=Int45_MeanPerBlock_corr_group;
    
    Group_Parameters.Int12_MedianPerBlock_corr_group=Int12_MedianPerBlock_corr_group;
    Group_Parameters.Int23_MedianPerBlock_corr_group=Int23_MedianPerBlock_corr_group;
    Group_Parameters.Int34_MedianPerBlock_corr_group=Int34_MedianPerBlock_corr_group;
    Group_Parameters.Int45_MedianPerBlock_corr_group=Int45_MedianPerBlock_corr_group;
    
    Group_Parameters.Int12_visual_corr_group=Int12_visual_corr_group;
    Group_Parameters.Int23_visual_corr_group=Int23_visual_corr_group;
    Group_Parameters.Int34_visual_corr_group=Int34_visual_corr_group;
    Group_Parameters.Int45_visual_corr_group=Int45_visual_corr_group;
end

 
Group_Parameters.interkey_mat=interkey_mat;
if SUJETOS(1).S.seq_results.flag_filt==1
    Group_Parameters.interkey_mat_corr=interkey_mat_corr;
end

% 6/1/23
Group_Parameters.media_MicroMOGS_group=media_MicroMOGS_group;
Group_Parameters.media_MicroMONGS_group=media_MicroMONGS_group;
Group_Parameters.media_MicroTL_group=media_MicroTL_group;

media_MicroMOGS_acum_group(media_MicroMOGS_acum_group==0)=NaN;
media_MicroMONGS_acum_group(media_MicroMONGS_acum_group==0)=NaN;
media_MicroTL_acum_group(media_MicroTL_acum_group==0)=NaN;
Group_Parameters.media_MicroMOGS_acum_group=media_MicroMOGS_acum_group;
Group_Parameters.media_MicroMONGS_acum_group=media_MicroMONGS_acum_group;
Group_Parameters.media_MicroTL_acum_group=media_MicroTL_acum_group;

Group_Parameters.mediana_MicroMOGS_group=mediana_MicroMOGS_group;
Group_Parameters.mediana_MicroMONGS_group=mediana_MicroMONGS_group;
Group_Parameters.mediana_MicroTL_group=mediana_MicroTL_group;

mediana_MicroMOGS_acum_group(mediana_MicroMOGS_acum_group==0)=NaN;
mediana_MicroMONGS_acum_group(mediana_MicroMONGS_acum_group==0)=NaN;
mediana_MicroTL_acum_group(mediana_MicroTL_acum_group==0)=NaN;
Group_Parameters.mediana_MicroMOGS_acum_group=mediana_MicroMOGS_acum_group;
Group_Parameters.mediana_MicroMONGS_acum_group=mediana_MicroMONGS_acum_group;
Group_Parameters.mediana_MicroTL_acum_group=mediana_MicroTL_acum_group;

% 18/1/23
if SUJETOS(1).S.seq_results.flag_norm==0 %grupo no normalizado
    Group_Parameters.flag_norm=0;
    Group_Parameters.titulo_norm='sin norm';
else %grupo normalizado
    Group_Parameters.flag_norm=1;
    Group_Parameters.flag_tipo_norm=SUJETOS(1).S.seq_results.flag_tipo_norm;
    Group_Parameters.titulo_norm=SUJETOS(1).S.seq_results.flag_tipo_norm;
end
% 20/5/2023
if SUJETOS(1).S.seq_results.flag_filt==0 %grupo no filtrado
    Group_Parameters.flag_filt=0;
    Group_Parameters.titulo_filt='sin filt';
else %grupo filtrado
    Group_Parameters.flag_filt=1;
    Group_Parameters.titulo_filt='filt';
end
Group_Parameters.titulo_analisis=[Group_Parameters.titulo_filt ' - ' Group_Parameters.titulo_norm];

if strcmp(paradigm_flag,'tiempo')==1
    %11/6/2023
    Group_Parameters.cant_SeqBlock=max_seq;
else
    Group_Parameters.cant_SeqBlock=SUJETOS(1).S.seq_results.cant_SeqBlock;
end


%% clear workspace
clear duration_group; clear seq_group; clear key_group; clear seq_group_num;
clear iki_group; clear iki_group_corr;

clear mogs_group; clear mogs_acumulado_group; clear mongs_acumulado_group; clear mongs_group;
clear TL_group; clear TL_acumulado_group;

clear mogs_group_corr; clear mogs_acumulado_group_corr; clear mongs_acumulado_group_corr; clear mongs_group_corr;
clear TL_group_corr; clear TL_acumulado_group_corr;

clear micro_mogs_group; clear micro_mogs_acumulado_group; clear micro_mongs_acumulado_group; clear micro_mongs_group;

%clear micro_mogs_group_corr; clear micro_mogs_acumulado_group_corr; clear micro_mongs_acumulado_group_corr; clear micro_mongs_group_corr;

% 6/1/23
clear media_MicroMOGS_group; clear media_MicroMONGS_group; clear mediana_MicroMOGS_acum_group; clear mediana_MicroMONGS_acum_group; 


%% Plot figures
Group_Results = Plot_Grupal_Figures(Group_Parameters,path,titulo,paradigm_flag,length(SUJETOS));

addpath('C:\Users\physi\Documents\Guada_2022\MSL guada\Task_MSL\stim-master\experiments\Funciones para analisis MSL\')

Group_Results.Group_Parameters = Group_Parameters;
%% Guardo los datos en resultados

if strcmp(paradigm_flag,'tiempo')==1
    Group_Results.max_seq=max_seq;
end

% 18/1/23
Group_Results.flag_norm=Group_Parameters.flag_norm;
Group_Results.flag_filt=Group_Parameters.flag_filt;

save([path 'Group_Results.mat'],'Group_Results');