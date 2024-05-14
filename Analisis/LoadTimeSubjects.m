%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    Guadaupe cutrera 31/8/2023                       %%%
%%% Se cargan los path de cada sujeto para poder darles un análisis de  %%% 
%%% forma grupal. Los sujetos deben estar analizados de forma individual%%%
%%% Es unicamente para los sujetos cuyo bloque termina por tiempo       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Piloto Without ITI - protocolo tiempo - teclado compu
clear;

S_Piloto_tiempo_wo_ITI(1).S='S002\David_S2_Task - sin_ITI_1';
S_Piloto_tiempo_wo_ITI(2).S='S003\claudia_S3_Task - sin_ITI_1';
S_Piloto_tiempo_wo_ITI(3).S='S005\Camila_S5_Task - sin_ITI_1';
S_Piloto_tiempo_wo_ITI(4).S='S006\Lucia_S6_Task - sin_ITI_1'; 
S_Piloto_tiempo_wo_ITI(5).S='S007\Ezequiel_S7_Task - sin_ITI_1';
S_Piloto_tiempo_wo_ITI(6).S='S008\Agatha_S8_Task - sin_ITI_1';
S_Piloto_tiempo_wo_ITI(7).S='S009\Patricio_S9_Task - sin_ITI_1'; 
S_Piloto_tiempo_wo_ITI(8).S='S010\Belen_S10_Task - sin_ITI_1';
S_Piloto_tiempo_wo_ITI(9).S='S011\Cristel_S11_Task - sin_ITI_1';
S_Piloto_tiempo_wo_ITI(10).S='S012\Sabrina_S12_Task - sin_ITI_1';
S_Piloto_tiempo_wo_ITI(11).S='S013\Joaquin_S13_Task - sin_ITI_1';
S_Piloto_tiempo_wo_ITI(12).S='S014\Matias_S14_Task - sin_ITI_1'; %el más rápido.
S_Piloto_tiempo_wo_ITI(13).S='S015\Florencia_S15_Task - sin_ITI_1';
S_Piloto_tiempo_wo_ITI(14).S='S016\Catalina_S16_Task - sin_ITI_1';
S_Piloto_tiempo_wo_ITI(15).S='S017\cristina_s17_Task - sin_ITI_1';
S_Piloto_tiempo_wo_ITI(16).S='S018\Beatriz_S18_Task - sin_ITI_1';

path_Piloto_tiempo_wo_ITI='C:\Users\guada\Desktop\LFA FMED\Codigos Guada\Guada_2023\AnalisisGuada\Behaviour\Piloto MSL time sin ITI\';

%Cargo los seq_results de los sujetos
for i=1:length(S_Piloto_tiempo_wo_ITI)
   Results(i).S=load([path_Piloto_tiempo_wo_ITI S_Piloto_tiempo_wo_ITI(i).S '_results.mat']);
end

%Analisis Grupal 
Resultados_Piloto_tiempo_wo_ITI=ID_RunAnalysisGroup2(Results,[path_Piloto_tiempo_wo_ITI 'Grupal\'],'tiempo','without ITI');


%% Piloto ITI 1s - protocolo tiempo -- teclado
clear;

S_Piloto_tiempo_ITI_1s(1).S='S002\Diego_S2_Task - con_ITI_1';
S_Piloto_tiempo_ITI_1s(2).S='S004\Sofia_S4_Task - con_ITI_1';
S_Piloto_tiempo_ITI_1s(3).S='S005\Gabriela_S5_Task - con_ITI_1';
S_Piloto_tiempo_ITI_1s(4).S='S006\andrea_S6_Task - con_ITI_1';
S_Piloto_tiempo_ITI_1s(5).S='S007\agustina_S7_Task - con_ITI_1';
S_Piloto_tiempo_ITI_1s(6).S='S009\alejandra_S9_Task - con_ITI_1';
S_Piloto_tiempo_ITI_1s(7).S='S011\ana_S11_Task - con_ITI_1';
S_Piloto_tiempo_ITI_1s(8).S='S012\violeta_S12_Task - con_ITI_1';
S_Piloto_tiempo_ITI_1s(9).S='S013\franco bagleone _S13_Task - con_ITI_1';
S_Piloto_tiempo_ITI_1s(10).S='S014\Gonzalo_S14_Task - con_ITI_1';
S_Piloto_tiempo_ITI_1s(11).S='S015\Ana_S15_Task - con_ITI_1';
S_Piloto_tiempo_ITI_1s(12).S='S016\Martina_S16_Task - con_ITI_1';
S_Piloto_tiempo_ITI_1s(13).S='S018\Stefania_S17_Task - con_ITI_1';
S_Piloto_tiempo_ITI_1s(14).S='S019\Sofia_S19_Task - con_ITI_1';
S_Piloto_tiempo_ITI_1s(15).S='S020\Emiliano_S20_Task - con_ITI_1';
S_Piloto_tiempo_ITI_1s(16).S='S021\Aldana_s21_Task - con_ITI_1';

path_Piloto_tiempo_ITI_1s='C:\Users\guada\Desktop\LFA FMED\Codigos Guada\Guada_2023\AnalisisGuada\Behaviour\Piloto MSL time 1s\';

%Cargo los seq_results de los sujetos
for i=1:length(S_Piloto_tiempo_ITI_1s)
   Results(i).S=load([path_Piloto_tiempo_ITI_1s S_Piloto_tiempo_ITI_1s(i).S '_results.mat']);
end

% 
%Analisis Grupal
Resultados_Piloto_tiempo_ITI_1s=ID_RunAnalysisGroup2(Results,[path_Piloto_tiempo_ITI_1s 'Grupal\'],'tiempo','with ITI');

