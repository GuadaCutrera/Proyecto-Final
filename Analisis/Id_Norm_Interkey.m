%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         Id_Norm_Interkey                            %%%
%%% Esta funcion permite normalizar la variable IKI_trial (Interkey     %%%
%%% Interval per trial). Con esto queremos evitar inducir el fenómeno de %%
%%% anticorrelación                                                     %%%
%%% entre MOGS y MONGS. Además, para el análisis grupal, permite que haya%%
%%% menos variablidad entre sujetos.                                    %%%
%%% Por otro lado, se normalizan a los intervalos de forma individual y %%%   
%%% por bloque con un zscore, 01 o un centrado.                         %%%
%%%                                                                     %%%
%%% Parámetros:                                                         %%%
%%%   -  IKI_trial: este ya se encuentra filtrado (bloque a bloque) por el% 
%%%   método seleccionado en Id_Filt_Interkey pero no está normalizado  %%%
%%%   - interval_ij: pueden estar filtrados o ser data cruda.           %%%
%%%                                                                     %%%
%%%                      Guadalupe Cutrera 13/1/23                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [interval12,interval23,interval34,interval45,IKI_trial_norm] = ...
    Id_Norm_Interkey(interval12,interval23,interval34,interval45,IKI_trial,fname,path,flag_tipo_norm,flag_media)
if nargin <9
    flag_media=0;
end

%% Normalizacion de los intervalos 
seq_per_block=size(interval12,2); %columnas
noBlock=size(interval12,1); %filas


%vectorizo para graficar
interval12=reshape(interval12',1, []);
interval23=reshape(interval23',1, []);
interval34=reshape(interval34',1, []);
interval45=reshape(interval45',1, []);

%plot
figure;
plot(interval12,'b.','MarkerSize',15); hold on;
plot(interval23,'r.','MarkerSize',15); hold on;
plot(interval34,'g.','MarkerSize',15); hold on;
plot(interval45,'m.','MarkerSize',15); hold on;

for i=1:seq_per_block:length(interval12)
    xline(i)
end
ylim([0 1.5])
title('crudo')
legend('Int-12: "41"','Int-23: "13"','Int-34: "32"','Int-45: "24"');

if ~strcmp(fname,'none')
    saveas(gcf,[path fname(1:end-4) '_Intervalos_crudo.' 'fig']);
    saveas(gcf,[path fname(1:end-4) '_Intervalos_crudo.' 'png']);
    %clf %clear figure
end
if strcmp(flag_tipo_norm,'01')
    % Devuelvo los intervalos a modo matricial
    interval12=reshape(interval12, seq_per_block, noBlock)';
    interval23=reshape(interval23, seq_per_block, noBlock)';
    interval34=reshape(interval34, seq_per_block, noBlock)';
    interval45=reshape(interval45, seq_per_block, noBlock)';

    aux12=interval12;
    aux23=interval23;
    aux34=interval34;
    aux45=interval45;

    for i=1:noBlock

        media12=nanmean(interval12(i,:));
        media23=nanmean(interval23(i,:));
        media34=nanmean(interval34(i,:));
        media45=nanmean(interval45(i,:));

        std12=nanstd(interval12(i,:));
        std23=nanstd(interval23(i,:));
        std34=nanstd(interval34(i,:));
        std45=nanstd(interval45(i,:));

        if std12~=0 && std23~=0 && std34~=0 && std45~=0
            %no uso la funcion zscore porque no contempla NaN's
            interval12(i,:)=(interval12(i,:)-media12)/std12;
            interval23(i,:)=(interval23(i,:)-media23)/std23;
            interval34(i,:)=(interval34(i,:)-media34)/std34;
            interval45(i,:)=(interval45(i,:)-media45)/std45; 
        end

    end

    %vectorizo para graficar
    interval12=reshape(interval12',1, []);
    interval23=reshape(interval23',1, []);
    interval34=reshape(interval34',1, []);
    interval45=reshape(interval45',1, []);
else % if normalización del tipo "cero": centrar la transición
    
    if flag_media==1
        %calculo el promedio por transicion
        delta12=nanmean(interval12);
        delta23=nanmean(interval23);
        delta34=nanmean(interval34);
        delta45=nanmean(interval45);
    else
        %calculo la mediana por transicion
        delta12=nanmedian(interval12);
        delta23=nanmedian(interval23);
        delta34=nanmedian(interval34);
        delta45=nanmedian(interval45);
    end
    
    interval12=interval12-delta12;
    interval23=interval23-delta23;
    interval34=interval34-delta34;
    interval45=interval45-delta45;
    
end

%plot
figure;
plot(interval12,'b.','MarkerSize',15); hold on;
plot(interval23,'r.','MarkerSize',15); hold on;
plot(interval34,'g.','MarkerSize',15); hold on;
plot(interval45,'m.','MarkerSize',15); hold on;

for i=1:seq_per_block:length(interval12)
    xline(i)
end
%ylim([0 1.5])
title('corrección')
legend('Int-12: "41"','Int-23: "13"','Int-34: "32"','Int-45: "24"');

if ~strcmp(fname,'none')
    saveas(gcf,[path fname(1:end-4) '_Intervalos_corr.' 'fig']);
    saveas(gcf,[path fname(1:end-4) '_Intervalos_corr.' 'png']);
    %clf %clear figure
end

% Devuelvo los intervalos a modo matricial
interval12=reshape(interval12, seq_per_block, noBlock)';
interval23=reshape(interval23, seq_per_block, noBlock)';
interval34=reshape(interval34, seq_per_block, noBlock)';
interval45=reshape(interval45, seq_per_block, noBlock)';

%% Normalización IKI_trial
IKI_trial=reshape(IKI_trial',1,[]);
if strcmp(flag_tipo_norm,'Z')
    % Normalización del IKI_trial con zscore (toda la tarea)
    IKI_trial_norm_Z=(IKI_trial-nanmedian(IKI_trial))/nanstd(IKI_trial);
    IKI_trial_norm=reshape(IKI_trial_norm_Z,seq_per_block,noBlock)';
elseif strcmp(flag_tipo_norm,'01')
    % Normalización del IKI_trial con 0-1 (toda la tarea)
  IKI_trial_norm_01=rescale(IKI_trial);
  IKI_trial_norm=reshape(IKI_trial_norm_01,seq_per_block,noBlock)';

else %cero
    for i=1:noBlock
        for j=1:seq_per_block
            if ~isnan(interval12(i,j))&& ~isnan(interval23(i,j))&& ~isnan(interval34(i,j))&& ~isnan(interval45(i,j)) %secuencia completa
                IKI_trial_norm(i,j)=nanmedian([interval12(i,j) interval23(i,j) interval34(i,j) interval45(i,j)]);
            else %secuencia incompleta
                if j>2 %no es la primera seq correcta del bloque, es la ultima y esta incompleta
                    if ~isnan(interval12(i,j))&& ~isnan(interval23(i,j))&& ~isnan(interval34(i,j))&& isnan(interval45(i,j)) %falta la ultima transicion
                        IKI_trial_norm(i,j)=nanmedian([interval12(i,j) interval23(i,j) interval34(i,j) interval45(i,j-1)]);
                    elseif ~isnan(interval12(i,j))&& ~isnan(interval23(i,j))&& isnan(interval34(i,j))&& isnan(interval45(i,j)) %faltan las dos ultimas transicion
                        IKI_trial_norm(i,j)=nanmedian([interval12(i,j) interval23(i,j) interval34(i,j-1) interval45(i,j-1)]);
                    elseif ~isnan(interval12(i,j))&& isnan(interval23(i,j))&& isnan(interval34(i,j))&& isnan(interval45(i,j)) %faltan las tres ultimas transicion
                        IKI_trial_norm(i,j)=nanmedian([interval12(i,j) interval23(i,j-1) interval34(i,j-1) interval45(i,j-1)]);
                    else %seq incorrecta
                        IKI_trial_norm(i,j)=NaN;
                    end
                else % si es la unica secuencia correcta (incompleta) del bloque
                    %van a ser nan las transiciones que no llegó a hacer.
                    IKI_trial_norm(i,j)=nanmedian([interval12(i,j) interval23(i,j) interval34(i,j) interval45(i,j)]);
                end
            end
        end
    end
    
end

