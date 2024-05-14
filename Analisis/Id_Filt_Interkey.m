%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 Analisis de los interkeys intervls: 
%%% Si hay un interkey interval muy distinto a los demas se le asigna NaN 
%%%
%%% La funcion recibe los intervalos correspondientes a las secuencias
%%% correctas, los de secuencias incorrectas tienen un NaN en su componente
%%% Devuelve los mismos intervalos pero sin los valores que son outliers.
%%% fname y path son parámetros que sirven para luego guardar las imágenes
%%% que pueden o no generarse. 
%%% Esto es para que cuando se realiza el promedio en Id_RunAnalysis o en
%%% Analisis_Seq_Individual estos valores no lo alteren. Este analisis 
%%% evita los saltos abruptos en los la curva de learning o en MOGS y MONGS
%%%
%%% NOTA: Este filtrado no se utilizó ya que la normalización fue
%%% suficiente. 
%%%                     
%%%                    Guadalupe Cutrera 26/12/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [filt12_block,filt23_block,filt34_block,filt45_block,cant_puntos_total]= Id_Filt_Interkey(interval12,interval23,interval34,interval45,fname, path)

seq_per_block=size(interval12,2); %columnas
noBlock=size(interval12,1); %filas

%% vectorizo los intervalos para poder graficarlos
interval12=reshape(interval12',1, []);
interval23=reshape(interval23',1, []);
interval34=reshape(interval34',1, []);
interval45=reshape(interval45',1, []);

figure;
%subplot(2,1,1)
plot(interval12,'b.','MarkerSize',15); hold on; 
plot(interval23,'r.','MarkerSize',15); hold on; 
plot(interval34,'g.','MarkerSize',15); hold on;
plot(interval45,'m.','MarkerSize',15); hold on;

for i=1:seq_per_block:length(interval12)
    xline(i)
end
%ylim([0 1.5])
title('Crudo');
legend('Int-12: "41"','Int-23: "13"','Int-34: "32"','Int-45: "24"');

if ~strcmp(fname,'none')
    saveas(gcf,[path fname(1:end-4) '_Intervalos_crudo.' 'fig']);
    saveas(gcf,[path fname(1:end-4) '_Intervalos_crudo.' 'png']);
    clf %clear figure
end
%% FILTRADO POR BLOQUE

filt12_block=interval12;
filt23_block=interval23;
filt34_block=interval34;
filt45_block=interval45;

%mad 3
%index=seq_per_block; %lo arranco en el segundo bloque
index=1;
for i=1:noBlock 
    
    %12
    aux12=filt12_block(index:(index+seq_per_block-1));
    mediana12=nanmedian(aux12);
    mad_12=mad(aux12);
    aux12(aux12>(mediana12+3*mad_12))=NaN;
    aux12(aux12<(mediana12-3*mad_12))=NaN;
    filt12_block(index:(index+seq_per_block-1))=aux12;

    %23
    aux23=filt23_block(index:(index+seq_per_block-1));
    mediana23=nanmedian(aux23);
    mad_23=mad(aux23);
    aux23(aux23>(mediana23+3*mad_23))=NaN;
    aux23(aux23<(mediana23-3*mad_23))=NaN;
    filt23_block(index:(index+seq_per_block-1))=aux23;

    %34
    aux34=filt34_block(index:(index+seq_per_block-1));
    mediana34=nanmedian(aux34);
    mad_34=mad(aux34);
    aux34(aux34>(mediana34+3*mad_34))=NaN;
    aux34(aux34<(mediana34-3*mad_34))=NaN;
    filt34_block(index:(index+seq_per_block-1))=aux34;
    
    %45
    aux45=filt45_block(index:(index+seq_per_block-1));
    mediana45=nanmedian(aux45);
    mad_45=mad(aux45);
    aux45(aux45>(mediana45+3*mad_45))=NaN;
    aux45(aux45<(mediana45-3*mad_45))=NaN;
    filt45_block(index:(index+seq_per_block-1))=aux45;
   
    index=index+seq_per_block;
    
end

figure;
plot(filt12_block,'b.','MarkerSize',15); hold on;
plot(filt23_block,'r.','MarkerSize',15); hold on;
plot(filt34_block,'g.','MarkerSize',15); hold on;
plot(filt45_block,'m.','MarkerSize',15); hold on;

for i=1:seq_per_block:length(interval12)
    xline(i)
end
%ylim([0 1.5])
title('filtro por bloque')
legend('Int-12: "41"','Int-23: "13"','Int-34: "32"','Int-45: "24"');

if ~strcmp(fname,'none')
    saveas(gcf,[path fname(1:end-4) '_Intervalos_filt.' 'fig']);
    saveas(gcf,[path fname(1:end-4) '_Intervalos_filt.' 'png']);
    clf %clear figure
end
%% Cuento la cantidad de puntos filtrados

%isnan devuleve un vector con 1 en las pos donde hay un nan y un 0 donde hay un numero real
nan12=isnan(interval12);
nan23=isnan(interval23);
nan34=isnan(interval34);
nan45=isnan(interval45);

nan_filt12=isnan(filt12_block);
nan_filt23=isnan(filt23_block);
nan_filt34=isnan(filt34_block);
nan_filt45=isnan(filt45_block);

%nan actuales - nan originales = nan nuevos
cant_puntos12=sum(nan_filt12)-sum(nan12);
cant_puntos23=sum(nan_filt23)-sum(nan23);
cant_puntos34=sum(nan_filt34)-sum(nan34);
cant_puntos45=sum(nan_filt45)-sum(nan45);

cant_puntos_total=cant_puntos12+cant_puntos23+cant_puntos34+cant_puntos45;

%me quedo con las posiciones de los puntos filtrados
pos_puntos12=(nan_filt12-nan12).*interval12;
pos_puntos23=(nan_filt23-nan23).*interval23;
pos_puntos34=(nan_filt34-nan34).*interval34;
pos_puntos45=(nan_filt45-nan45).*interval45;

pos_puntos12(pos_puntos12==0)=NaN;
pos_puntos23(pos_puntos23==0)=NaN;
pos_puntos34(pos_puntos34==0)=NaN;
pos_puntos45(pos_puntos45==0)=NaN;


figure;
%originales
plot(interval12,'b.','MarkerSize',15); hold on;
plot(interval23,'r.','MarkerSize',15); hold on;
plot(interval34,'g.','MarkerSize',15); hold on;
plot(interval45,'m.','MarkerSize',15); hold on;

%cruces negras son los puntos filtrados
plot(pos_puntos12,'x','Color',[0 0 0],'MarkerSize',15); hold on;
plot(pos_puntos23,'x','Color',[0 0 0],'MarkerSize',15); hold on;
plot(pos_puntos34,'x','Color',[0 0 0],'MarkerSize',15); hold on;
plot(pos_puntos45,'x','Color',[0 0 0],'MarkerSize',15); hold on;

for i=1:seq_per_block:length(interval12)
    xline(i)
end
%ylim([0 1.5])
hold on;
txt=['Quita ',num2str(cant_puntos_total), ' puntos'];
text(90,1.55,txt,'FontSize',10)
legend('Int-12: "41"','Int-23: "13"','Int-34: "32"','Int-45: "24"','filt');

if ~strcmp(fname,'none')
    saveas(gcf,[path fname(1:end-4) '_Intervalos_filt_completo.' 'fig']);
    saveas(gcf,[path fname(1:end-4) '_Intervalos_filt_completo.' 'png']);
end
clf %clear figure

%los vuelvo a poner en orden matricial para devolverlos en la función
filt12_block=reshape(filt12_block,seq_per_block, noBlock)';
filt23_block=reshape(filt23_block,seq_per_block, noBlock)';
filt34_block=reshape(filt34_block,seq_per_block, noBlock)';
filt45_block=reshape(filt45_block,seq_per_block, noBlock)';
end