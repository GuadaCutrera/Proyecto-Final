function [MicroMogs_acum,MicroMongs_acum,MicroTL_acum]= MicroMicroGains_Plot_Individual(MicroMogs,MicroMongs,MicroTL,titulo,seq_results)
%Guadalupe Cutrera 22/2/23
figure; set(gcf,'Color','white'); box OFF; hold on; sgtitle(titulo)

%%  ------------------ MICRO MICRO No acumulativo--------------------------
aux_mogs=reshape(MicroMogs',1, []);
aux_mongs=reshape(MicroMongs',1, []);
aux_TL=reshape(MicroTL',1, []);

% MicroMogs
subplot(2,3,1)
plot(aux_mogs,'r.')
for i=1:(max(seq_results.cant_SeqBlock)-1):length(aux_mogs)
    xline(i)
end
yline(0)
if seq_results(1,1).flag_norm==0 || contains(titulo,'crudo')
    ylim([-1.5 1.5])
else
    ylim([-15 15])
end
title('No Acum.')
ylabel('Micro MOGS')


% MicroMongs
subplot(2,3,2)
plot(aux_mongs,'b.')
for i=1:max(seq_results.cant_SeqBlock):length(aux_mongs)
    xline(i)
end
yline(0)
if seq_results(1,1).flag_norm==0 || contains(titulo,'crudo')
    ylim([-1.5 1.5])
else
    ylim([-15 15])
end
xlim([0 length(aux_mongs)])
title('No Acum.')
ylabel('Micro MONGS')

%Micro TL
subplot(2,3,3)
plot(aux_TL,'k.')
for i=1:max(seq_results.cant_SeqBlock):length(aux_TL)
    xline(i)
end
yline(0)
if seq_results(1,1).flag_norm==0 || contains(titulo,'crudo')
    ylim([-1.5 1.5])
else
    ylim([-15 15])
end
xlim([0 length(aux_TL)])
title('No Acum.')
ylabel('Micro TL')

%% ------------------------MICRO MICRO Acumulativo ------------------------
        %Se hace la acumulación por bloque, no de la tarea completa.

%MicroMogs
subplot(2,3,4)
MicroMogs_acum=[]
for i=1:size(MicroMogs,1)
     MicroMogs_acum(i,:)=cumsum(MicroMogs(i,:),'omitnan');
     plot(((i-1)*(seq_results.cant_SeqBlock-1)+1):((i-1)*(seq_results.cant_SeqBlock-1)+1+(seq_results.cant_SeqBlock-2)),MicroMogs_acum(i,:),'r.');
     xline((i-1)*(seq_results.cant_SeqBlock-1)+1); hold on; %la línea indica el comienzo del bloque
end
title('Acum.')
ylabel('Micro MOGS')
if seq_results(1,1).flag_norm==0 || contains(titulo,'crudo')
    ylim([-3 3])
else
    ylim([-15 15])
end
yline(0)

% MicroMongs
subplot(2,3,5)
MicroMongs_acum=[];
for i=1:size(MicroMongs,1)
     MicroMongs_acum(i,:)=cumsum(MicroMongs(i,:),'omitnan');
     plot(((i-1)*seq_results.cant_SeqBlock+1):((i-1)*seq_results.cant_SeqBlock+1+(seq_results.cant_SeqBlock-1)),MicroMongs_acum(i,:),'b.');
     xline((i-1)*seq_results.cant_SeqBlock+1); hold on; %la línea indica el comienzo del bloque
end
title('Acum.')
ylabel('Micro MONGS')
if seq_results(1,1).flag_norm==0 || contains(titulo,'crudo')
    ylim([-3 3])
else
    ylim([-15 15])
end
xlim([0 length(aux_mongs)])
yline(0)

%MicroTL
subplot(2,3,6)
MicroTL_acum=[];
for i=1:size(MicroTL,1)
     MicroTL_acum(i,:)=cumsum(MicroTL(i,:),'omitnan');
     plot(((i-1)*(seq_results.cant_SeqBlock-1)+1):((i-1)*(seq_results.cant_SeqBlock-1)+1+(seq_results.cant_SeqBlock-2)),MicroTL_acum(i,:),'k.');
     xline((i-1)*seq_results.cant_SeqBlock+1); hold on; %la línea indica el comienzo del bloque
end
title('Acum.')
ylabel('Micro TL')
if seq_results(1,1).flag_norm==0 || contains(titulo,'crudo')
    ylim([-3 3])
else
    ylim([-15 15])
end
xlim([0 length(aux_TL)])
yline(0)
