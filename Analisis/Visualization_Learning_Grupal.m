%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         FIGURE                                      %%%
%%% Learning Curve and Tapping Speed visualization                      %%%
%%%                                                                     %%%
%%%                      23/6/2023                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Visualization_Learning_Grupal(iki_mean,iki_std,variable,titulo)

durBlock=10;                                                              % cantidad de puntos dentro de un bloque
durRest=3;                                                                 %separacion visual entre bloques (cantidad de puntos a intercalar)

figure; set(gcf,'Color','white'); box OFF; hold on; 
set(gcf,'Position',get(0,'ScreenSize'));

%barra=max(iki_mean(:)+0.1);
contBlock=0;
for i=1:durBlock:length(iki_mean)
    %error estandar
    x_vector = [i+contBlock*durRest:i+(durBlock-1)+contBlock*durRest, fliplr(i+contBlock*durRest:i+(durBlock-1)+contBlock*durRest)];
    x_vector2=[iki_mean(i:i+(durBlock-1))+iki_std(i:i+(durBlock-1)),...
        fliplr(iki_mean(i:i+(durBlock-1))-iki_std(i:i+(durBlock-1)))];
    patch = fill(x_vector, x_vector2, [190 190 190]./255); %[128 193 219]./255
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', 1);
    hold on;

    %media
    plot(i+contBlock*durRest:i+(durBlock-1)+contBlock*durRest,iki_mean(i:i+(durBlock-1)),'k','LineWidth',0.5);
    hold on;
    contBlock=contBlock+1;

    %barras de rest
    %bar(i+(durBlock-1)+0.5,barra,1,'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7]);
    %hold on;
end
xlabel('Blocks');
%ylim([0 2.2])
if strcmp(variable,'speed')
    ylabel('Tapping Speed (keypresses/s)')
else 
    ylabel('Interkey Interval (s)')
end
sgtitle(titulo);
xticks([5:durBlock+durRest:395]); 
xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20',...
    '21','22','23','24','25','26','27','28'})