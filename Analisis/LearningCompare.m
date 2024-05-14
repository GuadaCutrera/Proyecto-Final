function LearningCompare(Matriz_Grupos,Matriz_std,flag_cantidad)

durBlock=10;                                                              % cantidad de puntos dentro de un bloque
durRest=3;     
contBlock=0;

figure; set(gcf,'Color','white'); box OFF; hold on; sgtitle('Curva de Aprendizaje')
set(gcf,'Position',get(0,'ScreenSize'));

for i=1:durBlock:size(Matriz_Grupos,2)
    %error estandar
    x_vector = [i+contBlock*durRest:i+(durBlock-1)+contBlock*durRest, fliplr(i+contBlock*durRest:i+(durBlock-1)+contBlock*durRest)];
    patch = fill(x_vector, [Matriz_Grupos(1,i:i+(durBlock-1)) + Matriz_std(1,i:i+(durBlock-1)) , fliplr(Matriz_Grupos(1,i:i+(durBlock-1)) - Matriz_std(1,i:i+(durBlock-1)))], [0.8500 0.3250 0.0980]);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', 0.3);
    hold on;
    %media
    plot(i+contBlock*durRest:i+(durBlock-1)+contBlock*durRest,Matriz_Grupos(1,i:i+(durBlock-1)),'Color',[0.8500 0.3250 0.0980]); 
    hold on;
    %error estandar
    x_vector = [i+contBlock*durRest:i+(durBlock-1)+contBlock*durRest, fliplr(i+contBlock*durRest:i+(durBlock-1)+contBlock*durRest)];
    patch = fill(x_vector, [Matriz_Grupos(2,i:i+(durBlock-1)) + Matriz_std(2,i:i+(durBlock-1)) , fliplr(Matriz_Grupos(2,i:i+(durBlock-1)) - Matriz_std(2,i:i+(durBlock-1)))],[0.4940 0.1840 0.5560]);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', 0.3);
    hold on;
    %media
    plot(i+contBlock*durRest:i+(durBlock-1)+contBlock*durRest,Matriz_Grupos(2,i:i+(durBlock-1)),'Color',[0.4940 0.1840 0.5560]);
    hold on;
    if flag_cantidad>2
         %error estandar
        x_vector = [i+contBlock*durRest:i+(durBlock-1)+contBlock*durRest, fliplr(i+contBlock*durRest:i+(durBlock-1)+contBlock*durRest)];
        patch = fill(x_vector, [Matriz_Grupos(3,i:i+(durBlock-1)) + Matriz_std(3,i:i+(durBlock-1)) , fliplr(Matriz_Grupos(3,i:i+(durBlock-1)) - Matriz_std(3,i:i+(durBlock-1)))], [169 243 114]./255);
        set(patch, 'edgecolor', 'none');
        set(patch, 'FaceAlpha', 0.3);
        hold on;
        %media
        plot(i+contBlock*durRest:i+(durBlock-1)+contBlock*durRest,Matriz_Grupos(3,i:i+(durBlock-1)),'g');
    end
    if flag_cantidad==4
        %error estandar
        x_vector = [i+contBlock*durRest:i+(durBlock-1)+contBlock*durRest, fliplr(i+contBlock*durRest:i+(durBlock-1)+contBlock*durRest)];
        patch = fill(x_vector, [Matriz_Grupos(4,i:i+(durBlock-1)) + Matriz_std(4,i:i+(durBlock-1)) , fliplr(Matriz_Grupos(4,i:i+(durBlock-1)) - Matriz_std(4,i:i+(durBlock-1)))], [[219 123 128]]./255);
        set(patch, 'edgecolor', 'none');
        set(patch, 'FaceAlpha', 0.3);
        hold on;
        %media
        plot(i+contBlock*durRest:i+(durBlock-1)+contBlock*durRest,Matriz_Grupos(4,i:i+(durBlock-1)),'m');
    end
    contBlock=contBlock+1;
end

if flag_cantidad==3
    %legend('Without ITI', 'With ITI 1s', 'With ITI 1.5s');
elseif flag_cantidad==4
    legend('','Int12','', 'Int23','','Int34','', 'Int45')
else
    %legend('With ITI 1.5s Key', 'With ITI 1.5s Time')
    legend('','Grupo sin ITI','', 'Grupo con ITI 1s','FontName','Arial','FontSize',12)
end

xlim([0 364]) 
xlabel('Bloque','FontName','Arial','FontSize',12)
ylabel('Intervalos entre teclas (s)','FontName','Arial','FontSize',12)

xticks([5:13:364]); 
xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20',...
    '21','22','23','24','25','26','27','28'})

