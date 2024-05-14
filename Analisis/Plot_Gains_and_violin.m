%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Función para proyecto. Imprime las ganancias juntas en un subplot y los
%%% violines en otro.                                                   %%%
%%%                         Guadalupe Cutrera 5/3/24                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plot_Gains_and_violin(Matriz_Gains_mean,Matriz_Gains_std,Matriz_Gains_acum)
%Matriz_Gains_mean=[TL; MOGS; MONGS ] filas
%Matriz_Gains_acum =[TL; MOGS; MONGS] filas

figure;set(gcf,'Color','white'); box OFF; hold on; %sgtitle(titulo)
set(gcf,'Position',get(0,'ScreenSize'));

subplot(1,2,1)
%MOGS
x_vector = [1:length(Matriz_Gains_mean(2,:)), fliplr(1:length(Matriz_Gains_mean(2,:)))];
patch = fill(x_vector, [Matriz_Gains_mean(2,:) + Matriz_Gains_std(2,:) , fliplr(Matriz_Gains_mean(2,:) - Matriz_Gains_std(2,:))], [219 123 128]./255);
set(patch, 'edgecolor', 'none'); set(patch, 'FaceAlpha', 0.5);
hold on;
plot(Matriz_Gains_mean(2,:),'r','LineWidth',1); 
hold on;

%MONGS
x_vector = [1:length(Matriz_Gains_mean(3,:)), fliplr(1:length(Matriz_Gains_mean(3,:)))];
patch = fill(x_vector, [Matriz_Gains_mean(3,:) + Matriz_Gains_std(3,:) , fliplr(Matriz_Gains_mean(3,:) - Matriz_Gains_std(3,:))], [128 193 219]./255);
set(patch, 'edgecolor', 'none'); set(patch, 'FaceAlpha', 0.5);
hold on;
plot(Matriz_Gains_mean(3,:),'b','LineWidth',1); 
hold on;
ylim([-1.5 1.5])
yline(0)
xlabel('Bloque','FontName','Arial','FontSize',15)
ylabel('Suma de deltas (s)','FontName','Arial','FontSize',15); 
xlim([0.7 11.3])

%TL
x_vector = [1:length(Matriz_Gains_mean(1,:)), fliplr(1:length(Matriz_Gains_mean(1,:)))];
patch = fill(x_vector, [Matriz_Gains_mean(1,:) + Matriz_Gains_std(1,:) , fliplr(Matriz_Gains_mean(1,:) - Matriz_Gains_std(1,:))],[160 160 160]./255);
set(patch, 'edgecolor', 'none'); set(patch, 'FaceAlpha', 0.5);
hold on;
plot(Matriz_Gains_mean(1,:),'k','LineWidth',1); 
hold on;



subplot(1,2,2)
%violines
violin(Matriz_Gains_acum','xlabel',{'TL','MOGS','MONGS'},'facecolor',[[160 160 160]./255;[219 123 128]./255; [128 193 219]./255],'edgecolor','k','mc','k','ylabel','Sum of deltas (s)','plotlegend',0,'medc',[],'facealpha',0.7);
plot(Matriz_Gains_acum,'k.','Markersize',15)
yline(0)
ylim([-2.5 2.5])