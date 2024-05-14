function Learning_print(iki_mean,IKI_per_subj,path,titulo)
% Guadalupe Cutrera 28/12/23
% Funcion para comparar cada sujeto con la media grupal
durBlock=10;
durRest=3;
for j=1:size(IKI_per_subj,1)
    figure; set(gcf,'Color','white'); box OFF; hold on;
    
    contBlock=0;
    for i=1:durBlock:length(iki_mean)
        %media
        plot(i+contBlock*durRest:i+(durBlock-1)+contBlock*durRest,iki_mean(i:i+(durBlock-1)),'k','LineWidth',0.5);
        hold on;
        
        
        %sujeto individual
        plot(i+contBlock*durRest:i+(durBlock-1)+contBlock*durRest,IKI_per_subj(j,i:i+(durBlock-1)),'b','LineWidth',0.5);
        ylim([0 2.2]);
        
        contBlock=contBlock+1;
    end
    
    saveas(gcf,[path '_Learning_visual_' titulo '_Sujeto_' num2str(j) '.fig']);
    saveas(gcf,[path '_Learning_visual_' titulo '_Sujeto_' num2str(j) '.png']);
end
