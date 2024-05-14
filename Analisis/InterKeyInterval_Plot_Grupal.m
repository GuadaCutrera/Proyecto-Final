%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     InterKeyInterval_Plot_Grupal                    %%%
%%% Esta funcion imprime una figura en la que se visualiza la media y el %%
%%% error estandar del in terkey interval. Recibe la media y error.     %%%
%%%                                                                     %%%
%%% paradigm_flag denota si el experimenti se realizó para sujetos que  %%%
%%% terminan el bloque por tiempo o por cantidad de teclas.             %%%
%%% Group_Parameters es un struct que entre otras cosas guarda la cantidad%
%%% máxima de secuencias en un bloque y el tipo de procesamiento        %%%
%%% titulo: Si el protocolo es con ITI o sin ITI, si la data a visualizar %
%%% es cruda o tiene algún procesamiento realizado. Ej: "With ITI - filt" %
%%%                                                                     %%%
%%%                   Guadalupe Cutrera 3/8/23                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function InterKeyInterval_Plot_Grupal(IKI_mean, IKI_std,paradigm_flag,Group_Parameters,titulo)


figure; set(gcf,'Color','white'); box OFF; hold on; sgtitle(titulo)
set(gcf,'Position',get(0,'ScreenSize'));

barra=max(IKI_mean(:)+0.1);

if strcmp(paradigm_flag,'teclas')==1
    
    for i=1:Group_Parameters.cant_SeqBlock:length(IKI_mean)
        %error estandar
        x_vector = [i:i+(Group_Parameters.cant_SeqBlock-1), fliplr(i:i+(Group_Parameters.cant_SeqBlock-1))];
        x_vector2=[IKI_mean(i:i+(Group_Parameters.cant_SeqBlock-1))+IKI_std(i:i+(Group_Parameters.cant_SeqBlock-1)),fliplr(IKI_mean(i:i+(Group_Parameters.cant_SeqBlock-1))-IKI_std(i:i+(Group_Parameters.cant_SeqBlock-1)))];
        patch = fill(x_vector, x_vector2, [128 193 219]./255);
        set(patch, 'edgecolor', 'none');
        set(patch, 'FaceAlpha', 1);
        hold on;

        %media
        plot(i:i+(Group_Parameters.cant_SeqBlock-1),IKI_mean(i:i+(Group_Parameters.cant_SeqBlock-1)),'k','LineWidth',0.5);
        hold on;

        %barras de rest
        bar(i+(Group_Parameters.cant_SeqBlock-1)+0.5,barra,1,'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7]);
        hold on;
    end
    
else %Cuando corta por tiempo
    
    for i=1:Group_Parameters.cant_SeqBlock:length(IKI_mean)
        %error estandar
        x_vector = [i:(i+Group_Parameters.cant_SeqBlock-1), fliplr(i:(i+Group_Parameters.cant_SeqBlock-1))];
        x_vector2=[IKI_mean(i:(i+Group_Parameters.cant_SeqBlock-1))+IKI_std(i:(i+Group_Parameters.cant_SeqBlock-1)),fliplr(IKI_mean(i:(i+Group_Parameters.cant_SeqBlock-1)) - IKI_std(i:(i+Group_Parameters.cant_SeqBlock-1)))];
        
        %este contador me va a permitir no graficar los errores de los NaNs y poder cerrar los patch para rellenarlos
        cont=1;
        while cont<=length(x_vector2) && ~isnan(x_vector2(cont)) 
            cont=cont+1; 
        end
        cont=cont-1;
        x_vector2(isnan(x_vector2))=[];
        if cont >= Group_Parameters.cant_SeqBlock
            cont=cont/2;
        end
        x_vector=[x_vector(1:cont) x_vector((end-cont+1):end)];
                
        patch = fill(x_vector, x_vector2, [128 193 219]./255);
        set(patch, 'edgecolor', 'none');
        set(patch, 'FaceAlpha', 1);
       
        hold on;

        %media
        plot(i:i+Group_Parameters.cant_SeqBlock-1,IKI_mean(i:(i+Group_Parameters.cant_SeqBlock-1)),'k','LineWidth',0.5);
        hold on;

        %barras de rest
        bar(i+(Group_Parameters.cant_SeqBlock-1)+0.5,barra,1,'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7]);
        hold on;
    end
    
end

xlabel('Blocks','FontName','Arial','FontSize',12)
ylabel('Interkeys interval','FontName','Arial','FontSize',12);
if  contains(titulo,'crudo') || Group_Parameters.flag_norm==0
       ylim([min(IKI_mean-0.1) max(IKI_mean+0.1)])
else
    if strcmp(Group_Parameters.flag_tipo_norm,'Z')
        ylim([-3 3])
    else %norm 01
        ylim([0 1])
    end
end
xlim([0 length(IKI_mean)+1]);