%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MOGS: MICRO OFFLINE GAINS. Diferencia entre el iki de la ultima     %%%
%%% secuencia del bloque con el iki de la primer secuencia del          %%%
%%% siguiente bloque                                                    %%%
%%%                                                                     %%%
%%% MONGS: MICRO ONLINE GAINS. Diferncia entre el iki de la primer      %%%
%%% secuencia del bloque con el iki de la ultima secuencia del mismo    %%%
%%% bloque                                                              %%%
%%%                                                                     %%%
%%% TOTAL LEARNING: diferencia entre el iki de la primer secuencia del  %%%
%%% bloque con el iki de la primera secuencia del siguiente bloque.     %%%
%%%                      Guadalupe cutrera 21/7/23                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MOGS,MONGS,Total_Learning]=Micro_gains_key(IKI_per_trial,noBlock)
%% calculo de los MOGS

for j=2:(noBlock)
    ultimo=size(IKI_per_trial,2); %cantidad de secuencias por bloque 
    primero=1;
    %los while son para que me de los IKI de secuencias correctas
    while isnan(IKI_per_trial(j,primero))
        primero=primero+1;
        if primero>size(IKI_per_trial,2) %si no hay secuencia correcta, se restará con el NaN  
            primero=primero-1;
            break;
        end
    end
    while isnan(IKI_per_trial(j-1,ultimo))
        ultimo=ultimo-1; 
        if 1>ultimo                     %si no hay secuencia correcta, se restará con el NaN
            ultimo=ultimo+1;
            break;
        end
    end
    
    mogs(j-1)=IKI_per_trial(j-1,ultimo)-IKI_per_trial(j,primero);
    
end

MOGS=mogs';

%% Calculo de MONGS

for j=1:(noBlock)
    ultimo=size(IKI_per_trial,2); %cantidad de secuencias por bloque 
    primero=1;
    %los while son para que me de los IKI de secuencias correctas
    while isnan(IKI_per_trial(j,primero))
        primero=primero+1;
        if primero>size(IKI_per_trial,2) %si no hay secuencia correcta, se restará con el NaN  
            primero=primero-1;
            break;
        end
    end
    while isnan(IKI_per_trial(j,ultimo))
        ultimo=ultimo-1; 
        if 1>ultimo                     %si no hay secuencia correcta, se restará con el NaN
            ultimo=ultimo+1;
            break;
        end
    end
    
    if ultimo==primero %una sola secuencia en ese bloque
        mongs(j)=NaN;
    else
        mongs(j)=IKI_per_trial(j,primero)-IKI_per_trial(j,ultimo);
    end
    
end

MONGS=mongs;

%% Calculo Total Learning

for j=2:(noBlock)
    PrimeroA=1;
    PrimeroB=1;
    %los while son para que me de los IKI de secuencias correctas
    while isnan(IKI_per_trial(j-1,PrimeroA))
        PrimeroA=PrimeroA+1;
        if PrimeroA>size(IKI_per_trial,2) %si no hay secuencia correcta, se restará con el NaN  
            PrimeroA=PrimeroA-1;
            break;
        end
    end
    while isnan(IKI_per_trial(j,PrimeroB))
        PrimeroB=PrimeroB+1;
        if PrimeroB>size(IKI_per_trial,2) %si no hay secuencia correcta, se restará con el NaN  
            PrimeroB=PrimeroB-1;
            break;
        end
    end
    
    total_L(j-1)=IKI_per_trial(j-1,PrimeroA)-IKI_per_trial(j,PrimeroB);
    
end

Total_Learning=total_L;

end

