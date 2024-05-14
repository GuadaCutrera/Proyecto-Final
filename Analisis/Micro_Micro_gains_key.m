%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         MICRO MICRO GAINS                           %%%
%%% This function allows cuantification of gains in between trials      %%%
%%% In "Wiht ITI" protocol it helps understand where the true gain is.  %%%
%%%                                                                     %%%
%%%                      Guadalupe cutrera 5/8/23                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MicroMOGS,MicroMONGS,MicroTL] = Micro_Micro_gains_key(interval12,interval45,flag_continuous_seq)

noBlock=size(interval12,1);
noSeq=size(interval12,2);
MicroMOGS=NaN(noBlock,noSeq-1); 
MicroTL=NaN(noBlock,noSeq-1); 
MicroMONGS=NaN(noBlock,noSeq); 

for i=1:noBlock
    for j=1:noSeq %la primera no tiene una secuencia previa y vale NaN siempre
        if  j>1 && flag_continuous_seq(i,j)== 1 
            MicroMOGS(i,j-1)=interval45(i,j-1)- interval12(i,j);
            MicroTL(i,j-1)=interval12(i,j-1)-interval12(i,j);
        end %END IF
        MicroMONGS(i,j)=interval12(i,j)-interval45(i,j);
    end %END FOR J
end %END FOR I


