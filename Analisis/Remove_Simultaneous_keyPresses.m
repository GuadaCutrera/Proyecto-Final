%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    Remove_Simultaneous_keyPresses                   %%%
%%% This function helps in cases where the subject presses more than one %%
%%% key on the keyboard at the same time. In this scenario, the script  %%%
%%% Analisis_Seq_Individual (time or key paradigm), cant tell whether the %
%%% sequence is correct or not.                                         %%%
%%% We remove simultaneous keypresses to prevent this issue.            %%%
%%%                                                                     %%%        
%%%     Input: logoriginal: cell that contains in 3 spaces the information%  
%%%     of what the subject did on their task.                          %%%
%%%                                                                     %%%
%%%     Output: same cell but with repeated keypresses removed          %%%
%%%                                                                     %%%    
%%%                                                                     %%%
%%%                     Guadalupe Cutrera 16/7/2023                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function new_logoriginal=Remove_Simultaneous_keyPresses(logoriginal)

index=1;

for nLine=1:length(logoriginal)
    if nLine<6 || str2double(logoriginal{nLine}{1})-str2double(logoriginal{nLine-1}{1})>0
        new_logoriginal{index}{1}=logoriginal{nLine}{1};
        new_logoriginal{index}{2}=logoriginal{nLine}{2};
        if size(logoriginal{nLine},2)==3
            new_logoriginal{index}{3}=logoriginal{nLine}{3};
        end
        index=index+1;
    end
end
    

