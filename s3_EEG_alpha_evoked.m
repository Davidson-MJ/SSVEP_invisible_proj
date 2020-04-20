% s3_EEG_alphalateralization
clear all, close all; clc
%% UPDATED for Wills data:
try cd('/Users/MDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA')
catch
end
basefol=pwd;
cd('EEG')
dirs = dir([pwd filesep '*_*' 'EEG']);

%Epoch by number of targets involved, and direction (disapp and reapp).
% also stores duration per event.
% also stores cumularive BP (single array) per event, for next job plots.

dbstop if error
%%

% epoch per side (left/right disappearance).

job.epochperppant_PFI=0; % epochs as disappearances by location
job.epochperppant_PMD=0; % as above.
   
job.participant_FilterandHilbert=0; % compute hilbert envelope.

job.combineAcrossppants=0; % combine per channel, across participants.


%% >>>>>> actual plot of hilbert results:
job.plotHILBresults =1;
%%
%load common variables:
window=[-3 3];
srate=250;
allppants=[1,2,4,6,7,9:19]; % new
getelocs;


%%%%%%% function calls:
                                % preprocessing
% epoch EEG left vs Right
if job.epochperppant_PFI==1
   epochperppantLeftvsRight(basefol,dirs, allppants)
end

if job.epochperppant_PMD==1
   epochperppantLeftvsRight_PMD;%(basefol,dirs, allppants)
end


%%                              % analysis
%crunch Time-frequency decomposition using multitaper method.

if job.participant_FilterandHilbert==1 %filter at alpha range     
    Ppant_hilbert_LeftvsRight(basefol,dirs,allppants)          
end
 

if job.combineAcrossppants==1       
       combineLRdataacrossppants(basefol,dirs,allppants)  
end
%% 
%%