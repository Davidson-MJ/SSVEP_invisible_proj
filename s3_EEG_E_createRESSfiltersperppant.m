% Motivation is to create unique RESS spatial filters per participant and
% frequency combo.


%we need all the trials in which the targets were present,


clear all

%% UPDATED for Wills data:

try cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA')
%     addpath(
catch
end
basefol=pwd;

clearvars -except basefol allppants
dbstop if error
%%

job.concatLiketrialsperppant=1; %also performs RESS, we are left with RESS filters per freq per ppant.

% SET UP RESS params.

peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; % don't change!




%set up directories and dependencies.
getelocs;

cd([basefol filesep 'EEG'])
eegfol=pwd;

pdirs = dir([pwd filesep '*_*' 'EEG']);

%% load data to determine physical catch timing.


% allppants=1:19;

%remaining participants after behavioral data analysis/exclusion
% allppants=[1,2,4,6,9:16,18]; %
allppants=[1,2,4,6,7,9:19]; % NEW ppants.
%newppants =[7,17,19]
window=[-3 3];

srate=250;

epochdur = sum(abs(window))*srate;

timeid = [0:1/srate:epochdur];
timeid= timeid-3;

onsetc = ceil(epochdur)/2;
% peakfreqsare=[20,40]; %hz
%timing

%%
if job.concatLiketrialsperppant==1
    
    
  
    %construct neighbourhood parameters.
    neighbour=[];
    for ihz = 1:length(peakfreqsare)
    
        hzAT = peakfreqsare(ihz);        
        %needs to be narrow to avoid other stim freqs.  
    neighbour(ihz).centrefreqs=[hzAT-3 , hzAT+3];   % a low and high centre freq for covR
    neighbour(ihz).sd=[1, 1]; %FWHM of those...

    end

    %%
    windowsmall=[];
    windowsmall(1,:) = [-3 -1] ;%  window targets present
    windowsmall(2,:) = [1 3] ;%  window targets present after BP
    %%
%     usechans=[51:64, 17:32]; %occipital
    usechans=[1:64]; %occipital
    
    
    for ifol =allppants
        
        
        %%
        cd(eegfol)
        cd(pdirs(ifol).name)
        
        
        %% load the relevant PFI data.
        load('ppant_PFI_Epoched');
        load('ppant_Catch_Epoched')
        
        
        %%       
                % collect relevant trials for each type of spatial
                % configuration/filter construction.
              
                %%
                for id=1:12% Use all epochs so as not to bias condition comparisons.                    
                    switch id
                        case 1
                            dataIN=ppant_SNREEG_PFI_0_1;
                            usewindow=1:2; 
                        case 2
                            dataIN=ppant_SNREEG_PFI_1_0;
                            usewindow=1:2; 
                        case 3
                            dataIN=ppant_SNREEG_PFI_1_2;
                            usewindow=1:2; 
                        case 4
                            dataIN=ppant_SNREEG_PFI_2_1;

                            usewindow=1:2; 
                        case 5
                            dataIN=ppant_SNREEG_PFI_3_2;

                            usewindow=1:2;
                        case 6
                            dataIN=ppant_SNREEG_PFI_2_3;                            

                            usewindow=1:2;
                        case 7
                            dataIN=ppant_SNREEG_PFI_4_3;

                            usewindow=1:2;
                        case 8
                            dataIN=ppant_SNREEG_PFI_3_4;                            

                            usewindow=1:2;
                        
                        case 9
                            dataIN=ppant_SNREEG_disapBPwithincatch;

                            usewindow=1; 
                        case 10
                            dataIN=ppant_SNREEG_reapBPaftercatch;

                            usewindow=2; 
                        case 11
                            dataIN=ppant_SNREEG_catchramponset;

                            usewindow=1;
                        case 12
                            dataIN=ppant_SNREEG_catchrampoffset;

                            usewindow=2;
                        
                        
                            
                    end
                    
                    trialtypesperTargPresent= 1:size(dataIN,1);
                    
                    % now that we have a datatype, get the right epochs that share stimulus configuration
                    
                    % also correct window (pre or post indication of target
                    % presence)
                    for windch=1:length(usewindow)
                        wind=usewindow(windch);
                    windcheck= windowsmall(wind,:);
                    
                    tidx=dsearchn(timeid', [windcheck]');
                    
                    
                    %reduce size.
                    datast= squeeze(dataIN(trialtypesperTargPresent,:,tidx(1):tidx(2)));
                    
                    durscheck=[];
                    %% check for bad trials (noisy)
                    %std per trial(average over electrodes)
                    tmp=[];
                    if ndims(datast)<3
                        tmp(1,:,:)= datast;
                        datast=tmp;
                    end
                    %SD per trial
                    datastSD = std(squeeze(mean(datast,2)),0,2);
                    %SD per chan.
                    datastSDchan = mean(std(datast,0 ,3),1);
                    
                    %remove those with 2.5*std from mean.                    
                    mSD=nanmean(datastSD);
                    trialSD=nanstd(datastSD);                    
                    
                    %remove the trials with excessive noise.
                    badtrials=find(datastSD>mSD+2.5*trialSD)';
                                        
                    % also skip the trials which were only transient button
                    % presses. (less than one second).                    
                    shorttrials = find(durscheck<60);
                    
                    %also find out if there were any nan!
                    nantrials = find(isnan(datastSD))';
                                       
                    allbadtrials = [badtrials, shorttrials,nantrials];
                    
                    % remove these from consideration.
                    keeptrials=1:size(datast,1);                    
                    keeptrials(allbadtrials)=[];
                    datast=datast(keeptrials,:,:);
                    
                    %sanity check this is clean data:
%                     figure(1); subplot(211); plot(squeeze(datast(:,62,:))');
%                     subplot(212); topoplot(datastSDchan, elocs(1:64)); c=colorbar; ylabel(c, 'SD');
                    %shows how noisy the channels are!
                    
                    
                    % concatenate all these types to increase quality of cov matrix.
                    if id==1 && windch==1%start here
                        if ndims(datast)<3
                        dataOUT(1,:,:)=datast;
                        else
                            dataOUT=datast;
                        end
                        
                    else
                        % append
                         if ndims(datast)<3
                        tmp(1,:,:)=datast;
                        dataOUT=cat(1,dataOUT,tmp);
                        
                        else
                            dataOUT=cat(1,dataOUT,datast);
                         end                        
                        
                    end
                    
                    end
                    
                end % concatenating all data sets.
                
                % We want a single component for this type of config, per ppant.
                % trials as last dimension
                %%
                ressd=permute(dataOUT,[2 3 1]);
                
                 
                ressEVEC_byHz = zeros(length(peakfreqsare), length(usechans));
                ressEVEC_byHz_MAPS=zeros(length(peakfreqsare),length(usechans));
                clf
                for ifreq=1:length(peakfreqsare)
                    
                    
                    
                    [ressEVEC, ~, ~, ressMAPS, hz]= constrRESSacrosstrials_willsData(ressd, peakfreqsare, usechans, neighbour, window,srate, ifreq);
                    
                    %
                    
                    % now we can multiple each type of data by the appropriate
                    %                 % ressEVEC
                    %%
                    figure(ifol); subplot(2,5,ifreq)
                    peakfreq1=peakfreqsare(ifreq);
                    topoplot(ressMAPS, elocs(usechans));
                    title(num2str(num2str(peakfreq1)));
                    %%
                    ressEVEC_byHz(ifreq,:) = ressEVEC;
                    ressEVEC_byHz_MAPS(ifreq,:)= ressMAPS;
                    
                    
                    
                end
        %%
            savename='RESSfilterdata';
        
        
        save(savename, 'ressEVEC_byHz', 'ressEVEC_byHz_MAPS', 'usechans', 'neighbour', 'windowsmall', 'peakfreqsare')
       %%
        cd([basefol filesep 'Figures' filesep 'RESS filters per participant'])
        suptitle(['Participant ' num2str(ifol) ])
       %%
        print('-dpng', ['RESS filters for  ppant ' num2str(ifol)]);
        disp(['fin for ppant' num2str(ifol)])
    end
end



