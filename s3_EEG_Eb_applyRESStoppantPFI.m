
clear all

addpath('/Users/MattDavidson/Desktop/SSVEP-Wills')
cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA')


basefol=pwd;
clearvars -except basefol allppants
dbstop if error

%%%%%% NOTE that at this stage, we sort to include ONLY relevant
%%%%%% disappearances and reappearances, where the RESS filter (HzxLoc) matches
%%%%%% the disap type (HzxLoc).

 
cd('EEG');
pdirs = dir([pwd filesep '*EEG']);

%% load data to determine physical catch timing.


allppants=[1,2,4,6,7,9:19]; %

% SET UP params.
window=[-3 3];
srate=250;
epochdur = sum(abs(window))*250;
onsetc = ceil(epochdur)/2;

%timing
tt = 0:1/250:60;


% which frequencies to analyze?/apply?
peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; % don't change!        

    for ifol =allppants
        
        
        %%
        cd(basefol)
        cd('EEG')
        cd(pdirs(ifol).name)
        
        %% load the relevant PFI data.
        load('ppant_PFI_Epoched');
        
        load('RESSfilterdata')
        
%         load('TrialIndicesbyLocationandHz.mat')

        dataIN=[];
        
        %%
        for id=1:8% Use all epochs so as not to bias condition comparisons.
            
            trialsbyHz=[];
            switch id
                case 1
                    dataIN=ppant_SNREEG_PFI_0_1;                    
                    BPstokeep= BPs0_1;
                case 2
                    dataIN=ppant_SNREEG_PFI_1_0;                    
                    BPstokeep= BPs1_0;
                    
                case 3
                    dataIN=ppant_SNREEG_PFI_1_2;                    
                    BPstokeep= BPs1_2;
                case 4
                    dataIN=ppant_SNREEG_PFI_2_1;                    
                    BPstokeep= BPs2_1;
                case 5
                    dataIN=ppant_SNREEG_PFI_2_3;                    
                    BPstokeep= BPs2_3;
                    
                case 6                    
                    dataIN=ppant_SNREEG_PFI_3_2;
                    BPstokeep= BPs3_2;
                case 7
                    dataIN=ppant_SNREEG_PFI_3_4;
                    BPstokeep= BPs3_4;
                case 8
                    
                    dataIN=ppant_SNREEG_PFI_4_3;
                    BPstokeep= BPs4_3;
                    
                    
            end
            
            for ifreq=1:length(peakfreqsare)
                usehz=peakfreqsare(ifreq);
                
                
                trialtypesperTargPresent = 1:size(dataIN,1);
                
               
                
                    datast=dataIN;
                    
                    % remove bad trials.
                    % check for bad trials (noisy)
                                        %std per trial(average over electrodes)
                                        tmp=[];
                                        if ndims(datast)<3
                                            tmp(1,:,:)= datast;
                                            datast=tmp;
                                        end
                                        datastSD = nanstd(squeeze(nanmean(datast,2)),0,2);
                    
                                        %remove those with 2.5*std from mean.
                                        trialSD=nanstd(datastSD);
                                        mSD=nanmean(datastSD);
                                        keeptrials=1:size(datast,1);
                    
                                        %remove the trials with excessive noise.
                                        badtrials=find(datastSD>mSD+2.5*trialSD)';
                    
                                        % also skip the trials which were only transient button
                                        % presses. (less than one second).
%                                         shorttrials = find(durscheck<60);
%                                         badtrials = [badtrials, shorttrials];
                    
                                        % remove these from consideration.
                                        keeptrials(badtrials)=[];
                                        datast=datast(keeptrials,:,:);
                    
                                        BPstosave = BPstokeep(keeptrials,:);
                    
                    %now we have the correct trials, get the appropriate
                    %filter
                    evecs = squeeze(ressEVEC_byHz(ifreq,:)); %
                    
                    
                    
                    ress_ts1=zeros(size(datast,1), size(datast,3));
                    for ti=1:size(datast,1)
                        ress_ts1(ti,:) = evecs*squeeze(datast(ti,:,:));
                    end
                    
                
                %Save the ress data per type.
                %also keep BP data which we will need for plotting.
                switch id
                    case 1
                        ress_PFI_0_1_Hz(ifreq).ressTS= ress_ts1;
                        ress_PFI_0_1_Hz(ifreq).BPs = BPstosave;
                        
                        
                        
                    case 2
                        
                        ress_PFI_1_0_Hz(ifreq).ressTS= ress_ts1;
                        ress_PFI_1_0_Hz(ifreq).BPs = BPstosave;
                        
                        
                        
                    case 3

                        ress_PFI_1_2_Hz(ifreq).ressTS= ress_ts1;
                        ress_PFI_1_2_Hz(ifreq).BPs = BPstosave;
                        
                        
                        
                        
                    case 4
                        
                        ress_PFI_2_1_Hz(ifreq).ressTS= ress_ts1;
                        ress_PFI_2_1_Hz(ifreq).BPs = BPstosave;
                        
                        
                    case 5
                        
                        ress_PFI_2_3_Hz(ifreq).ressTS= ress_ts1;
                        ress_PFI_2_3_Hz(ifreq).BPs = BPstosave;
                        
                        
                    case 6
                        
                        ress_PFI_3_2_Hz(ifreq).ressTS= ress_ts1;
                        ress_PFI_3_2_Hz(ifreq).BPs = BPstosave;
                        
                        
                    case 7
                        ress_PFI_3_4_Hz(ifreq).ressTS= ress_ts1;
                        ress_PFI_3_4_Hz(ifreq).BPs = BPstosave;
                    case 8
                        ress_PFI_4_3_Hz(ifreq).ressTS= ress_ts1;
                        ress_PFI_4_3_Hz(ifreq).BPs = BPstosave;
                        
                        
                end
            end
            clearvars dataIN 
            
                end
        
            savename='ppant_PFI_Epoched_RESS';
        
            
        save(savename, 'ress_PFI_0_1_Hz', 'ress_PFI_1_0_Hz',...
            'ress_PFI_1_2_Hz','ress_PFI_2_1_Hz', 'ress_PFI_2_3_Hz', 'ress_PFI_3_2_Hz','ress_PFI_4_3_Hz','ress_PFI_3_4_Hz')
        display([' finished ppant ' num2str(ifol)])
    end
    clearvars ress_*
