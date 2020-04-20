% clear all

addpath('/Users/MattDavidson/Desktop/SSVEP-Wills')
cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA')


basefol=pwd;
clearvars -except basefol allppants
dbstop if error

cd('EEG');
pdirs = dir([pwd filesep '*EEG']);


%% load data to determine physical catch timing.

allppants=[1,2,4,6,7,9:19]; 

window=[-3 3];
srate=250;

epochdur = sum(abs(window))*srate;
onsetc = ceil(epochdur)/2;
% peakfreqsare=[20,40]; %hz
%timing
tt = 0:1/srate:60;

% which frequencies to analyze?/apply?
peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; % don't change!  
%so that correct RESS filter gets applied.
%%

  
    for ifol = allppants
        
        
        %%
           cd(basefol)
        cd('EEG')
        cd(pdirs(ifol).name)
        
        
        %% load the relevant BEH/EEGdata.
        load('ppant_Catch_Epoched');
        
        
        load('RESSfilterdata');
        
        load('TrialIndicesbyCatchLocationandNum');
        %%
        
        %%%%%% POOL all conditions to increase the quality of the
        %%%%%% covariance matrix.
        dataIN=[];
        
        %%
        
        for id=1:5
            
            switch id
                case 1
                    dataIN=ppant_SNREEG_catchramponset;
                    
                    
                    % note the number of freqs,  then trials per location relevant. 
                    ress_catchonsetTGs=zeros(3,48,size(dataIN,3));                    
                    ress_catchonsetBGs=zeros(3,48,size(dataIN,3));
                    ress_catchonsetIMs=zeros(3,48,size(dataIN,3));
                case 2
                    dataIN=ppant_SNREEG_catchrampoffset;
                    ress_catchoffsetTGs=zeros(3,48,size(dataIN,3));
                    ress_catchoffsetBGs=zeros(3,48,size(dataIN,3));
                    ress_catchoffsetIMs=zeros(3,48,size(dataIN,3));
                    
                case 3
                    dataIN=ppant_SNREEG_disapBPwithincatch;
                    ress_BPcatchonsetTGs=zeros(3,48,size(dataIN,3));
                    ress_BPcatchonsetBGs=zeros(3,48,size(dataIN,3));
                    ress_BPcatchonsetIMs=zeros(3,48,size(dataIN,3));
                case 4
                    dataIN=ppant_SNREEG_reapBPaftercatch;
                    ress_BPcatchoffsetTGs=zeros(3,48,size(dataIN,3));
                    ress_BPcatchoffsetBGs=zeros(3,48,size(dataIN,3));
                    ress_BPcatchoffsetIMs=zeros(3,48,size(dataIN,3));
                case 5
                    dataIN=ppant_SNREEG_invisiblecatchonset;
                    ress_invisiblecatchonsetTGs=zeros(3,48,size(dataIN,3));
                    ress_invisiblecatchonsetBGs=zeros(3,48,size(dataIN,3));
                    ress_invisiblecatchonsetIMs=zeros(3,48,size(dataIN,3));
                                        
                    
            end
            %since all catch, no durs, 
            durscheck=[];
            searchtrials=1:24;
            
            for ifreq=1:length(peakfreqsare)
               
                    
                     %reduce size?
                    datast= dataIN;
                    
                    %now we have the correct trials, get the appropriate
                    %filter                    
                    evecs = squeeze(ressEVEC_byHz(ifreq,:)); %
                    
                    ress_ts1=zeros(size(datast,1), size(datast,3));
                    for ti=1:size(datast,1)
                        ress_ts1(ti,:) = evecs*squeeze(datast(ti,:,:));
                    end
                    
                    
                    %%
                    
                    %Save the ress data per type.

                       dimplacer=[1,1,2,2,3,3,1,2,3];
                       thisfreq=dimplacer(ifreq);
                    switch id
                        case 1
                            if ifreq<7  %TG or BG
                                if mod(ifreq,2)~=0 % TG
                                    ress_catchonsetTGs(thisfreq,:,:)= ress_ts1;
                                else %BG
                                    ress_catchonsetBGs(thisfreq,:,:)= ress_ts1;
                                end 
                                   
                            else%IMs in that case.
                                ress_catchonsetIMs(thisfreq,:,:)= ress_ts1;

                            end
                        case 2
                            
                             if ifreq<7  %TG or BG
                                if mod(ifreq,2)~=0 % TG
                                    ress_catchoffsetTGs(thisfreq,:,:)= ress_ts1;
                                else %BG
                                    ress_catchoffsetBGs(thisfreq,:,:)= ress_ts1;
                                end 
                                   
                            else%IMs in that case.
                                ress_catchoffsetIMs(thisfreq,:,:)= ress_ts1;

                            end
                            
                        case 3
                             if ifreq<7  %TG or BG
                                if mod(ifreq,2)~=0 % TG
                                    ress_BPcatchonsetTGs(thisfreq,:,:)= ress_ts1;
                                else %BG
                                    ress_BPcatchonsetBGs(thisfreq,:,:)= ress_ts1;
                                end 
                                   
                            else%IMs in that case.
                                ress_BPcatchonsetIMs(thisfreq,:,:)= ress_ts1;

                            end
                            
                        case 4
                            if ifreq<7  %TG or BG
                                if mod(ifreq,2)~=0 % TG
                                    ress_BPcatchoffsetTGs(thisfreq,:,:)= ress_ts1;
                                else %BG
                                    ress_BPcatchoffsetBGs(thisfreq,:,:)= ress_ts1;
                                end 
                                   
                            else%IMs in that case.
                                ress_BPcatchoffsetIMs(thisfreq,:,:)= ress_ts1;

                            end
                            
                            
                        case 5
                            if ifreq<7  %TG or BG
                                if mod(ifreq,2)~=0 % TG
                                    ress_invisiblecatchonsetTGs(thisfreq,:,:)= ress_ts1;
                                else %BG
                                    ress_invisiblecatchonsetBGs(thisfreq,:,:)= ress_ts1;
                                end 
                                   
                            else%IMs in that case.
                                ress_invisiblecatchonsetIMs(thisfreq,:,:)= ress_ts1;

                            end
                            
                    end
                
            end
            clearvars dataIN
            
        end
        
            savename='ppant_Catch_Epoched_RESS';

        save(savename, 'ress_catchonsetTGs', 'ress_catchonsetBGs','ress_catchonsetIMs',...
            'ress_catchoffsetTGs', 'ress_catchoffsetBGs','ress_catchoffsetIMs',...
            'ress_BPcatchonsetTGs', 'ress_BPcatchonsetBGs','ress_BPcatchonsetIMs',...
            'ress_BPcatchoffsetTGs', 'ress_BPcatchoffsetBGs','ress_BPcatchoffsetIMs',...
            'ress_invisiblecatchonsetTGs','ress_invisiblecatchonsetBGs','ress_invisiblecatchonsetIMs',...
        'peakfreqsare')
        
    end
    
    

