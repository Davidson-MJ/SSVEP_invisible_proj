% Calculate new participant PFI data, outside of catch events.
% formatted for new experiment, based on previous PFI-SSSVEP paper code


%%%%%%%%%%
%%%%%%%%
%%%%%%%%
%MD July 2018
%%%%%%%%
%%%%%%%%
%Update 2020.
clear all;  clc;

cd('/Users/mdavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA/Behaviour')
basedir=pwd;
dbstop if error
%%

%load BEH data
load('MD_AllBP_perppant.mat');
%load Catch Performance structure
load('Catchperformance.mat')
%work with retained participants.
allppants=[1,2,4,6,7,9:19]; %

fontsize=25;

job.removeCatchperiodsfromBPtrace=0;

%Sort PFI by number of targets invisible.
%participant level:
job.calcPFIdataperNum=0; %resaves into PFI_only with new fields in structure..
job.concatPFIacrossPpants_num=0;

job.createShufflePFIdata_pernum=0; % now with nshuff=1000
job.calcPFIdataperNum_shuffled=0;
job.concatPFIacrossPpants_num_shuffled=0;



%%%%%% plotting

% plot these together:
job.plotBehaviouraldata_num_with_shuffled=1;
job.compareSlopesofShuffledvsObserveddata=1;


nppants=19;
excludeTransientlength = 0; %minimum frames for counted PFI. (Fs=60),

%%
if job.removeCatchperiodsfromBPtrace==1
    %%
    PFI_only=[];
    for ippant=1:nppants %do all ppants, just analyze / plot retained, in case we change criteria at some point.
        
        BPtmp= ppantTrialDatawithDetails(ippant).AllBPData;
        catchremovedallBPtmp=BPtmp;
        for itrial=1:48
            %some trials had catch start / end ,yet no actual target removal.
            if ppantTrialDatawithDetails(ippant).TrialDetails(itrial).TotalCatchTargetsRemoved > 0
                %new analysis, replacing the catch target with NaNs, to remove from averaging.
                
                
                catchstart= ppantTrialDatawithDetails(ippant).TrialDetails(itrial).Catchstart_frames;
                catchend=ppantTrialDatawithDetails(ippant).TrialDetails(itrial).Catchend_frames;
                
                
                
                for iloc=1:4
                    switch iloc
                        case 1 %check TL
                            checkloc=ppantTrialDatawithDetails(ippant).TrialDetails(itrial).TL_Catch;
                        case 2
                            checkloc=ppantTrialDatawithDetails(ippant).TrialDetails(itrial).TR_Catch;
                        case 3
                            checkloc=ppantTrialDatawithDetails(ippant).TrialDetails(itrial).BL_Catch;
                        case 4
                            checkloc=ppantTrialDatawithDetails(ippant).TrialDetails(itrial).BR_Catch;
                    end
                    
                    if checkloc==1 %if this BP location contained catch event
                        %convert section of that BP trace to nan
                        nantrain = nan(1,(catchend-catchstart)+1);
                        tmptrace=squeeze(BPtmp(itrial,iloc,:))';
                        tmptrace(1,catchstart:catchend)=nantrain;
                        catchremovedallBPtmp(itrial,iloc, :)=tmptrace;
                    end
                end
                
                catchdur= ppantTrialDatawithDetails(ippant).TrialDetails(itrial).totalCatchdur;
                
                PFI_only(ippant).Trial(itrial).allBPs= squeeze(catchremovedallBPtmp(itrial,:,:));
                PFI_only(ippant).Trial(itrial).resultTriallength= length(BPtmp) - catchdur;
                
            else
                
                PFI_only(ippant).Trial(itrial).allBPs= squeeze(BPtmp(itrial,:,:));
                PFI_only(ippant).Trial(itrial).resultTriallength= length(BPtmp);
            end
        end
    end
    save('PFI_data', 'PFI_only')
end
%%

%%
if job.calcPFIdataperNum==1
    %%
    load('PFI_data.mat')
    %%
    
    
    for ippant = 1:nppants
        
        
        for itrial=1:48
            allBPtmp= PFI_only(ippant).Trial(itrial).allBPs;
            accumBP = squeeze(nansum(allBPtmp,1)); %combine each locations BP.
            accumBP(1,1)=0;
            dur0Disap=0;
            dur1Disap=0; %total time PFI
            dur2Disap=0;
            dur3Disap=0; %3
            dur4Disap=0;
            %tally disappearances.
            
            ZeroTargetDisap=0;
            OneTargetDisap=0;
            TwoTargetDisap=0;
            ThreeTargetDisap=0;
            FourTargetDisap=0;
            
            framestamp0_begin=[];
            framestamp0_end=[];
            framestamp1_begin=[];
            framestamp1_end=[];
            framestamp2_begin=[];
            framestamp2_end=[];
            framestamp3_begin=[];
            framestamp3_end=[];
            framestamp4_begin=[];
            framestamp4_end=[];
            
            
            
            rowis = (itrial+1) + 48*(ippant-1);
            
            %6th column is logical index for trial reject
            if strcmp(catchStruct{rowis,6},'0');
                
                for timeDur=2:(length(accumBP)-1)
                    
                    if accumBP(1,timeDur)==1 %if single button is pressed.
                        
                        
                        if accumBP(1,timeDur-1)~=1 %single button onset
                            if framestamp1_begin>0
                                framestamp1_begin= [framestamp1_begin, timeDur];
                            else
                                framestamp1_begin= timeDur;
                            end
                            OneTargetDisap=OneTargetDisap+1;
                            
                        end
                        
                        
                        if accumBP(1,timeDur+1)~=1 || timeDur==length(accumBP)-1
                            %singlebutton release
                            
                            if framestamp1_end>0
                                framestamp1_end= [framestamp1_end, timeDur+1];
                            else
                                framestamp1_end= timeDur+1;
                            end
                        end
                        
                        
                    elseif accumBP(1,timeDur)==2
                        
                        
                        
                        if accumBP(1,timeDur-1)~=2 %double button onset
                            if framestamp2_begin>0
                                framestamp2_begin= [framestamp2_begin, timeDur];
                            else
                                framestamp2_begin= timeDur;
                            end
                            TwoTargetDisap=TwoTargetDisap+1;
                            
                        end
                        
                        if accumBP(1,timeDur+1)~=2 || timeDur==length(accumBP)-1
                            %singlebutton offset
                            
                            if framestamp2_end>0
                                framestamp2_end= [framestamp2_end, timeDur+1];
                            else
                                framestamp2_end= timeDur+1;
                            end
                        end
                        
                        
                    elseif  accumBP(1,timeDur)==3
                        
                        
                        
                        if accumBP(1,timeDur-1)~=3 %double button onset
                            if framestamp3_begin>0
                                framestamp3_begin= [framestamp3_begin, timeDur];
                            else
                                framestamp3_begin= timeDur;
                            end
                            
                            ThreeTargetDisap=ThreeTargetDisap+1;
                            
                        end
                        
                        if accumBP(1,timeDur+1)~=3 || timeDur==length(accumBP)-1
                            %singlebutton offset
                            
                            if framestamp3_end>0
                                framestamp3_end= [framestamp3_end, timeDur+1];
                            else
                                framestamp3_end= timeDur+1;
                            end
                        end
                        
                        %%%%%%%%%% 4 targets disappeared.
                    elseif  accumBP(1,timeDur)==4
                        
                        
                        
                        if accumBP(1,timeDur-1)~=4
                            if framestamp4_begin>0
                                framestamp4_begin= [framestamp4_begin, timeDur];
                            else
                                framestamp4_begin= timeDur;
                            end
                            
                            FourTargetDisap=FourTargetDisap+1;
                            
                        end
                        
                        if accumBP(1,timeDur+1)~=4 || timeDur==length(accumBP)-1
                            %singlebutton offset
                            
                            if framestamp4_end>0
                                framestamp4_end= [framestamp4_end, timeDur+1];
                            else
                                framestamp4_end= timeDur+1;
                            end
                        end
                        
                    elseif accumBP(1,timeDur)==0 % now included for new figs.
                        if accumBP(1,timeDur-1) ~=0 %start of a disappearance.
                            if framestamp0_begin>0
                                framestamp0_begin = [framestamp0_begin, timeDur];
                            else
                                framestamp0_begin = timeDur;
                            end
                            
                            ZeroTargetDisap=ZeroTargetDisap+1;
                            
                            
                        end
                        
                        if accumBP(1,timeDur+1)~=0 || timeDur==length(accumBP)-1 %accounting for end of trial.
                            %singlebutton offset
                            
                            if framestamp0_end>0
                                framestamp0_end= [framestamp0_end, timeDur+1];
                            else
                                framestamp0_end= timeDur+1;
                            end
                            
                        end
                        
                    end
                end
                %in case a button was held until the end of a trial, we
                %need to remove these last disapperances
                
                
                %adjust for trial beginning.
                if accumBP(1,1)==0
                    framestamp0_begin = [1 framestamp0_begin];
                    ZeroTargetDisap = length(framestamp0_begin);
                end
                
                if length(framestamp0_begin) > length(framestamp0_end)
                    %adjust length
                    framestamp0_begin= framestamp0_begin(1:length(framestamp0_end));
                    
                end
                
                
                if length(framestamp1_begin) > length(framestamp1_end)
                    %adjust length
                    framestamp1_begin= framestamp1_begin(1:length(framestamp1_end));
                    OneTargetDisap = length(framestamp1_begin);
                end
                %same for 2 and three target instances.
                if  length(framestamp2_begin) > length(framestamp2_end)
                    
                    framestamp2_begin= framestamp2_begin(1:length(framestamp2_end));
                    TwoTargetDisap = length(framestamp2_begin);
                    
                end
                
                if  length(framestamp3_begin) > length(framestamp3_end)
                    
                    framestamp3_begin= framestamp3_begin(1:length(framestamp3_end));
                    ThreeTargetDisap = length(framestamp3_begin);
                end
                
                if  length(framestamp4_begin) > length(framestamp4_end)
                    
                    framestamp4_begin= framestamp4_begin(1:length(framestamp4_end));
                    FourTargetDisap= length(framestamp4_begin);
                end
                
                
                
                
                alldurs0 = framestamp0_end - framestamp0_begin;
                alldurs1 = framestamp1_end - framestamp1_begin;
                alldurs2 = framestamp2_end - framestamp2_begin;
                alldurs3= framestamp3_end - framestamp3_begin;
                alldurs4= framestamp4_end - framestamp4_begin;
                
                
                
                %for each type, remove the 'disappearances' which were only
                %transients.
                for in=1:5
                    switch in
                        case 1
                            ddurs=alldurs0;
                            nPFI=ZeroTargetDisap;
                            frame_end=framestamp0_end;
                            frame_begin=framestamp0_begin;
                        case 2
                            ddurs=alldurs1;
                            nPFI=OneTargetDisap;
                            frame_end=framestamp1_end;
                            frame_begin=framestamp1_begin;
                        case 3
                            ddurs=alldurs2;
                            nPFI=TwoTargetDisap;
                            frame_end=framestamp2_end;
                            frame_begin=framestamp2_begin;
                        case 4
                            ddurs=alldurs3;
                            nPFI=ThreeTargetDisap;
                            frame_end=framestamp3_end;
                            frame_begin=framestamp3_begin;
                        case 5
                            ddurs=alldurs4;
                            nPFI=FourTargetDisap;
                            frame_end=framestamp4_end;
                            frame_begin=framestamp4_begin;
                    end
                    
                    keep= find(ddurs>excludeTransientlength);
                    
                    %now adjust outgoings to only keep the PFI that were
                    %not transients.
                    durDisap = sum(ddurs(keep));
                    nPFI=length(keep);
                    alldurs = ddurs(keep);
                    
                    frame_begs = frame_begin(keep);
                    frame_ends = frame_end(keep);
                    %store
                    switch in
                        case 1
                            %mean duration per individual disap (in seconds).
                            PFI_only(ippant).Trial(itrial).meanPFIduration_perdisap_0target = durDisap/nPFI/60;
                            %  total duration per trial (in seconds)
                            PFI_only(ippant).Trial(itrial).PFItotalduration_disap_0target = durDisap/60;
                            %frequency of PFI per trial (how many times PFI occurs)
                            PFI_only(ippant).Trial(itrial).PFI_num0target_unadjusted = nPFI;
                            %record time stamps (in frames) of these events.
                            PFI_only(ippant).Trial(itrial).PFI_disap_0target_framestart= frame_begs;
                            PFI_only(ippant).Trial(itrial).PFI_disap_0target_frameend= frame_ends;
                            PFI_only(ippant).Trial(itrial).PFI_disap_0target_durs= alldurs;
                            
                        case 2
                            PFI_only(ippant).Trial(itrial).meanPFIduration_perdisap_1target = durDisap/nPFI/60;
                            PFI_only(ippant).Trial(itrial).PFItotalduration_disap_1target = durDisap/60;
                            PFI_only(ippant).Trial(itrial).PFI_num1target_unadjusted = nPFI;
                            PFI_only(ippant).Trial(itrial).PFI_disap_1target_framestart= frame_begs;
                            PFI_only(ippant).Trial(itrial).PFI_disap_1target_frameend= frame_ends;
                            PFI_only(ippant).Trial(itrial).PFI_disap_1target_durs= alldurs;
                        case 3
                            PFI_only(ippant).Trial(itrial).meanPFIduration_perdisap_2target = durDisap/nPFI/60;
                            PFI_only(ippant).Trial(itrial).PFItotalduration_disap_2target = durDisap/60;
                            PFI_only(ippant).Trial(itrial).PFI_num2target_unadjusted = nPFI;
                            PFI_only(ippant).Trial(itrial).PFI_disap_2target_framestart= frame_begs;
                            PFI_only(ippant).Trial(itrial).PFI_disap_2target_frameend= frame_ends;
                            PFI_only(ippant).Trial(itrial).PFI_disap_2target_durs= alldurs;
                        case 4
                            PFI_only(ippant).Trial(itrial).meanPFIduration_perdisap_3target = durDisap/nPFI/60;
                            PFI_only(ippant).Trial(itrial).PFItotalduration_disap_3target = durDisap/60;
                            PFI_only(ippant).Trial(itrial).PFI_num3target_unadjusted = nPFI;
                            PFI_only(ippant).Trial(itrial).PFI_disap_3target_framestart= frame_begs;
                            PFI_only(ippant).Trial(itrial).PFI_disap_3target_frameend= frame_ends;
                            PFI_only(ippant).Trial(itrial).PFI_disap_3target_durs= alldurs;
                        case 5
                            PFI_only(ippant).Trial(itrial).meanPFIduration_perdisap_4target = durDisap/nPFI/60;
                            PFI_only(ippant).Trial(itrial).PFItotalduration_disap_4target = durDisap/60;
                            PFI_only(ippant).Trial(itrial).PFI_num4target_unadjusted = nPFI;
                            PFI_only(ippant).Trial(itrial).PFI_disap_4target_framestart= frame_begs;
                            PFI_only(ippant).Trial(itrial).PFI_disap_4target_frameend= frame_ends;
                            PFI_only(ippant).Trial(itrial).PFI_disap_4target_durs= alldurs;
                            
                    end
                end
                
                %collect
                PFI_only(ippant).Trial(itrial).Goodtrial=1;
            else
                PFI_only(ippant).Trial(itrial).Goodtrial=0;
            end
            
            
        end
        
    end
    
    save('PFI_data', 'PFI_only')
end
%%

%%
if job.concatPFIacrossPpants_num==1
    load('PFI_data')
    Freq_NumPFI_acrossTrials = nan(nppants,48,5);
    mDurperNumPFI_acrossTrials=nan(nppants,48,5);
    totalDurperNumPFI_acrossTrials=nan(nppants,48,5);
    
    for ippant = 1:nppants
        
        for itrial=1:48
            
            if PFI_only(ippant).Trial(itrial).Goodtrial==1 %not rejected
                usedata=PFI_only(ippant).Trial(itrial);
                %zero targets first
                Freq_NumPFI_acrossTrials(ippant,itrial,1) = usedata.PFI_num0target_unadjusted;
                mDurperNumPFI_acrossTrials(ippant,itrial,1)=usedata.meanPFIduration_perdisap_0target;
                totalDurperNumPFI_acrossTrials(ippant,itrial,1)=usedata.PFItotalduration_disap_0target;
                
                Freq_NumPFI_acrossTrials(ippant,itrial,2) = usedata.PFI_num1target_unadjusted;
                mDurperNumPFI_acrossTrials(ippant,itrial,2)=usedata.meanPFIduration_perdisap_1target;
                totalDurperNumPFI_acrossTrials(ippant,itrial,2)=usedata.PFItotalduration_disap_1target;
                %double
                
                Freq_NumPFI_acrossTrials(ippant,itrial,3) = usedata.PFI_num2target_unadjusted;
                mDurperNumPFI_acrossTrials(ippant,itrial,3)=usedata.meanPFIduration_perdisap_2target;
                totalDurperNumPFI_acrossTrials(ippant,itrial,3)=usedata.PFItotalduration_disap_2target;
                
                %'3
                
                Freq_NumPFI_acrossTrials(ippant,itrial,4) = usedata.PFI_num3target_unadjusted;
                mDurperNumPFI_acrossTrials(ippant,itrial,4)=usedata.meanPFIduration_perdisap_3target;
                totalDurperNumPFI_acrossTrials(ippant,itrial,4)=usedata.PFItotalduration_disap_3target;
                
                %4
                Freq_NumPFI_acrossTrials(ippant,itrial,5) = usedata.PFI_num4target_unadjusted;
                mDurperNumPFI_acrossTrials(ippant,itrial,5)=usedata.meanPFIduration_perdisap_4target;
                totalDurperNumPFI_acrossTrials(ippant,itrial,5)=usedata.PFItotalduration_disap_4target;
                
            end
            
        end
    end
    
    save('PFI_data_concat', 'Freq_NumPFI_acrossTrials', 'mDurperNumPFI_acrossTrials', 'totalDurperNumPFI_acrossTrials')
    
    
end
%%

%% based on previous code (reshuffle_analysis)
if job.createShufflePFIdata_pernum==1
    load('MD_AllBP_perppant.mat')
    nshuff = 1000;
    allRandomAllPP =zeros(length(allppants),nshuff,5,3600);
    
    
    for ippant=1:length(allppants)
        goodPP=allppants(ippant);
        
        cd(basedir)
        
        
        BPdata=ppantTrialDatawithDetails(goodPP).AllBPData;
        
        allRandom = zeros(nshuff,5, 3600);
        for countRand=1:nshuff
            %index of four random trials.
            randTrial=randi([1 48],1,4);
            
            
            %take a location on screen from separate trials, and concatenate.
            plotTable(1,:)=BPdata(randTrial(1),1,:);
            plotTable(2,:)=BPdata(randTrial(2),2,:);
            plotTable(3,:)=BPdata(randTrial(3),3,:);
            plotTable(4,:)=BPdata(randTrial(4),4,:);
            %accumulative sum
            plotTable(5,:)=sum(plotTable(1:4,:),1);
            
            allRandom(countRand,:,:)= plotTable;
            
        end
        
        allRandomAllPP(ippant,:,:,:)=allRandom;
        
        
    end
    
    try save('ShuffledData', 'allRandomAllPP')
    catch
        save('ShuffledData', 'allRandomAllPP', '-V7.3')
    end
    
    
    %%
end

if job.calcPFIdataperNum_shuffled==1 %resaves into ShuffledData
    
    cd(basedir)
    load('ShuffledData')
    for ippant = 1:size(allRandomAllPP,1)
        
        
        for itrial=1:1000
            
            accumBP = squeeze(allRandomAllPP(ippant,itrial,5,:))'; %combine each locations BP.
            dur0Disap=0;
            dur1Disap=0; %total time PFI
            dur2Disap=0;
            dur3Disap=0; %3
            dur4Disap=0; %4
            
            ZeroTargetDisap=0;
            
            OneTargetDisap=0;
            TwoTargetDisap=0;
            ThreeTargetDisap=0;
            FourTargetDisap=0;
            
            
            framestamp0_begin=[];
            framestamp0_end=[];
            
            framestamp1_begin=[];
            framestamp1_end=[];
            framestamp2_begin=[];
            framestamp2_end=[];
            framestamp3_begin=[];
            framestamp3_end=[];
            
            framestamp4_begin=[];
            framestamp4_end=[];
            
            
            
            accumBP(1)=0;
            for timeDur=2:(length(accumBP)-1)
                
                if accumBP(1,timeDur)==1
                    
                    
                    if accumBP(1,timeDur-1)~=1 %single button onset
                        if framestamp1_begin>0
                            framestamp1_begin= [framestamp1_begin, timeDur];
                        else
                            framestamp1_begin= timeDur;
                        end
                        OneTargetDisap=OneTargetDisap+1;
                        
                    end
                    
                    
                    if accumBP(1,timeDur+1)~=1 || timeDur==length(accumBP)-1
                        %singlebutton offset
                        
                        if framestamp1_end>0
                            framestamp1_end= [framestamp1_end, timeDur+1];
                        else
                            framestamp1_end= timeDur+1;
                        end
                    end
                    
                    
                elseif accumBP(1,timeDur)==2
                    
                    
                    
                    if accumBP(1,timeDur-1)~=2 %double button onset
                        if framestamp2_begin>0
                            framestamp2_begin= [framestamp2_begin, timeDur];
                        else
                            framestamp2_begin= timeDur;
                        end
                        TwoTargetDisap=TwoTargetDisap+1;
                        
                    end
                    
                    if accumBP(1,timeDur+1)~=2 || timeDur==length(accumBP)-1
                        %singlebutton offset
                        
                        if framestamp2_end>0
                            framestamp2_end= [framestamp2_end, timeDur+1];
                        else
                            framestamp2_end= timeDur+1;
                        end
                    end
                    
                    
                elseif  accumBP(1,timeDur)==3
                    
                    
                    
                    if accumBP(1,timeDur-1)~=3 %
                        if framestamp3_begin>0
                            framestamp3_begin= [framestamp3_begin, timeDur];
                        else
                            framestamp3_begin= timeDur;
                        end
                        ThreeTargetDisap=ThreeTargetDisap+1;
                        
                    end
                    
                    if accumBP(1,timeDur+1)~=3 || timeDur==length(accumBP)-1
                        %singlebutton offset
                        
                        if framestamp3_end>0
                            framestamp3_end= [framestamp3_end, timeDur+1];
                        else
                            framestamp3_end= timeDur+1;
                        end
                    end
                elseif  accumBP(1,timeDur)==4
                    
                    
                    
                    if accumBP(1,timeDur-1)~=4 %
                        if framestamp4_begin>0
                            framestamp4_begin= [framestamp4_begin, timeDur];
                        else
                            framestamp4_begin= timeDur;
                        end
                        FourTargetDisap=FourTargetDisap+1;
                        
                    end
                    
                    if accumBP(1,timeDur+1)~=4 || timeDur==length(accumBP)-1
                        %singlebutton offset
                        
                        if framestamp4_end>0
                            framestamp4_end= [framestamp4_end, timeDur+1];
                        else
                            framestamp4_end= timeDur+1;
                        end
                    end
                    
                elseif accumBP(1,timeDur)==0 % now included for new figs.
                    if accumBP(1,timeDur-1) ~=0 %start of a disappearance.
                        if framestamp0_begin>0
                            framestamp0_begin = [framestamp0_begin, timeDur];
                        else
                            framestamp0_begin = timeDur;
                        end
                        ZeroTargetDisap=ZeroTargetDisap+1;
                    end
                    
                    if accumBP(1,timeDur+1)~=0 || timeDur==length(accumBP)-1 %accounting for end of trial.
                        %singlebutton offset
                        
                        if framestamp0_end>0
                            framestamp0_end= [framestamp0_end, timeDur+1];
                        else
                            framestamp0_end= timeDur+1;
                        end
                        
                    end
                    
                    
                end
                
            end
            
            %in case a button was held until the end of a trial, we
            %need to remove these last disapperances
            
            %adjust for trial beginning.
            if accumBP(1,1)==0
                framestamp0_begin = [1 framestamp0_begin];
                
                ZeroTargetDisap = length(framestamp0_begin);
                
            end
            
            if length(framestamp0_begin) > length(framestamp0_end)
                %adjust length
                framestamp0_begin= framestamp0_begin(1:length(framestamp0_end));
                
            end
            
            if length(framestamp1_begin) > length(framestamp1_end)
                %adjust length
                framestamp1_begin= framestamp1_begin(1:length(framestamp1_end));
                OneTargetDisap = length(framestamp1_begin);
            end
            %same for 2 and three target instances.
            if  length(framestamp2_begin) > length(framestamp2_end)
                
                framestamp2_begin= framestamp2_begin(1:length(framestamp2_end));
                TwoTargetDisap = length(framestamp2_begin);
                
            end
            
            if  length(framestamp3_begin) > length(framestamp3_end)
                
                framestamp3_begin= framestamp3_begin(1:length(framestamp3_end));
                ThreeTargetDisap = length(framestamp3_begin);
            end
            
            if  length(framestamp4_begin) > length(framestamp4_end)
                
                framestamp4_begin= framestamp4_begin(1:length(framestamp4_end));
                FourTargetDisap= length(framestamp4_begin);
            end
            
            
            
            
            alldurs0 = framestamp0_end - framestamp0_begin;
            alldurs1 = framestamp1_end - framestamp1_begin;
            alldurs2 = framestamp2_end - framestamp2_begin;
            alldurs3= framestamp3_end - framestamp3_begin;
            alldurs4= framestamp4_end - framestamp4_begin;
            
            
            dur0Disap= sum(alldurs0);
            dur1Disap = sum(alldurs1);
            dur2Disap = sum(alldurs2);
            dur3Disap = sum(alldurs3);
            dur4Disap = sum(alldurs4);
            
            %for each type, remove the 'disappearances' which were only
            %transients.
            for in=1:5
                switch in
                    case 1
                        ddurs=alldurs0;
                        nPFI=ZeroTargetDisap;
                        frame_end=framestamp0_end;
                        frame_begin=framestamp0_begin;
                    case 2
                        ddurs=alldurs1;
                        nPFI=OneTargetDisap;
                        frame_end=framestamp1_end;
                        frame_begin=framestamp1_begin;
                    case 3
                        ddurs=alldurs2;
                        nPFI=TwoTargetDisap;
                        frame_end=framestamp2_end;
                        frame_begin=framestamp2_begin;
                    case 4
                        ddurs=alldurs3;
                        nPFI=ThreeTargetDisap;
                        frame_end=framestamp3_end;
                        frame_begin=framestamp3_begin;
                    case 5
                        ddurs=alldurs4;
                        nPFI=FourTargetDisap;
                        frame_end=framestamp4_end;
                        frame_begin=framestamp4_begin;
                end
                
                keep= find(ddurs>excludeTransientlength);
                
                %now adjust outgoings to only keep the PFI that were
                %not transients.
                durDisap = sum(ddurs(keep));
                nPFI=length(keep);
                alldurs = ddurs(keep);
                
                frame_begs = frame_begin(keep);
                frame_ends = frame_end(keep);
                %store
                switch in
                    case 1
                        %mean duration per individual disap (in seconds).
                        shuffledPFI(ippant).Trial(itrial).meanPFIduration_perdisap_0target = durDisap/nPFI/60;
                        %  total duration per trial (in seconds)
                        shuffledPFI(ippant).Trial(itrial).PFItotalduration_disap_0target = durDisap/60;
                        %frequency of PFI per trial (how many times PFI occurs)
                        shuffledPFI(ippant).Trial(itrial).PFI_num0target_unadjusted = nPFI;
                        %record time stamps (in frames) of these events.
                        shuffledPFI(ippant).Trial(itrial).PFI_disap_0target_framestart= frame_begs;
                        shuffledPFI(ippant).Trial(itrial).PFI_disap_0target_frameend= frame_ends;
                        shuffledPFI(ippant).Trial(itrial).PFI_disap_0target_durs= alldurs;
                        
                    case 2
                        shuffledPFI(ippant).Trial(itrial).meanPFIduration_perdisap_1target = durDisap/nPFI/60;
                        shuffledPFI(ippant).Trial(itrial).PFItotalduration_disap_1target = durDisap/60;
                        shuffledPFI(ippant).Trial(itrial).PFI_num1target_unadjusted = nPFI;
                        shuffledPFI(ippant).Trial(itrial).PFI_disap_1target_framestart= frame_begs;
                        shuffledPFI(ippant).Trial(itrial).PFI_disap_1target_frameend= frame_ends;
                        shuffledPFI(ippant).Trial(itrial).PFI_disap_1target_durs= alldurs;
                    case 3
                        shuffledPFI(ippant).Trial(itrial).meanPFIduration_perdisap_2target = durDisap/nPFI/60;
                        shuffledPFI(ippant).Trial(itrial).PFItotalduration_disap_2target = durDisap/60;
                        shuffledPFI(ippant).Trial(itrial).PFI_num2target_unadjusted = nPFI;
                        shuffledPFI(ippant).Trial(itrial).PFI_disap_2target_framestart= frame_begs;
                        shuffledPFI(ippant).Trial(itrial).PFI_disap_2target_frameend= frame_ends;
                        shuffledPFI(ippant).Trial(itrial).PFI_disap_2target_durs= alldurs;
                    case 4
                        shuffledPFI(ippant).Trial(itrial).meanPFIduration_perdisap_3target = durDisap/nPFI/60;
                        shuffledPFI(ippant).Trial(itrial).PFItotalduration_disap_3target = durDisap/60;
                        shuffledPFI(ippant).Trial(itrial).PFI_num3target_unadjusted = nPFI;
                        shuffledPFI(ippant).Trial(itrial).PFI_disap_3target_framestart= frame_begs;
                        shuffledPFI(ippant).Trial(itrial).PFI_disap_3target_frameend= frame_ends;
                        shuffledPFI(ippant).Trial(itrial).PFI_disap_3target_durs= alldurs;
                    case 5
                        shuffledPFI(ippant).Trial(itrial).meanPFIduration_perdisap_4target = durDisap/nPFI/60;
                        shuffledPFI(ippant).Trial(itrial).PFItotalduration_disap_4target = durDisap/60;
                        shuffledPFI(ippant).Trial(itrial).PFI_num4target_unadjusted = nPFI;
                        shuffledPFI(ippant).Trial(itrial).PFI_disap_4target_framestart= frame_begs;
                        shuffledPFI(ippant).Trial(itrial).PFI_disap_4target_frameend= frame_ends;
                        shuffledPFI(ippant).Trial(itrial).PFI_disap_4target_durs= alldurs;
                        
                end
            end
            
            
            %collect
            shuffledPFI(ippant).Trial(itrial).Goodtrial=1;
            
            
            
        end
        
    end
    
    save('ShuffledData', 'shuffledPFI', '-append')
    
end
%%
if job.concatPFIacrossPpants_num_shuffled==1
    cd(basedir)
    load('ShuffledData')
    
    Freq_NumPFI_acrossTrials_shuffled = nan(length(allppants),nshuff,5);
    mDurperNumPFI_acrossTrials_shuffled=nan(length(allppants),nshuff,5);
    totalDurperNumPFI_acrossTrials_shuffled=nan(length(allppants),nshuff,5);
    
    for ippant = 1:length(shuffledPFI)
        
        for itrial=1:nshuff
            
            if shuffledPFI(ippant).Trial(itrial).Goodtrial==1 %not rejected
                usedata=shuffledPFI(ippant).Trial(itrial);
                
                %zero targets first
                Freq_NumPFI_acrossTrials_shuffled(ippant,itrial,1) = usedata.PFI_num0target_unadjusted;
                
                
                tmp=usedata.meanPFIduration_perdisap_0target;
                if isinf(tmp)
                    error('check')
                else
                    mDurperNumPFI_acrossTrials_shuffled(ippant,itrial,1)=usedata.meanPFIduration_perdisap_0target;
                end
                totalDurperNumPFI_acrossTrials_shuffled(ippant,itrial,1)=usedata.PFItotalduration_disap_0target;
                
                
                Freq_NumPFI_acrossTrials_shuffled(ippant,itrial,2) = usedata.PFI_num1target_unadjusted;
                mDurperNumPFI_acrossTrials_shuffled(ippant,itrial,2)=usedata.meanPFIduration_perdisap_1target;
                totalDurperNumPFI_acrossTrials_shuffled(ippant,itrial,2)=usedata.PFItotalduration_disap_1target;
                %double
                
                Freq_NumPFI_acrossTrials_shuffled(ippant,itrial,3) = usedata.PFI_num2target_unadjusted;
                mDurperNumPFI_acrossTrials_shuffled(ippant,itrial,3)=usedata.meanPFIduration_perdisap_2target;
                totalDurperNumPFI_acrossTrials_shuffled(ippant,itrial,3)=usedata.PFItotalduration_disap_2target;
                
                %3
                
                Freq_NumPFI_acrossTrials_shuffled(ippant,itrial,4) = usedata.PFI_num3target_unadjusted;
                mDurperNumPFI_acrossTrials_shuffled(ippant,itrial,4)=usedata.meanPFIduration_perdisap_3target;
                totalDurperNumPFI_acrossTrials_shuffled(ippant,itrial,4)=usedata.PFItotalduration_disap_3target;
                
                %4
                Freq_NumPFI_acrossTrials_shuffled(ippant,itrial,5) = usedata.PFI_num4target_unadjusted;
                mDurperNumPFI_acrossTrials_shuffled(ippant,itrial,5)=usedata.meanPFIduration_perdisap_4target;
                totalDurperNumPFI_acrossTrials_shuffled(ippant,itrial,5)=usedata.PFItotalduration_disap_4target;
                
            end
            
        end
        
        %         mean(mDurperNumPFI_acrossTrials_shuffled(ippant,:,:),2)
    end
    
    save('ShuffledData', 'Freq_NumPFI_acrossTrials_shuffled', 'mDurperNumPFI_acrossTrials_shuffled', 'totalDurperNumPFI_acrossTrials_shuffled', '-append')
    
end


%% Plot Behavioural output.


if job.plotBehaviouraldata_num_with_shuffled==1
    %slim version of the above.
    cd(basedir)
    load('PFI_data_concat')
    load('ShuffledData')
    cd ../
    cd('Figures')
    cd('Behavioural Bar Plots')
    %%
    
    
    plotallNUM=1; % change to plot 0-4 targets (5 possible responses).
    
    
    %compare shuffled to number absent
    p1= Freq_NumPFI_acrossTrials;
    p2= mDurperNumPFI_acrossTrials;
    p3= totalDurperNumPFI_acrossTrials;
    xlabelis = 'nPFI';
    if plotallNUM==1
        xticks = [{'0'}, {'1'} {'2'}, {'3'}, {'4'}];
    else
        xticks = [{'1'} {'2'}, {'3'}, {'4'}];
    end
    
    % allocate shuffled data.
    s1=Freq_NumPFI_acrossTrials_shuffled;
    s2= mDurperNumPFI_acrossTrials_shuffled;
    s3= totalDurperNumPFI_acrossTrials_shuffled;
    
    
    
    clf
    for iDV=1:3
        %take mean within, then across ppants.
        %within
        
        subplot(2,3,iDV)
        switch iDV
            case 1
                datatoplot=p1;
                shuffledtoplot=s1;
                yis='PFI per trial';
            case 2
                
                datatoplot=p2;
                shuffledtoplot=s2;
                yis='PFI duration [s]';
                
            case 3
                datatoplot=p3;
                shuffledtoplot=s3;
                yis= 'Total duration [s]';
                
        end
        %%
        %only look at good ppants.
        datanowreal=squeeze(nanmean(datatoplot(allppants,:,:),2));
        datanowshuffled=squeeze(nanmean(shuffledtoplot,2));
        
        %across
        mDataReal=squeeze(nanmean(datanowreal,1));
        mDataShuff=squeeze(nanmean(datanowshuffled,1));
        %rearrange for plots
        
        if plotallNUM==1
            mData = [mDataReal(1,1), mDataShuff(1,1); mDataReal(1,2), mDataShuff(1,2); mDataReal(1,3), mDataShuff(1,3);mDataReal(1,4), mDataShuff(1,4);mDataReal(1,5), mDataShuff(1,5)];
        else
            mData = [mDataReal(1,2), mDataShuff(1,2); mDataReal(1,3), mDataShuff(1,3); mDataReal(1,4), mDataShuff(1,4); mDataReal(1,5), mDataShuff(1,5)];
        end
        %standard err m
        %standard err m  %note different calculations for data and
        %shuffled.
        for id=1:2
            switch id
                case 1
                    x=datanowreal;
                case 2
                    x=shuffledtoplot;
                    
            end
            
            if id==1
                mXppant =squeeze( nanmean(x,2)); %mean across conditions we are comparing (within ppant ie. RMfactors)
                mXgroup = mean(nanmean(mXppant)); %mean overall (remove b/w sub differences
                
                %for each observation, subjtract the subj average, add
                %the group average.
                NEWdata = x - repmat(mXppant, 1,size(x,2)) + repmat(mXgroup, size(x));
                
                std_tmp = nanstd(NEWdata);
                %careful here, as we need to divide by the number of data
                %points
                stErr=zeros(size(std_tmp));
                for idv=1:length(std_tmp);
                    %numberobs
                    n_isnan = length(find(isnan(NEWdata(:,idv))));
                    %correct for stErr.
                    stErr(idv)= std_tmp(idv)/ sqrt(size(NEWdata,1)-n_isnan);
                end
                
            else
                
                %loop over all shuffles.
                allstE = zeros(size(x,2), size(x,3));
                for ishuff=1:size(x,2)
                    newx = squeeze(x(:,ishuff,:));
                    mXppant =squeeze( nanmean(newx,2)); %mean across conditions we are comparing (within ppant ie. RMfactors)
                    mXgroup = mean(nanmean(mXppant)); %mean overall (remove b/w sub differences
                    
                    %for each observation, subjtract the subj average, add
                    %the group average.
                    NEWdata = newx - repmat(mXppant, 1,size(newx,2)) + repmat(mXgroup, size(newx));
                    
                    std_tmp = nanstd(NEWdata);
                    %careful here, as we need to divide by the number of data
                    %points
                    stErr=zeros(size(std_tmp));
                    for idv=1:length(std_tmp);
                        %numberobs
                        n_isnan = length(find(isnan(NEWdata(:,idv))));
                        %correct for stErr.
                        stErr(idv)= std_tmp(idv)/ sqrt(size(NEWdata,1)-n_isnan);
                    end
                    %             stErr = nanstd(NEWdata)/sqrt(size(newx,1));
                    
                    
                    allstE(ishuff, :) = stErr;
                end
                
                %take mean across shuffles.
                stErr = nanmean(allstE,1);
            end
            
            switch id
                case 1
                    stErrReal = stErr;
                case 2
                    stErrShuffled= stErr;
            end
        end
        
        %             for id=1:2
        %                 switch id
        %                     case 1
        %                         x=datanowreal;
        %                     case 2
        %                         x=datanowshuffled;
        %
        %                 end
        %                 mXppant =squeeze( nanmean(x,2)); %mean across conditions we are comparing (within ppant ie. RMfactors)
        %             mXgroup = mean(nanmean(mXppant)); %mean overall (remove b/w sub differences
        %
        %             %for each observation, subjtract the subj average, add
        %             %the group average.
        %             NEWdata = x - repmat(mXppant, 1,size(x,2)) + repmat(mXgroup, size(x));
        %
        %             stErr = nanstd(NEWdata)/sqrt(size(x,1));
        %             switch id
        %                 case 1
        %                     stErrReal = stErr;
        %                 case 2
        %                     stErrShuffled= stErr;
        %             end
        %             end
        
        offset=.15;
        
        
        if plotallNUM==1
            stErr=[stErrReal(1), stErrShuffled(1);stErrReal(2), stErrShuffled(2);stErrReal(3), stErrShuffled(3) ;stErrReal(4), stErrShuffled(4);stErrReal(5), stErrShuffled(5) ];
            xis=[1-offset, 1+offset; 2-offset, 2+offset; 3-offset,3+offset; 4-offset,4+offset; 5-offset, 5+offset];
        else
            stErr=[stErrReal(2), stErrShuffled(2);stErrReal(3), stErrShuffled(3) ;stErrReal(4), stErrShuffled(4);stErrReal(5), stErrShuffled(5) ];
            xis=[1-offset, 1+offset; 2-offset, 2+offset; 3-offset,3+offset; 4-offset, 4+offset];
        end
        %
        %%
        %remove real for ASSC
        %             mData(:,1) = nan;
        bh=bar(mData);
        bh(1).FaceColor= ['b'];
        bh(2).FaceColor= [.5 .5 .5];
        %%
        hold on
        
        errorbar(xis,mData,stErr,'linestyle', ['none'], 'linewidth', 2 , 'Color', 'k')
        
        
        %%
        if iDV==3
            legend 'Observed' 'Shuffled'
        end
        
        axis('tight')
        switch iDV
            case 1
                ylabel([yis]);
                ylimsa=[0 10];
                set(gca, 'ytick', [0:2:ylimsa(2)])
            case 2
                ylabel({[yis ]});
                ylimsa=[0 10];
                set(gca, 'ytick', [0:1:ylimsa(2)])
                
            case 3
                ylabel({[yis ]});
                if plotallNUM==1
                    ylimsa=[0 40];
                else
                    ylimsa=[0 25];
                end
                set(gca, 'ytick', [0:5:ylimsa(2)])
                
                
        end
        axis tight
        ylim([ylimsa])
        
        
        set(gca, 'xticklabel', xticks)
        set(gca, 'fontsize', fontsize-2)
        
        
        xlabel(xlabelis)
        xlim([.5 5.5])
    end
    
    shg
    %%
    set(gcf, 'color', 'w')
    if plotallNUM==1
        printfilename = ['PFI data by ' xlabelis ' shuffled'];
    else
        printfilename = ['PFI data by ' xlabelis ' shuffled, no 0 PFI'];
    end
    %%
    print('-dpng', printfilename)
end


if job.compareSlopesofShuffledvsObserveddata==1
    
    
    cd(basedir)
    load('PFI_data_concat')
    load('ShuffledData')
    
    
    testnPFIinclZERO=1;
    
    p1= Freq_NumPFI_acrossTrials;
    p2= mDurperNumPFI_acrossTrials;
    p3= totalDurperNumPFI_acrossTrials;
    
    % allocate shuffled data.
    s1=Freq_NumPFI_acrossTrials_shuffled;
    s2= mDurperNumPFI_acrossTrials_shuffled;
    s3= totalDurperNumPFI_acrossTrials_shuffled;
    %
    
    pvalswere=[];
    figure(1);
    for iDV=1:3
        %take mean within, then across ppants.
        %within
        
        subplot(2,3,iDV)
        switch iDV
            case 1
                datatoplot=p1;
                shuffledtoplot=s1;
                yis='PFI per trial';
            case 2
                
                datatoplot=p2;
                shuffledtoplot=s2;
                yis='PFI duration [s]';
                
            case 3
                datatoplot=p3;
                shuffledtoplot=s3;
                yis= 'Total duration [s]';
                if testnPFIinclZERO==1
                    xticks = [{'0'},{'1'} {'2'}, {'3'}, {'4'}];
                end
        end
        if testnPFIinclZERO==1
            xticks = [{'0'},{'1'} {'2'}, {'3'}, {'4'}];
            usedata=1:5;
            xrange=usedata;
        else
            usedata=2:5;
            xrange=1:4;
        end
        %%
        %only look at good ppants.
        pl=[];
        for iregres=1:2
            if iregres==1
                datanowreal=squeeze(nanmean(shuffledtoplot(:,:,usedata),2));
                Bcolor=[.5 .5 .5];
                dcolor= [.5 .5 .5];
            else
                
                datanowreal=squeeze(nanmean(datatoplot(allppants,:,usedata),2));
                Bcolor='b';
                dcolor='b';
            end
            
            
            %% observed slope =
            scAx = 1:size(datanowreal,2);
            if iregres==1
                scAx=scAx+.15;
            end
            scAx = repmat(scAx, [size(datanowreal,1),1]);
            
            %reshape into vector for plot
            scAxm = reshape(scAx, 1, size(scAx,1)*size(scAx,2));
            datanowSc = reshape(datanowreal, 1, size(scAxm,1)*size(scAxm,2));
            
            
            
            %         take  regression and plot coeff. and slope
            %         p1 is the slope, p2 the intercept
            
            %         [p, S] = polyfit(scAxm, datanowSc, 1); %linear fit
            [p, S] = polyfit(scAxm, datanowSc, 2); %2nd order polinomial
            
            
            if iregres==2 % store the observed slopes for shuffled data.
                ObservedSlope=p(1);
            end
            
            f1=polyval(p,scAxm);
            hold on;
            
        end
        %%
        subplot(2,3,iDV+3)
        allslopes = zeros(1,size(shuffledtoplot,2));
        for ishuff=1:size(shuffledtoplot,2)
            
            datanow=squeeze(shuffledtoplot(:,ishuff,usedata));
            
            %% observed slope =
            scAx = 1:size(datanow,2);
            scAx = repmat(scAx, [size(datanow,1),1]);
            
            %reshape into vector for plot
            scAxm = reshape(scAx, 1, size(scAx,1)*size(scAx,2));
            datanowSc = reshape(datanow, 1, size(scAxm,1)*size(scAxm,2));
            
            %we do need to remove Nans to use polyfit however:
            nanid=find(isnan(datanowSc));
            scAxm(nanid)=[];
            datanowSc(nanid)=[];
            
            %take linear regression and plot coeff. and slope
            %         [p, S] = polyfit(scAxm, datanowSc, 2); %linear fit
            [p, S] = polyfit(scAxm, datanowSc, 2); %2nd order
            %p1 is the slope, p2 the intercept
            
            %         f1=polyval(p,scAxm);
            %         hold on; plot(scAxm, f1)
            allslopes(ishuff)=p(1);
        end
        %plot all.
        shg
        %%
        hb=histogram(sort(allslopes),25);
        hb.FaceColor=[.5 .5 .5];
        
        % fit CDF
        %note that if we have positive and negative numbers, this will
        %error.
        if any(hb.Data>0) && any(hb.Data<0)
            % shift all values so CDF is correct
            hb4cdf= hb.Data+1000;
            
        else
            hb4cdf=hb.Data;
        end
        cdf= cumsum(hb4cdf)/ sum(hb4cdf);
        %the X values (actual CV) corresponding to .05
        
        [~,cv05uncorr] = (min(abs(cdf-.95)));
        [~, c2] = min(abs(hb.Data-ObservedSlope)); %closest value.
        pvalis= 1-cdf(c2);
        
        % height of plot?
        yt=get(gca, 'ylim');
        %plot
        hold on;
        p05=plot([hb.Data(cv05uncorr) hb.Data(cv05uncorr)], [0, yt(2)*.7], ['k:'], 'linew', 3);
        
        text(hb.Data(cv05uncorr), yt(2)*.75, '95%', 'fontsize', 15, 'color', 'k');
        
        po=plot([ObservedSlope,ObservedSlope], [0, yt(2)*.7], ['-b'], 'linew', 3);
        
        if iDV==3
            legend([po], [{'Observed'}])
        end
        %adjust xlim.
        mX=min(hb.Data); maxX= ObservedSlope;
        xlim([ mX ObservedSlope + (ObservedSlope - mX)*.2])
        
        xlabel(['quadratic fit']);
        ylabel('count')
        
        set(gca, 'fontsize', 23)
        pvalswere(iDV) = pvalis;
        
    end
    %%
    set(gcf, 'color', 'w')
    
    disp(['p values are ' num2str(pvalswere )])
    %%
    cd([basedir])
    %%
    cd ../
    cd(['Figures' filesep 'Behavioural Bar Plots'])
    %%
    
    print('-dpng', 'PFI Fit and shuffle')
end

