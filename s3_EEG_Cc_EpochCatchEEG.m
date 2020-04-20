
addpath('/Users/MattDavidson/Desktop/SSVEP-Wills')
cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA')

%
basefol=pwd;


dbstop if error

%% job list

job.sortTrialindicesbyFreqandLocation=0; %store relevant catch/target loc data.

job.EpochperppantCATCH=0; %Epoch each participants Catch windows
%also stores the BP trace for next job.


%%




job.erpimageCATCHusingSNR=0; % saves at ppant level. % done. for 1f TG and BG
job.ppantCATCH_topotime=0;

job.concaterpdataacrossppants=1;
job.concater_invis_pdataacrossppants=0;

job.erpimageacrossppants=0;
job.BPandSSVEPtimecourseacrossppants=0;


% job.calcRESSforcatchtrials=1; % also performs a baseline removal.
% job.concatacrossppants=1;
% job.plotacrossppants=1;

ntrials=48;

%%
param_spcgrm.tapers = [1 1];
param_spcgrm.Fs= [250];
param_spcgrm.Fpass= [0 50];
param_spcgrm.trialave=0;
param_spcgrm.pad = 2;
% movingwin=[2,.15]; % increase to 2.5?
movingwin=[2.5,.25]; % 

%% based on/replaces Mkat_MT_taper from Irene.
cd(basefol)
cd('Behaviour')
%% load data to determine physical catch timing.
load('MD_AllBP_perppant')
load('Catchperformance') ; %to identify failed catches. (no BP).
%%
cd([basefol filesep 'EEG'])
pdirs = dir([pwd filesep '*_*' 'EEG']);
%%
%remaining participants after behavioral data analysis/exclusion
% allppants=[1,2,4,6,9:16,18]; %
 allppants=[1,2,4,6,7,9:19]; % NEW ppants.

% window=[-2 4];
window=[-3 3];

srate=250;

rmvbaseEPOCH=0; %0 for no removal,  1 for trial baseline, 2 for pooled freq baseline (compare across conds).
baserem1 = [-2.99 -2]; % baseline to remove in secs for ONSET of catch
baserem2 = baserem1;%[2 2.99]; % baseline to remove in secs for offset of catch
epochdur = sum(abs(window))*srate;
onsetc = ceil(epochdur)/2;
% peakfreqsare=[20,40]; %hz
%timing
tt = 0:1/srate:60;


if job.sortTrialindicesbyFreqandLocation==1
    
    if job.sortTrialindicesbyFreqandLocation==1
        %% This stores the 'like' target location x Hz combinations, for RESS analysis across like trials,
        % which needs spatial stability.
        %%
        cd([basefol filesep 'Behaviour'])
        
        load('Raw_Behavioral_AllPP.mat')
        
        %%
        %save in  EEG folder since that's where we will use it.
        cd([basefol filesep 'EEG'])
        pdirs = dir([pwd filesep '*_*' 'EEG']);
        for ifol=1%:19
            
            cd([basefol filesep 'EEG'])
            %redirect to ppant folder.
            cd(pdirs(ifol).name)
            
            
            %preallocate
            catch_location_x_index=zeros(4,30);  
            
            for Loc = 1:4 %for each location, (TL,TR,BL,BR) , store trial indices
                %by Hz, (8, 13,15,18)
                
                % 48 rows per participant.
                corr_rows = [1:48] + (48*(ifol-1));
                
                
                    %the columns for target location start at 6. Find
                    %wheneach location was missing.
                    
                    ftrials=find(AllTrialsAllPPBottomLeft(corr_rows,(Loc+5))== 1);
                    catch_location_x_index(Loc,:)= ftrials;
                    
                
            end
            catch_index_x_nremoved= sum(AllTrialsAllPPBottomLeft(corr_rows,6:9),2);
            save('TrialIndicesbyCatchLocationandNum','catch_location_x_index','catch_index_x_nremoved')
        end
    end
    
    
end
if job.EpochperppantCATCH==1
    %%
    cd([basefol filesep 'Behaviour'])
    
    load('Catchperformance.mat', 'catchStruct') % note that for some trials no catch occured (and screen froze!).
    %%
    
    for ifol = allppants
        
        cd([basefol filesep 'EEG'])        
        cd(pdirs(ifol).name)
        
        %load new rereferenced data.
        %         load('ppantRESSwholetrial')
        
        load(['P' num2str(ifol) '_autopreprocd.mat'])
        
        %there are four separate types of catch we can  store:
        %.1=all catches at on/offset (regardless of response),
        %.2=failed catches at on/offset (target absence no response),
        %.3= BP centred catches (changes timing).
        
        
        %.1
        ppant_SNREEG_catchramponset=nan(ntrials,64,epochdur+1);
        ppant_SNREEG_catchrampoffset = nan(ntrials,64,epochdur+1);
        
        %.2
        ppant_SNREEG_Disapfailedcatches=nan(ntrials,64,epochdur+1);
        ppant_SNREEG_Reapfailedcatches=nan(ntrials,64,epochdur+1);
        
        %.3
        ppant_SNREEG_disapBPwithincatch = nan(ntrials,64,epochdur+1);
        ppant_SNREEG_reapBPaftercatch = nan(ntrials,64,epochdur+1);
        
        %.4 
        ppant_SNREEG_invisiblecatchonset= nan(ntrials,64,epochdur+1);
        ppant_SNREEG_invisiblecatchoffset = nan(ntrials,64,epochdur+1);
        
        
        %how many did this ppant pass/fail
        ppantrows = [1:48] + (48*(ifol-1)) + 1;
        %         tmp=str2num(cell2mat(catchStruct(ppantrows,6)));
        %         goodtrials = ntrials- sum(tmp);
        %         badtrials = sum(tmp);
        
        
        goodcounter=1;
        badcounter=1;
        
        catchonsetBPs=nan(ntrials,361);
        catchoffsetBPs=nan(ntrials,361);
        withincatchonsetBPs=nan(ntrials,361);
        postcatchoffsetBPs=nan(ntrials,361);
        invisibleonsetBPs= nan(ntrials, 361);
        
        firstBPwithincatch_RTs=nan(ntrials,1);
        firstBPaftercatch_RTs=nan(ntrials,1);
        
        for itrial=1:ntrials
            
            
            ppant_SNREEG_catchramponset_tmp=[];
            ppant_SNREEG_catchrampoffset_tmp=[];
            ppant_SNREEG_disapBPwithincatch_tmp=[];
            ppant_SNREEG_reapBPaftercatch_tmp=[];
            
            ppant_SNREEG_invisibleonset_tmp=[];
            ppant_SNREEG_invisibleoffset_tmp=[];
            
            
            
            
            
            
            trialdata= ppantTrialDatawithDetails(ifol).TrialDetails(itrial);
            
            
            catchrampstart = trialdata.Catchstart_frames/60;
            catchrampend =  trialdata.Catchend_frames/60; %perhaps remove the ramp?
            
            
            %which target frequencies were the 'catch' in this trial?
            catchlocs = [trialdata.TL_Catch, trialdata.TR_Catch, trialdata.BL_Catch, trialdata.BR_Catch];
            
            %any invisible catches this trial?
            invislocs= trialdata.invisiblecatchlocations;
            
            
            
            % take 64 channel data for this EPOCH.
            tmpEpoch = squeeze(EEG_preprocd(:,:,itrial)); %both 20 and 40 hz.
            %%
            
            %initialize
            tmpCatchONsetOFFset=[];
            if sum(catchlocs)~=0 %skip trial.
            %%%%%%
            %%%%%%
            % First collect onset and offset EEG epochs. (not accounting for BP)
            %%%%%%%
            %%%%%%%
              catchlocs= find(catchlocs);
            for iwin=1:2
                switch iwin
                    case 1
                        centrep = catchrampstart;
                    case 2
                        centrep = catchrampend;
                end
                
                %prep data length
                windowstart = centrep-abs(window(1));                                
                [~, rampst] = min(abs(tt-windowstart));                
                rampend=rampst+epochdur; %add epoch length.
                
                %collect EEG
                switch iwin
                    case 1 %onsets
                        
                        
                        % was the catch invisible?
                        if ~isempty(invislocs) % if it was invisible:
                            ppant_SNREEG_invisibleonset_tmp=tmpEpoch(:,rampst:rampend);
                            ppant_SNREEG_catchramponset_tmp =  nan(64,epochdur+1);
                            
                            
                            % store the BP also
                            try invisibleonsetBPs(itrial,:)= squeeze(nansum(trialdata.BPData(catchlocs,ceil(catchrampstart*60-3*60):ceil(catchrampstart*60+3*60)),1));
                            catch
                                invisibleonsetBPs(itrial,:)= squeeze(nansum(trialdata.BPData(catchlocs,ceil(catchrampstart*60-3*60):ceil(catchrampstart*60+3*60+1)),1));
                            end
                            
                            catchonsetBPs(itrial,:) = nan(1,361);
                            
                            
                            
                        else % if it was not.
                            ppant_SNREEG_catchramponset_tmp = tmpEpoch(:,rampst:rampend);
                            ppant_SNREEG_invisibleonset_tmp =  nan(64,epochdur+1);
                            
                            % store the BP also
                            try catchonsetBPs(itrial,:)= squeeze(nansum(trialdata.BPData(catchlocs,ceil(catchrampstart*60-3*60):ceil(catchrampstart*60+3*60)),1));
                            catch
                                catchonsetBPs(itrial,:)= squeeze(nansum(trialdata.BPData(catchlocs,ceil(catchrampstart*60-3*60):ceil(catchrampstart*60+3*60+1)),1));
                            end
                            
                            invisibleonsetBPs(itrial,:) = nan(1,361);
                            
                            
                        end
                        
                        
                        
                    case 2
                        %adjust for late catch events
                        if rampend>length(tmpEpoch)
                            nantrain = nan(64,rampend-length(tmpEpoch));
                            tmpcat=horzcat(tmpEpoch, nantrain);
                        
                        ppant_SNREEG_catchrampoffset_tmp = tmpcat(:,rampst:rampend);
                        else
                            ppant_SNREEG_catchrampoffset_tmp = tmpEpoch(:,rampst:rampend);
                        end
                        
                end
            end
            
            % store the BPs associated
                %catch aligned
              
                try catchoffsetBPs(itrial,:)= squeeze(nansum(trialdata.BPData(catchlocs,ceil(catchrampend*60-3*60):ceil(catchrampend*60+3*60)),1));
                catch
                    %extend in case a late catch occurred.
                    if (catchrampend*60 + 3*60 +1) > length(trialdata.BPData)
                        nantrain = nan(4,length(trialdata.BPData)-catchrampend*60);
                        tmpcat=horzcat(trialdata.BPData, nantrain);                        
                    catchoffsetBPs(itrial,:)= squeeze(nansum(tmpcat(catchlocs,ceil(catchrampend*60-3*60):ceil(catchrampend*60+3*60)),1));
                    else
                        catchoffsetBPs(itrial,:)= squeeze(nansum(trialdata.BPData(catchlocs,ceil(catchrampend*60-3*60):ceil(catchrampend*60+3*60+1)),1));
                    end
                end
            
            
            
           
              
                %%%%%%%
                % now epoch around BP , first within catch, then after.
                % Was there a BP for catch? ie catch (pass/fail)?
                %%%%%%
                
                rowis=itrial + (ntrials*(ifol-1)) + 1;
                badcatch = str2num(catchStruct{rowis, 6}); %6th column for reject or no.
                
                if badcatch~=1 %then collect BP centred data also
                    
                    CatchBPs = nansum(trialdata.CatchBPs_exact,1);
                    
                    %if BPdata during catch removal
                    firstBPwithincatch= min(find(diff(CatchBPs)>0));
                    
                    if firstBPwithincatch>0 %change in BP state (not already pressed)
                        
                        
                        
                        
                        %Epoch around this point (notice catch gone)
                        
                        %                     correct for difference in timing units.
                        BP_disapwithincatch_secs = firstBPwithincatch/60 + catchrampstart;
                        
                        [~, BPdisap1] = min(abs(tt-(BP_disapwithincatch_secs-abs(window(1)))));
                        BPdisap2 = BPdisap1+epochdur;
                        
                        %adjust for late events
                        if BPdisap2 > length(tmpEpoch)
                            nantrain=nan(64, BPdisap2-length(tmpEpoch));
                            tmpcat = horzcat(tmpEpoch, nantrain);
                        ppant_SNREEG_disapBPwithincatch_tmp= tmpcat(:,BPdisap1:BPdisap2);
                        else
                            ppant_SNREEG_disapBPwithincatch_tmp= tmpEpoch(:,BPdisap1:BPdisap2);
                        end
                        
                        firstBPwithincatch_RTs(itrial)=firstBPwithincatch;
                        
                        
                        %store associated button press Epoch.
                        try withincatchonsetBPs(itrial,:)= squeeze(nansum(trialdata.BPData(catchlocs,ceil(BP_disapwithincatch_secs*60-3*60):ceil(BP_disapwithincatch_secs*60+3*60)),1));
                        catch
                            if BP_disapwithincatch_secs*60+3*60 > length(trialdata.BPData)
                                     nantrain = nan(4,length(trialdata.BPData)-BP_disapwithincatch_secs*60+3*60);
                                     tmpcat=horzcat(trialdata.BPData, nantrain);
                                     withincatchonsetBPs(itrial,:)= squeeze(nansum(tmpcat(catchlocs,ceil(BP_disapwithincatch_secs*60-3*60):ceil(BP_disapwithincatch_secs*60+3*60)),1));
                            else
                                withincatchonsetBPs(itrial,:)= squeeze(nansum(trialdata.BPData(catchlocs,ceil(BP_disapwithincatch_secs*60-3*60):ceil(BP_disapwithincatch_secs*60+3*60+1)),1));
                            end
                        end
                        
                        
                        
                        
                    else %no response to catch?
                        ppant_SNREEG_disapBPwithincatch_tmp= nan(64,epochdur+1);
                        firstBPwithincatch_RTs(itrial)=nan;
                    end
                    
                    %%%%% also collect and epoch around first BP after catch
                    %%%%% return.
                    
                    %note that the 'long epoch' is 2s before catch start and
                    %8seconds after catch START
                    
                    catchdur=trialdata.totalCatchdur;
                    
                    longcatchEpoch=nansum(trialdata.CatchBPs_longepoch,1);
                    
                    %catch end frames, within this epoch = diff +2s for the
                    %onset
                    postOffset=longcatchEpoch(1,(catchdur+120):end);
                    
                    %first change in BP, (could also use return to zero?)
                    
                    firstBPaftercatch= min(find(diff(postOffset)<0));
                    
                    %                     firstBPaftercatch= min(find(postOffset~=0));
                    
                    %ajust this time, to account for seconds in whole trial EEG
                    
                    catchreturn_secs = trialdata.Catchend_frames/60 + firstBPaftercatch/60;
                    
                    try %will fail if BPdisap2 is beyond the end of the trial.
                        [~, BPdisap1] = min(abs(tt-(catchreturn_secs-abs(window(1)))));
                        BPdisap2 = BPdisap1+epochdur;
                        
                        ppant_SNREEG_reapBPaftercatch_tmp= tmpEpoch(:,BPdisap1:BPdisap2);
                        
                        firstBPaftercatch_RTs(itrial)=firstBPaftercatch;
                        
                        
                        
                        %store BP aligned EPOCH
                        
                        
                        try postcatchoffsetBPs(itrial,:)= squeeze(nansum(trialdata.BPData(catchlocs,ceil(catchreturn_secs*60-3*60):ceil(catchreturn_secs*60+3*60)),1));
                        catch
                            postcatchoffsetBPs(itrial,:)= squeeze(nansum(trialdata.BPData(catchlocs,ceil(catchreturn_secs*60-3*60):ceil(catchreturn_secs*60+3*60+1)),1));
                        end
                        
                        
                    catch
                        ppant_SNREEG_reapBPaftercatch_tmp= nan(64,epochdur+1);
                        firstBPaftercatch_RTs(itrial)=nan;
                    end
                    
                    
                    %%%%%%%%
                    
                    
                    goodcounter=goodcounter+1;
                else %just record the failed catch EEG, no BP window
                    
                    ppant_SNREEG_reapBPaftercatch_tmp= nan(64,epochdur+1);
                    ppant_SNREEG_disapBPwithincatch_tmp= nan(64,epochdur+1);
                    %                 badcounter=badcounter+1;
                end
                
                
                % STORE and save EEG per ppant.
                
                
                ppant_SNREEG_catchramponset(itrial,:,:)=ppant_SNREEG_catchramponset_tmp;
                
                
                ppant_SNREEG_catchrampoffset(itrial,:,:)=ppant_SNREEG_catchrampoffset_tmp;
                
                ppant_SNREEG_disapBPwithincatch(itrial,:,:)= ppant_SNREEG_disapBPwithincatch_tmp;
                
                ppant_SNREEG_reapBPaftercatch(itrial,:,:)= ppant_SNREEG_reapBPaftercatch_tmp;
                
                ppant_SNREEG_invisiblecatchonset(itrial,:,:)=ppant_SNREEG_invisibleonset_tmp;
                
                
                
            end
            
        end
        
        savename = 'ppant_Catch_Epoched';
        
        save(savename, ...
            'ppant_SNREEG_catchrampoffset',...
            'ppant_SNREEG_catchramponset',...
            'ppant_SNREEG_disapBPwithincatch',...
            'ppant_SNREEG_reapBPaftercatch',...
            'ppant_SNREEG_invisiblecatchonset',...            
            'window','tt', 'firstBPaftercatch_RTs', 'firstBPwithincatch_RTs', ...
            'catchonsetBPs', 'catchoffsetBPs',...
            'invisibleonsetBPs',...
            'withincatchonsetBPs', 'postcatchoffsetBPs')
        
                disp(['Finished ppant ' num2str(ifol)])
    end
end



%%





if job.erpimageCATCHusingSNR==1
    getelocs
    window=[-3 3];
    
    srate=250;
    
    epochdur = sum(abs(window))*srate;
    
    timeid = [0:1/srate:epochdur];
    timeid= timeid-3;
    
    onsetc = ceil(epochdur)/2;
    
    %timing
    tt = 0:1/srate:60;
    
    %%
    param_spctrm.tapers = [1 1];
    param_spctrm.Fs= [250];
    param_spctrm.Fpass= [0 50];
    param_spctrm.trialave=0;
%     
%     param_spcgrm.tapers = [1 1];
%     param_spcgrm.Fs= [250];
%     param_spcgrm.Fpass= [0 50];
%     param_spcgrm.trialave=0;
%     movingwin=[1,.15];
    
    rmvbase=0;
    snrmethod=2; % 1 = division (slow), 2=dot point kernel.
    %
    for ifol = allppants
        
        
        for hzis=5:6
            usehz=peakfreqsare(hzis);
        
            for alignment=1:2
                
                cd(basefol)
                cd('EEG')
                cd(pdirs(ifol).name)
                
                load(['ppant_Catch_Epoched'])
%                 load('ProposedCATCHonsettiralindextorej')
                for itimezero = 1:3
                    switch itimezero
                        case 1
                        
                        
                        if alignment==1
                            datatouse = ppant_SNREEG_catchramponset;
                            RTstouse = firstBPwithincatch_RTs;
                            BPstouse = catchonsetBPs;
                            ctype = 'onset visible';
                        else
                            
                            datatouse = ppant_SNREEG_disapBPwithincatch;
                            RTstouse = firstBPwithincatch_RTs;
                            BPstouse = withincatchonsetBPs;
                            ctype = 'onset, BP aligned';
                        end
                        bsrem = [-3 -1]; %seconds
                        
                        case 2
                        
                        if alignment==1
                            datatouse = ppant_SNREEG_catchrampoffset;
                            RTstouse = firstBPaftercatch_RTs;
                            BPstouse = catchoffsetBPs;
                            ctype = 'offset';
                        else
                            datatouse = ppant_SNREEG_reapBPaftercatch;
                            RTstouse = firstBPaftercatch_RTs;
                            BPstouse = postcatchoffsetBPs;
                            ctype = 'offset, BP aligned';
                            
                            
                        end
                        
                        
                        bsrem = [1 3]; %seconds
                        
                        case 3 
                            % only perform once.
                           if alignment==2
                               break
                           else
                               
                                datatouse = ppant_SNREEG_invisiblecatchonset;
                                BPstouse = invisibleonsetBPs;
                                ctype = 'onset invisible';
                           end
                            
                    end
                    
                    
                    %plot the topo for sanity check: (pre disap)
                    windcheck= [-3 -0.1];
                    tidx=dsearchn(timeid', [windcheck]');
                    
                    %reduce size.
                    topod= squeeze(datatouse(:,:,tidx(1):tidx(2)));
                    %calc statSNR per channel, over trials
                    
                    
                    %%
                    snr20=zeros(64,1);
                    
                    for ichan= 1:64
                        
                        [s,f]=mtspectrumc(squeeze(topod(:,ichan,:))', param_spctrm)    ;
                        
                        %comp SNR;
                        kernelw = [-.25 -.25 0 0 1 0  0 -.25 -.25];
                        
                        sNOW=zeros(size(s));
                        for tmptr=1:size(s,2)
                            
                            sNOW(:,tmptr)=conv(log(s(:,tmptr)), kernelw, 'same');
                        end
                        
                        %store just stim freq.
                        [~, hzid]= min(abs(f-usehz));
                        
                        snr20(ichan)=squeeze(nanmean(sNOW(hzid,:),2));
                    end
                    %%
                    %             [~,usechan]=max(snr20);
                    usechan=62; %POz for all
                    
                    
                    %
                    figure(1)
                    clf
                    subplot(4,2,1:2)
                    topoplot(snr20, elocs(1:64), 'emarker2', {usechan, 'o', 'w'});
                    c=colorbar;
                    caxis([0 2])
                    ylabel(c, 'SNR')
                    title({[num2str(usehz) 'Hz SNR '];['-3000:-100ms']})
                    
                    
                    
                    %plot BP for comparison
                    subplot(4,2, [3 5])
                    %sort trial index by longest RT first.
                    %sort by duration?
                    
                    %                 [sortedRTs, cid] = sort(RTstouse, 'descend');
                    %
                    %sort by accum Targ gone.
                    %sort by accum PFI in window after onset.
                    if itimezero==1 || itimezero==3
                        accumPeriod = sum(BPstouse(:,180:end),2);
                    else
                        accumPeriod = sum(BPstouse(:,1:180),2);
                    end
                    
                    
                    
                    %%
                    [checkBP, cid] = sort(accumPeriod, 'descend');
                    
                    % MOVE nan to bottom of order (not top)
                    nanind= find(isnan(checkBP));
                    %new start point
                    cid=[cid(length(nanind)+1:end) ; cid(nanind)];
                    %                     checkBP=[checkBP(length(nanind)+1:end); checkBP(nanind)]
                    
                    %rearrange.
                    BPstousesorted=BPstouse(cid,:);
                    dataistmp=squeeze(datatouse(cid,usechan,:));
                    
                    %make sure no Nans in SNR either.
                    nanind= find(isnan(mean(dataistmp,2)));
                    %move to bottom.
                    tmpD= dataistmp;
                   
                    
                   allt=1:size(dataistmp,1);
                   keept= setdiff(allt,nanind);
                   rmv = tmpD([nanind],:);
                   kpD=tmpD(keept,:);
                   newD = [kpD;rmv];     
                   datais=newD;
                    %%
                    figure(1)
                    
                    imagesc(-3:1/60:3,1:ntrials,BPstousesorted)
                    c=colorbar;
                    ylabel(c, 'buttons pressed')
                    set(gca, 'ytick', 1:ntrials, 'yticklabel', cid)
                    ylabel('Trialind')
                    xlabel('time [sec]')
                    hold on
                    plot([0 0] , ylim, ['w:'], 'linewidth', 3)
                    title({['sorted BP data for catch ' ctype ','];['ppant' num2str(ifol)]})
                    set(gca, 'fontsize', 10)
                    %% show mean over time.
                    subplot(4,2,7)
                    plot([-3:1/60:3], nanmean(BPstousesorted,1),['k-'], 'linewidth', 3)
                    %         shadedErrorBar([-3:1/60:3], nanmean(acrBP,1),nanstd(acrBP,1))
                    set(gca, 'fontsize', 15)
                    hold on;
                    ylabel('nanmean BP ')
                    %
                    xlabel('Time secs')
                    axis tight
                    plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                    
                    
                    set(gca, 'fontsize', 10)
                    %%
                    % now plot spectrogram SNR.
                    
                    
                    %
                    %             %rmvbaseline from EEG.
                    bsdata = zeros(size(datais));
                    
                    for itrial = 1:size(datais,1)
                        td = detrend(datais(itrial,:), 'linear');
                        tdrm= mean(td(1,1:250));
                        rmb= repmat(tdrm, [1 length(td)]);
                        bsdata(itrial,:) = td-rmb;
                    end
                    datais=bsdata;
                    %%
                    [sgrm ,tgrm, fgrm] = mtspecgramc(datais', movingwin, param_spcgrm);
                    %conv SNR
                snr_sgrm =zeros(size(sgrm));
                
                 tmps=sgrm;
                            %adjust for HBW 
                                k = ((param_spcgrm.tapers(1,1)));%
                                hbw = (k+1)./ (2.*[movingwin(1,1)]);
                                neighb = 2; %hz
                                
                                %space out kernel appropriately
                                distz = dsearchn(fgrm', [hbw, neighb, neighb*2+hbw]');
                                
                                % we don't want odd numbers, messes with
                                % calculations, so round allup to even.
                                ods = find(mod(distz,2)~=0);
                                %adjust.
                                distz(ods)= distz(ods)+1;
                                
                                %make kernel
                                kernelneighb = -1*(repmat(1/(distz(2)*2), [1, distz(2)]));
                                kernelskip = zeros(1,distz(1));
                                
                                kernelw= [kernelneighb, kernelskip, 1, kernelskip, kernelneighb];
                                if sum(kernelw)~=0
                                    error('')
                                end
                            if snrmethod==1%method 1, using division.
                                
                               
                                
                                
                                
                                % loop over frequencies and compute SNR
                                numbins = distz(2);
                                skipbins = distz(1);
                                snr_sgrm= zeros(size(sgrm));
                                for itrial = 1:size(sgrm,3)
                                for itime=1:size(tmps,1)
                                    for hzi=numbins+1:length(fgrm)-numbins-1
                                        numer = tmps(itime,hzi, itrial);
                                        denom = nanmean( tmps(itime,[hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
                                        snr_sgrm(itime,hzi,itrial) = numer./denom;
                                        
                                    end
                                end
                                end
                                
                            else
                                %                             %conv SNR
                                %                             snr_sgrm =zeros(size(sgrm));
                                %
                                
                                %%
%                                 kernelw= [-1/8 -1/8 -1/8 -1/8 0 0 1 0 0 -1/8 -1/8 -1/8 -1/8];
                                
                                %snr on trials or
                                for itrial=1:size(sgrm,3)
                                    %compute SNR
                                    tmps= squeeze(sgrm(:,:,itrial));
                                    
                                    for itime= 1:size(tmps,1)
                                        checkput = conv(log(tmps(itime,:)), kernelw,'same');
                                        if ~isreal(checkput)
                                            snr_sgrm(itime,:,itrial)= nan(1, size(tmps,2));
                                        else
                                            snr_sgrm(itime,:,itrial)= conv(log(tmps(itime,:)), kernelw,'same');
                                        end
                                    end
                                end
                            
                             %% sanity check:
%                               figure(1);
%                               %plot specgrm and spctrm.
%                               subplot(221)
%                             imagesc(tgrm, fgrm, squeeze(mean(log(sgrm),3))'); colorbar; caxis([0 2])
%                             subplot(223)
%                             mtrials=squeeze(mean(log(sgrm),3));
%                             plot(fgrm, squeeze(mean(mtrials,1)));
%                             %plot specgrm and spctrm.
%                               subplot(222)
%                             imagesc(tgrm, fgrm, squeeze(mean((snr_sgrm),3))'); colorbar; caxis([0 2])
%                             subplot(224)
%                             mtrials=squeeze(mean((snr_sgrm),3));
%                             plot(fgrm, squeeze(mean(mtrials,1)));
%                             hold on; 
%                                % center kernel at 20 Hz
%     hlfkern= (length(kernelw)-1)/2;
%     [~,hzid]= min(abs(fgrm-20)); 
%     fvector = (hzid- hlfkern):(hzid+hlfkern);
%     plot(fgrm(fvector), kernelw*2, 'r')
                            end
                    %%
                    %store just stim freq.
                        [~, hzid]= min(abs(fgrm-usehz));
                    %reduce size.
                    snrgrm20=squeeze(snr_sgrm(:,hzid,:))';
                    %
                    
                    tbase = tgrm-3;
                    tidFREEZE = dsearchn(tbase', [-.15 .35]');
                    
                    if rmvbase==1
                        %rmv baseline ?
                        acrSNRsort_rmv=zeros(size(snrgrm20));
                        for itrial=1:size(snrgrm20,1);
                            
                            tmp=snrgrm20(itrial,:);
                            %which baseline?
                            tidx= dsearchn(tbase', [bsrem]');
                            
                            %rmvbase
                            bs= mean(tmp(1,tidx(1):tidx(2)));
                            bs=repmat(bs, 1, length(tmp));
                            acrSNRsort_rmv(itrial,:)= tmp-bs;
                            %                         acrSNRsort_rmv(itrial,:)= tmp./bs;
                            
                        end
                        snrgrm20=acrSNRsort_rmv;
                        
                        
                        
                        
                        
                        
                    end
                    %
                    
                    % smooth across trials
                    sm_snrgrm20=zeros(size(snrgrm20));
                    for itime=1:size(snrgrm20,2)
                        sm_snrgrm20(:,itime)= smooth(snrgrm20(:,itime),3, 'moving');
                    end
                    
                    subplot(4,2,[4 6]);
                    imagesc(tgrm-3,  1:size(snr_sgrm,3), sm_snrgrm20);
%                     imagesc(tgrm-3,  1:size(snr_sgrm,3), snrgrm20);
                    c=colorbar;
                    ylabel(c, 'log(SNR)')
                    xlabel('Time secs')
                    caxis([-1*max(max(sm_snrgrm20)) max(max(sm_snrgrm20))])
                    %                     caxis([-10 10])
                    set(gca, 'ytick', 1:ntrials, 'yticklabel', cid)
                    hold on
                    %%
                    plot([0 0] , ylim, ['w:'], 'linewidth', 3)
                    title({['sorted ' num2str(usehz) 'Hz SNR data (smoothed) for catch ' ctype ','];['ppant' num2str(ifol)]})
                    set(gca, 'fontsize', 10)
                    
                    %% show mean over time.
                    subplot(4,2,8)
                    plot(tbase, nanmean(snrgrm20,1),['k-'], 'linewidth', 3)
                    %         shadedErrorBar([-3:1/60:3], nanmean(acrBP,1),nanstd(acrBP,1))
                    set(gca, 'fontsize', 15)
                    hold on;
                    ylabel('nanmean SNR (basesub) ')
                    
                    axis tight
                    plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                    %
                    xlabel('Time secs')
                    set(gcf,'color', 'w');
                    
                    %         ylim([0 20])
                    
                    %%
                    cd(basefol)
                    cd('Figures')
                    cd('Ppant rawSNR summary (Catch)')
                    
                    %%
                    print('-dpng', ['figure_catch_' ctype '_summary_BPandSSVEP_' num2str(usehz) ' ppant ' num2str(ifol) '.png']);
                    
                    %sanitychecks
%                     figure(2)
%                     clf
%                         for ipl=1:ntrials
%                             subplot(6,4,ipl), 
%                             plot(tgrm-3, snrgrm20(ipl,:));
%                             title([num2str(cid(ipl))]);
%                             ylim([-20 20]); hold on;
%                             plot([0 0], ylim, ['k:']);
%                         end
                    %%
                    switch itimezero
                        case 1 %store for across ppant plots:
                            % store based on alignment
                            
                            Ppant_onsetBP=BPstousesorted;
                            Ppant_onsetSNR=snrgrm20; %sorted.
                            Ppant_onsetRTs=RTstouse;
                            Ppant_onsetTOPO=snr20;
                        case 2
                            
                            Ppant_offsetBP=BPstousesorted;
                            Ppant_offsetSNR=snrgrm20; %always sorted in descending order of PFI.
                            Ppant_offsetRTs=RTstouse;
                            Ppant_offsetTOPO=snr20;
                            
                        case 3 % invisible case
                            Ppant_invisibleonsetBP=BPstousesorted;
                            Ppant_invisibleonsetSNR=snrgrm20; %sorted.
                    end
                end
                
                
                        savename=['Catchperformance_withSNR_' num2str(usehz)];
                
                if alignment==2
                    savename = [savename ',BPaligned'];
                end
                cd(basefol)
                cd('EEG')
                cd(pdirs(ifol).name)
                
                save(savename,...
                    'Ppant_onsetBP','Ppant_offsetBP',...
                    'Ppant_onsetSNR', 'Ppant_offsetSNR', ...
                    'Ppant_onsetRTs','Ppant_offsetRTs',...
                    'Ppant_invisibleonsetBP', 'Ppant_invisibleonsetSNR',...
                    'Ppant_onsetTOPO', 'Ppant_offsetTOPO', 'tgrm')
                
                
                
                
            end
            
        end
        
        
        
    end
end


if job.concaterpdataacrossppants==1
    
    ppantsmoothing=1; % average across participants, after smoothing, or no.
    
    for ihz=1:7
        loadname =['Catchperformance_withSNR_' num2str(peakfreqsare(ihz))];
       
        for ialignment=1:2
            
            storeacrossPpant_onsetBP=[];
            storeacrossPpant_onsetSNR=[];
            storeacrossPpant_onsetRTs=[];
            storeacrossPpant_onsetTOPO=[];
            storeacrossPpant_offsetBP=[];
            storeacrossPpant_offsetSNR=[];
            storeacrossPpant_offsetRTs=[];
            storeacrossPpant_offsetTOPO=[];
            
            icounter=1;
            
            
            
            if ialignment==2 
                loadname= [loadname ',BPaligned'];
            end
            
            for ippant = allppants;
                cd(basefol)
                cd('EEG')
                cd(pdirs(ippant).name);
                %onset types
                
                 
                load(loadname)
                
                
                
                if any(isnan(Ppant_onsetSNR(:,1)))
                    lc= min(find(isnan(Ppant_onsetSNR(:,1))));
                    lc=lc-1;
                else
                    lc=size(Ppant_onsetSNR,1);
                end
                
                % we need to resample the BP and SNR data, to equate across
                % trial types    (ignoring nan)
                % currently at 100 fs.
                %%
                Ppant_onsetBP= resample(Ppant_onsetBP(1:lc,:),100,lc);
                Ppant_onsetSNR= resample(Ppant_onsetSNR(1:lc,:),100,lc);
                
                 
                %offsets too
                if any(isnan(Ppant_offsetSNR(:,1)));
                    lc= min(find(isnan(Ppant_offsetSNR(:,1))));
                    lc=lc-1;
                else
                    lc=size(Ppant_offsetSNR,1);
                end
                
                % we need to resample the BP and SNR data, to equate across
                % trial types    (ignoring nan)
                % currently at 100 fs.
                %%
                Ppant_offsetBP= resample(Ppant_offsetBP(1:lc,:),100,lc);
                Ppant_offsetSNR= resample(Ppant_offsetSNR(1:lc,:),100,lc);
                
               % also apply participant level smoothing.
                if ppantsmoothing==1
                    %smooth across trials, per ppant for SNR
                    
                    
                    % note, need to adjust for edge artefacts in smoothig
                    % process, repmat the end trials so smoothing doesnt
                    % contaminate.
                    
                    
                    for ionoff=1:4
                        
                        switch ionoff
                            case 1
                                pData = Ppant_onsetSNR;
                            case 2
                                pData = Ppant_offsetSNR;
                            case 3
                                pData = Ppant_onsetBP;
                            case 4
                                pData = Ppant_offsetBP;
                        end
                        
                        
                           %new data with expanded width.
                            pdataN = zeros(size(pData,1)+32, size(pData,2));
                            %top of image
%                             pdataN(1:16,:) = repmat(pData(1,:), [16,1]);
                            pdataN(1:16,:) = flipud(pData(1:16,:));
                            %and bottom;
%                             pdataN(end-15:end,:) = repmat(pData(end,:), [16,1]);
                            pdataN(end-15:end,:) = flipud(pData(end-15:end,:));
                            %fill centre
                            pdataN(17:size(pData,1)+16,:) = pData;
                            
                            
                            %now safe to smooth
                            
                            pdataOUT= zeros(size(pdataN));
                            for itime=1:size(pdataN,2)
                                pdataOUT(:,itime)= smooth(pdataN(:,itime),15);
                            end
                            
                    
                    %now collapse to original size.
                    pdataOUT = pdataOUT(17:end-16,:);
                    
                    switch ionoff
                        case 1
                            Ppant_onsetSNR = pdataOUT;
                        case 2
                            Ppant_offsetSNR = pdataOUT;
                            case 3
                            Ppant_onsetBP = pdataOUT;
                        case 4
                            Ppant_offsetBP = pdataOUT;
                            
                    end
                    
                    end
                
                end
                
                
                %%
                figure(1); clf;
                
                subplot(221); 
                imagesc(Ppant_onsetBP); 
                
                subplot(223)
                imagesc(Ppant_onsetSNR);
                subplot(222)
                imagesc(Ppant_offsetBP); subplot(224)
                imagesc(Ppant_offsetSNR);
                %%
                
                
                storeacrossPpant_onsetBP(icounter,:,:)=Ppant_onsetBP;
                storeacrossPpant_onsetSNR(icounter,:,:)=Ppant_onsetSNR;
                storeacrossPpant_onsetRTs(icounter,:)=Ppant_onsetRTs;
                storeacrossPpant_onsetTOPO(icounter,:)=Ppant_onsetTOPO;
               
                
                 storeacrossPpant_offsetBP(icounter,:,:)=Ppant_offsetBP;
                storeacrossPpant_offsetSNR(icounter,:,:)=Ppant_offsetSNR;
                storeacrossPpant_offsetRTs(icounter,:)=Ppant_offsetRTs;
                storeacrossPpant_offsetTOPO(icounter,:)=Ppant_offsetTOPO;
                
                
                
                
                
                icounter=icounter+1;
            end
            
            %save appropriately
            cd(basefol)
            cd('EEG')
            cd('GFX_Pre-RESS')
            
            savename=['GFX_Catchperformance_withSNR_' num2str(peakfreqsare(ihz)) '_min0_RESS'];
            
            
            if ialignment==2
                savename = [savename ',BPaligned'];
            end
            save(savename,...
                'storeacrossPpant_onsetBP','storeacrossPpant_offsetBP',...
                'storeacrossPpant_onsetSNR', 'storeacrossPpant_offsetSNR', ...
                'storeacrossPpant_onsetRTs','storeacrossPpant_offsetRTs',...
                'storeacrossPpant_onsetTOPO', 'storeacrossPpant_offsetTOPO', 'tgrm')
            
            
            
        end
    end
    
    
end

if job.concater_invis_pdataacrossppants==1
    
    ppantsmoothing=1; % average across participants, after smoothing, or no.
    
    for ihz=1:2
        
        switch ihz
            case 1
            loadname='Catchperformance_withSNR_15';
            case 2
            loadname='Catchperformance_withSNR_20';
            case 3
                loadname='Catchperformance_withSNR_30';
            case 4
                loadname='Catchperformance_withSNR_40';                
        end
        for ialignment=1% just onset type for now.
            
            storeacrossPpant_onsetBP_invisible=[];
            storeacrossPpant_onsetSNR_invisible=[];
%             storeacrossPpant_onsetRTs=[];
%             storeacrossPpant_onsetTOPO=[];
%             storeacrossPpant_offsetBP=[];
%             storeacrossPpant_offsetSNR=[];
%             storeacrossPpant_offsetRTs=[];
%             storeacrossPpant_offsetTOPO=[];
%             
            icounter=1;
            for ippant = allppants;
                cd(basefol)
                cd('EEG')
                cd(pdirs(ippant).name);
                %onset types
                
                load(loadname)
                
                %%%%%%% %%%%% 
                Ppant_onsetSNR=Ppant_invisibleonsetSNR;
                Ppant_onsetBP=Ppant_invisibleonsetBP;
                
                if any(isnan(Ppant_onsetSNR(:,1)))
                    lc= min(find(isnan(Ppant_onsetSNR(:,1))));
                    lc=lc-1;
                else
                    lc=size(Ppant_onsetSNR,1);
                end
                
                % we need to resample the BP and SNR data, to equate across
                % trial types    (ignoring nan)
                % currently at 100 fs.
                %%
                Ppant_onsetBP= resample(Ppant_onsetBP(1:lc,:),100,lc);
                Ppant_onsetSNR= resample(Ppant_onsetSNR(1:lc,:),100,lc);
                
%skipping offsets.                 
%                 %offsets too
%                 if any(isnan(Ppant_offsetBP(:,1)));
%                     lc= min(find(isnan(Ppant_offsetBP(:,1))));
%                     lc=lc-1;
%                 else
%                     lc=size(Ppant_offsetBP,1);
%                 end
                
                % we need to resample the BP and SNR data, to equate across
                % trial types    (ignoring nan)
                % currently at 100 fs.
%                 %%
%                 Ppant_offsetBP= resample(Ppant_offsetBP(1:lc,:),100,lc);
%                 Ppant_offsetSNR= resample(Ppant_offsetSNR(1:lc,:),100,lc);
                
               % also apply participant level smoothing.
                if ppantsmoothing==1
                    %smooth across trials, per ppant for SNR
                    
                    
                    % note, need to adjust for edge artefacts in smoothig
                    % process, repmat the end trials so smoothing doesnt
                    % contaminate.
                    
                    
                    for ionoff=1:2
                        
                        switch ionoff
                            case 1
                                pData = Ppant_onsetSNR;
%                             case 2
%                                 pData = Ppant_offsetSNR;
                            case 2
                                pData = Ppant_onsetBP;
%                             case 4
%                                 pData = Ppant_offsetBP;
                        end
                        
                        
                           %new data with expanded width.
                            pdataN = zeros(size(pData,1)+32, size(pData,2));
                            %top of image
%                             pdataN(1:16,:) = repmat(pData(1,:), [16,1]);
                            pdataN(1:16,:) = flipud(pData(1:16,:));
                            %and bottom;
%                             pdataN(end-15:end,:) = repmat(pData(end,:), [16,1]);
                            pdataN(end-15:end,:) = flipud(pData(end-15:end,:));
                            %fill centre
                            pdataN(17:size(pData,1)+16,:) = pData;
                            
                            
                            %now safe to smooth
                            
                            pdataOUT= zeros(size(pdataN));
                            for itime=1:size(pdataN,2)
                                pdataOUT(:,itime)= smooth(pdataN(:,itime),15);
                            end
                            
                    
                    %now collapse to original size.
                    pdataOUT = pdataOUT(17:end-16,:);
                    
                    switch ionoff
                        case 1
                            Ppant_onsetSNR = pdataOUT;
%                         case 2
%                             Ppant_offsetSNR = pdataOUT;
                            case 2
                            Ppant_onsetBP = pdataOUT;
%                         case 4
%                             Ppant_offsetBP = pdataOUT;
                            
                    end
                    
                    end
                
                end
                
                
                
                
                storeacrossPpant_onsetBP_invisible(icounter,:,:)=Ppant_onsetBP;
                storeacrossPpant_onsetSNR_invisible(icounter,:,:)=Ppant_onsetSNR;
%                 storeacrossPpant_onsetRTs(icounter,:)=Ppant_onsetRTs;
%                 storeacrossPpant_onsetTOPO(icounter,:)=Ppant_onsetTOPO;
               
                
%                  storeacrossPpant_offsetBP(icounter,:,:)=Ppant_offsetBP;
%                 storeacrossPpant_offsetSNR(icounter,:,:)=Ppant_offsetSNR;
%                 storeacrossPpant_offsetRTs(icounter,:)=Ppant_offsetRTs;
%                 storeacrossPpant_offsetTOPO(icounter,:)=Ppant_offsetTOPO;
                
                
                
                
                
                icounter=icounter+1;
            end
            
            %save appropriately
            cd(basefol)
            cd('EEG')
            cd('GFX_Pre-RESS')
            
            switch ihz
                case 1
                    savename='GFX_Catchperformance_withSNR_15';
                case 2
                    savename='GFX_Catchperformance_withSNR_20';
                    case 3
                    savename='GFX_Catchperformance_withSNR_30';
                    case 4
                    savename='GFX_Catchperformance_withSNR_40';
            end
            
            save(savename,...
                'storeacrossPpant_onsetBP_invisible',...
                'storeacrossPpant_onsetSNR_invisible','tgrm', '-append')
            
            
            
        end
    end
    
    
end


if job.ppantCATCH_topotime==1
 getelocs
    window=[-3 3];
    snrmethod=2;
    srate=250;
    rmvbase=0;
    epochdur = sum(abs(window))*srate;
    
    timeid = [0:1/srate:epochdur];
    timeid= timeid-3;
    
    onsetc = ceil(epochdur)/2;
    % peakfreqsare=[20,40]; %hz
    %timing
    tt = 0:1/srate:60;
    
    %%
    param_spctrm.tapers = [1 1];
    param_spctrm.Fs= [250];
    param_spctrm.Fpass= [0 50];
    param_spctrm.trialave=0;
     peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; % don't change!        

    %%
    for ifol = allppants
        
        for hzis=[1,2,3,4,7];
            usehz = peakfreqsare(hzis);
            
            
            
            icounter=1;            
            cd(basefol)
            cd('EEG')
            cd(pdirs(ifol).name)
            
            
            load(['ppant_Catch_Epoched'])
            
            for itimezero = 1:2
                if itimezero==1
                    
                    %append them all. TARG-> Disappearing (more buttons
                    %pressed).
                    datatouse = ppant_SNREEG_disapBPwithincatch;
                   
                    ctype = 'report catch onset';
%                     bsrem = [-3 -1]; %seconds
                    
                else
                    datatouse = ppant_SNREEG_reapBPaftercatch;
%                     datatouse = cat(1, ppant_SNREEG_PFI_1_0,ppant_SNREEG_PFI_2_1,ppant_SNREEG_PFI_3_2,ppant_SNREEG_PFI_4_3);
%                     RTstouse = [durs1_0'; durs2_1'; durs3_2'];
%                     BPstouse = cat(1, BPs1_0, BPs2_1, BPs3_2);
                    ctype = 'report catch offset';
                    
                    
                    bsrem = [1 3]; %seconds
                    
                    
                end
                
                %             %rmvbaseline from EEG.
                datais=datatouse;
                
                
                % remove any NANs
                chNan = squeeze(mean(datais,2));
                chNan = squeeze(mean(chNan,2));
                nanid=find(isnan(chNan));
                
                keeptrials= setdiff(1:48, nanid);
                
                datais = datais(keeptrials,:,:);
                
                bsdata = zeros(size(datais));
                
                for ichan=1:64
                for itrial = 1:size(datais,1)
                    td = detrend(squeeze(datais(itrial,ichan,:)), 'linear');
%                     tdrm= mean(td(1,1:250));
%                     rmb= repmat(tdrm, [1 length(td)]);
                    bsdata(itrial,ichan,:) = td;%-rmb;
                end
                end
                datais=bsdata;
                %%
                snr_grmout=zeros(64,14); % chans by time points.
%                 snr_grmout=zeros(64,size(datais,1),14);
                for ichan=1:64
                    datacr=squeeze(datais(:,ichan,:));
                    
                [sgrm ,tgrm, fgrm] = mtspecgramc(datacr', movingwin, param_spcgrm);
                %%
                %conv SNR
                 tmps=sgrm;
                            %adjust for HBW 
                                k = ((param_spcgrm.tapers(1,1)));%
                                hbw = (k+1)./ (2.*[movingwin(1,1)]);
                                neighb = 2; %hz
                                
                                %space out kernel appropriately
                                distz = dsearchn(fgrm', [hbw, neighb, neighb*2+hbw]');
                                
                                % we don't want odd numbers, messes with
                                % calculations, so round allup to even.
                                ods = find(mod(distz,2)~=0);
                                %adjust.
                                distz(ods)= distz(ods)+1;
                                
                                %make kernel
                                kernelneighb = -1*(repmat(1/(distz(2)*2), [1, distz(2)]));
                                kernelskip = zeros(1,distz(1));
                                
                                kernelw= [kernelneighb, kernelskip, 1, kernelskip, kernelneighb];
                                if sum(kernelw)~=0
                                    error('')
                                end
                                 snr_sgrm= zeros(size(sgrm));
                            if snrmethod==1%method 1, using division.
                                
                               
                                
                                
                                
                                % loop over frequencies and compute SNR
                                numbins = distz(2);
                                skipbins = distz(1);
                                snr_sgrm= zeros(size(sgrm));
                                for itrial = 1:size(sgrm,3)
                                for itime=1:size(tmps,1)
                                    for hzi=numbins+1:length(fgrm)-numbins-1
                                        numer = tmps(itime,hzi, itrial);
                                        denom = nanmean( tmps(itime,[hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
                                        snr_sgrm(itime,hzi,itrial) = numer./denom;
                                        
                                    end
                                end
                                end
                                
                            else
                                %                             %conv SNR
                                %                             snr_sgrm =zeros(size(sgrm));
                                %
                                
                                %%
%                                 kernelw= [-1/8 -1/8 -1/8 -1/8 0 0 1 0 0 -1/8 -1/8 -1/8 -1/8];
                                
                                %snr on trials or
                                for itrial=1:size(sgrm,3)
                                    %compute SNR
                                    tmps= squeeze(sgrm(:,:,itrial));
                                    
                                    for itime= 1:size(tmps,1)
                                        checkput = conv(log(tmps(itime,:)), kernelw,'same');
                                        if ~isreal(checkput)
                                            snr_sgrm(itime,:,itrial)= nan(1, size(tmps,2));
                                        else
                                            snr_sgrm(itime,:,itrial)= conv(log(tmps(itime,:)), kernelw,'same');
                                        end
                                    end
                                end
                            
                             %% sanity check:
%                               figure(1);
%                               %plot specgrm and spctrm.
%                               subplot(221)
%                             imagesc(tgrm, fgrm, squeeze(mean(log(sgrm),3))'); colorbar; caxis([0 2])
%                             subplot(223)
%                             mtrials=squeeze(mean(log(sgrm),3));
%                             plot(fgrm, squeeze(mean(mtrials,1)));
%                             %plot specgrm and spctrm.
%                               subplot(222)
%                             imagesc(tgrm, fgrm, squeeze(mean((snr_sgrm),3))'); colorbar; caxis([0 2])
%                             subplot(224)
%                             mtrials=squeeze(mean((snr_sgrm),3));
%                             plot(fgrm, squeeze(mean(mtrials,1)));
%                             hold on; 
%                                % center kernel at 20 Hz
%     hlfkern= (length(kernelw)-1)/2;
%     [~,hzid]= min(abs(fgrm-20)); 
%     fvector = (hzid- hlfkern):(hzid+hlfkern);
%     plot(fgrm(fvector), kernelw*2, 'r')
                            end
                
                %store just stim freq.
                [~, hzid]= min(abs(fgrm-usehz));
                %%
                %reduce size.
                snrgrm20=squeeze(snr_sgrm(:,hzid,:));
                tbase= tgrm-3;
                %%
                if rmvbase==1
                    %rmv baseline ?
                    acrSNRsort_rmv=zeros(size(snrgrm20));
                    for itrial=1:size(snrgrm20,2);
                        
                        tmp=snrgrm20(:,itrial)';
                        %which baseline?
                        tidx= dsearchn(tbase', [bsrem]');
                        
                        %rmvbase
                        bs= mean(tmp(1,tidx(1):tidx(2)));
                        bs=repmat(bs, 1, length(tmp));
                        acrSNRsort_rmv(:,itrial)= tmp-bs;
                        %                         acrSNRsort_rmv(itrial,:)= tmp./bs;
                        
                    end
                    snrgrm20=acrSNRsort_rmv;
                end
                
                snr_grmout(ichan,:)=squeeze(nanmean(snrgrm20,2));
%                 snr_grmout(ichan,:,:)=snrgrm20';
                end
                %%
                
                switch itimezero
                    case 1 %store for across ppant plots:
                        
                        Ppant_catchBPonsetSNR_allchan=snr_grmout; 
                        
                    case 2
                        
                        
                        Ppant_catchBPoffsetSNR_allchan=snr_grmout; 
                        
                end
            end
            
            
            %%
            
                    savename=['Catchperformance_withSNR_' num2str(usehz)];
            
            
            
            save(savename,...
                'Ppant_catchBPoffsetSNR_allchan','Ppant_catchBPonsetSNR_allchan','-append');
                
            disp([ 'fin all chan append for ' num2str(ifol)])
        end
    end
end


if job.concatTOPOTIMEacrossppants==1
    %%
    for ihz=[1:4,7]
        
        usehz=peakfreqsare(ihz);
        loadname=['Catchperformance_withSNR_' num2str(usehz)];

        storeacrossPpant_onsetSNR_chans=zeros(length(allppants),64,14);
        storeacrossPpant_offsetSNR_chans=zeros(length(allppants),64,14);
        
        icounter=1;
        
        for ippant = allppants
            cd(basefol)
            cd('EEG')
            cd(pdirs(ippant).name);
            %onset types
            
            load(loadname)
%             storeacrossPpant_onsetSNR_chans(icounter,:,:)=squeeze(nanmean(Ppant_onsetSNR_allchan,2));
%             storeacrossPpant_offsetSNR_chans(icounter,:,:)=squeeze(nanmean(Ppant_offsetSNR_allchan,2));
            storeacrossPpant_onsetSNR_chans(icounter,:,:)=Ppant_catchBPonsetSNR_allchan;
            storeacrossPpant_offsetSNR_chans(icounter,:,:)=Ppant_catchBPoffsetSNR_allchan;
            
            icounter=icounter+1;
        end
        %%
        cd(basefol)
        cd('EEG')
        cd('GFX_Pre-RESS')
        %%
        
        savename=['GFX_Catchperformance_withSNR_' num2str(usehz) '_allchan'];
        
        save(savename, 'storeacrossPpant_onsetSNR_chans',...
            'storeacrossPpant_offsetSNR_chans', 'tgrm');
        
    end
    
end
% %%%%% %%% %% % % %
% %%%%% %%% %% % % %% %%%%% %%% %% % % %
% %%%%% %%% %% % % %   TOPO time plotted in a sep script for Will's data).
% %%%%% %%% %% % % %
% %%%%% %%% %% % % %
% %%%%% %%% %% % % %

%%
if job.erpimageacrossppants==1
    %% now plot across ppants.
    %reshape manually.
    
    %%
    getelocs
    rmvbase=0;
    
    
    plotinvisible=0; % if plotting invisible catch onset =1;
    for alignment=1%:2%1:2 % 1 for physical change at time zero, 2 for BP noticing catch change at time zero.
        clf
        
        for hzis=1:2
            cd([basefol filesep 'EEG'])
            cd('GFX_Pre-RESS')
            
            switch hzis
                case 1
                    if alignment==1
                        load('GFX_Catchperformance_withSNR_15')
                    else
                        load('GFX_Catchperformance_withSNR_15,BPaligned')
                    end
                    clims=[-.2 .2];
                    usehz=15;
                    
                case 2
                    if alignment==1
                        load('GFX_Catchperformance_withSNR_20')
                    else
                        load('GFX_Catchperformance_withSNR_20,BPaligned')
                    end
                    clims=[1 1.2];
                    usehz=20;
                case 3
                    if alignment==1
                        load('GFX_Catchperformance_withSNR_30')
                    else
                        load('GFX_Catchperformance_withSNR_30,BPaligned')
                    end
                    
                    usehz=30;
                    
                case 4
                    if alignment==1
                        load('GFX_Catchperformance_withSNR_40')
                    else
                        load('GFX_Catchperformance_withSNR_40,BPaligned')
                    end
                    usehz=40;
                    
                    
            end
            if hzis==3
                clf
            end
            tbase = tgrm-3;
            %%
            for itimezero=2%:2
                
                
                switch itimezero
                    case 1
                        if plotinvisible==1
                        useBP=storeacrossPpant_onsetBP_invisible;
                        useSNR=storeacrossPpant_onsetSNR_invisible;
                        
                        else
                        useBP=storeacrossPpant_onsetBP;
                        useSNR=storeacrossPpant_onsetSNR;
                        useRTs=storeacrossPpant_onsetRTs;
                        useTOPO=storeacrossPpant_onsetTOPO;
                        end
                        if alignment==1
                            ctype= 'onset';
                            chtype='catch onset';
                        else
                            ctype= 'onset, BPaligned';
                            chtype='reporting catch onset';
                        end
                        bsrem = [-3 -2]; %seconds
                    case 2
                        
                        useBP=storeacrossPpant_offsetBP;
                        useSNR=storeacrossPpant_offsetSNR;
                        useRTs=storeacrossPpant_offsetRTs;
                        useTOPO=storeacrossPpant_offsetTOPO;
                        
                        if alignment==1
                            ctype= 'offset';
                            chtype= 'catch offset';
                        else
                            ctype= 'offset, BPaligned';
                            chtype='reporting catch offset';
                        end
                        
                        bsrem = [2 3]; %seconds
                        
                end
                
                
                
                %                 for ip= 1:length(allppants)
                %                     trialind = [1:ntrials] + ntrials*(ip-1);
                %                     acrBP(trialind,:,:) = useBP(ip,:,:);
                %                     acrSNR(trialind,:,:)=    useSNR(ip,:,:);
                %                     acrRTs(trialind)=    useRTs(ip,:);
                %
                %                 end
                
                %take mean across ppants for erp images..
                acrBP= squeeze(nanmean(useBP,1));
                
                acrSNR=squeeze(nanmean(useSNR,1));
                acrTOPO=squeeze(nanmean(useTOPO,1));
                
                
                
                %sort across all
                %sort trial index by longest RT first.
                %%
                figure(1)
%                 clf
%% not topoplotting anymore.                
%                 usechan=30;
%                 subplot(2,2,1:2)
%                 topoplot(acrTOPO, elocs(1:64), 'emarker2', {usechan, 'o', 'w'});
%                 c=colorbar;
%                 caxis([0 15])
%                 ylabel(c, {['10log10(SNR)'];['-3000:-100ms']})
%                 hold on
%                 title({[num2str(usehz) ' Hz SSVEP']})
%                 set(gca, 'fontsize', 20)
                %%
                subplot(3,1,1)
                %         [sortedRTs, cid] = sort(acrRTs, 'descend');
                
                %
                imagesc([-3:1/60:3], 1:size(acrBP,1), acrBP);%(cid,:));
                %
                title({['Buttons Pressed'] }, 'fontsize', 25)
                c=colorbar
                ylabel(c, 'Total')
                %                 ylabel('resampled catch trials')
                xlim([-2.5 2.5])
                set(gca, 'fontsize', 25, 'ytick', [])
                hold on
                plot([0 0] , ylim, ['w:'], 'linewidth', 3)
                plot([0 0] , ylim, ['k:'], 'linewidth', 2)
                %                 subplot(4,2,7)
                
                ppantMeanBP= squeeze(mean(useBP,2));
                %adjust standard error as per COusineau(2005)
                %confidence interval for within subj designs.
                % y = x - mXsub + mXGroup,
                x = ppantMeanBP;
                
                mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
                mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
                
                %for each observation, subjtract the subj average, add
                %the group average.
                NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
                
                %compute new stErr %which version?
                %             stE = std(NEWdata)/sqrt(size(x,1));
                %                 shadedErrorBar([-3:1/60:3], mean(ppantMeanBP,1),stE)
                %                 set(gca, 'fontsize', 15)
                %                 hold on;
                                ylabel({['Normalized'];['trial count']})
                %                 plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                %%
                xlabel(['Time from ' chtype])
                
                %SNR
                try
                subplot(3, 1, hzis+1);
                catch
                    subplot(3, 1, hzis-1);
                end
                
                %reorder SNR
                acrSNRsort=acrSNR;
                
                %%
                %                 %smooth across trials
                sm_snrgrm20=zeros(size(acrSNRsort));
                
                for itime=1:size(acrSNRsort,2)
                    sm_snrgrm20(:,itime)= smooth(acrSNRsort(:,itime),15);
                end
                
                
                %     imagesc(acrSNR);
                imagesc(tgrm-3, 1:size(acrSNRsort,1), acrSNRsort);
                c=colorbar;
                ylabel(c, {['log(SNR)']})
                %
%                 caxis([-1*max(max(sm_snrgrm20))/3 max(max(sm_snrgrm20))/3])
                caxis([clims])
%                 caxis([-1 1])
%                 caxis([-5 5])
                %         caxis([0 10])
                title([num2str(usehz) 'Hz dyn-SSVEP SNR'])
                %             set(gca, 'ytick', 1:ntrials, 'yticklabel', round(sortedRTs./60,2))
                hold on
                plot([0 0] , ylim, ['k:'], 'linewidth', 3)
                plot([0 0] , ylim, ['w:'], 'linewidth', 2)
                xlim([-2.5 2.5])
                set(gca, 'fontsize', 25, 'ytick', [])
                %                 subplot(4,2,8)
                %
                %plot across ppant trace
                ppantMeanSNR= squeeze(mean(useSNR,2));
                
                %adjust standard error as per COusineau(2005)
                %confidence interval for within subj designs.
                % y = x - mXsub + mXGroup,
%                 x = ppantMeanSNR;
%                 
%                 mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
%                 mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
%                 
%                 %for each observation, subjtract the subj average, add
%                 %the group average.
%                 NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
%                 
                %             %compute new stErr %which version?
                %             stE = std(NEWdata)/sqrt(size(x,1));
                %
                %                 shadedErrorBar(tgrm-3, mean(ppantMeanSNR,1),stE,'k',[])
                %                 set(gca, 'fontsize', 15)
                %                 hold on;
                %                 %%
                %
                %                 ylabel('mean SNR')
                %                 axis tight
                %                 plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                xlabel(['Time from ' chtype])
                ylabel({['Normalized'];['trial count']})
                set(gca, 'fontsize', 25)
                cd(basefol)
                cd('Figures')
                set(gcf, 'color', 'w')
                shg
                %%
%                 if itimezero==1 && usehz==40 && alignment==1
%                     yt=ylim;
%                     xV = [-.5 -.5 .5 .5]; %x coordinates of vertices.
%                     yV = [yt(1) yt(2) yt(2) yt(1)]; %y coordinates of vertices.
%                     pch=patch(xV,yV,[.5 .5 .5]);
%                     pch.FaceAlpha =(.75);
%                     shg
%                 end
                %%
                
%                 figure(2)
%                 for ippant=1:length(allppants);
%                    subplot(4,4,ippant)
%                    imagesc(squeeze(useSNR(ippant,:,:)))
%                    shg
%                     
%                 end
            end
            %%
        colormap('viridis')
        %%
        end
        cd(basefol)
        cd('Figures')
        cd('SNR trial-by-trial')
        if hzis==2
           print('-dpng', ['Catch SNR summary,'  chtype ' 1stharms.png'])
        elseif hzis==4
            print('-dpng', ['Catch SNR summary,'  chtype ' 2ndharms.png'])
        end
             
    end
end


if job.BPandSSVEPtimecourseacrossppants==1;
    %%
    getelocs
    rmvbase=0;
    checksigON=0;
    
    interpCatch=0; % to plot with interpolated region for catch onsets
    
    
    for alignment=2%1:2%1:2 % 1 for physical change at time zero, 2 for BP noticing catch change at time zero.
        figure(1);
        clf
        plcount=1;
        legendprint=[];
        legendnames=[];
        lc=1;
        for hzis=3:4
            cd(basefol)
            
            cd('EEG')
            
            switch hzis
                case 1
                    if alignment==1
                        load('GFX_Catchperformance_withSNR_15')
                    else
                        load('GFX_Catchperformance_withSNR_15,BPaligned')
                    end
                    lint='-';
                    usehz=15;
                    sigheight= 4.1; %where on figure to place '*' if significant.
                case 2
                    if alignment==1
                        load('GFX_Catchperformance_withSNR_20')
                    else
                        load('GFX_Catchperformance_withSNR_20,BPaligned')
                    end
                    usehz=20;
                    lint=':';
                    sigheight= 3.9; %where on figure to place '*' if significant.
                case 3
                    if alignment==1
                        load('GFX_Catchperformance_withSNR_30')
                    else
                        load('GFX_Catchperformance_withSNR_30,BPaligned')
                    end
                    usehz=30;
                    lint='-';
                    sigheight= 3.9; %where on figure to place '*' if significant.
                case 4
                    if alignment==1
                        load('GFX_Catchperformance_withSNR_40')
                    else
                        load('GFX_Catchperformance_withSNR_40,BPaligned')
                    end
                    usehz=40;
                    lint=':';
                    sigheight= 3.9; %where on figure to place '*' if significant.
            end
            tbase = tgrm-3;
            %%
            %             clf
            ttestdata=[];
            for itimezero=1:2%1%:2
                
                
                switch itimezero
                    case 1
                        useBP=storeacrossPpant_onsetBP;
                        useSNR=storeacrossPpant_onsetSNR;
                        useRTs=storeacrossPpant_onsetRTs;
                        useTOPO=storeacrossPpant_onsetTOPO;
                        
                        if alignment==1
                            ctype= 'onset';
                            chtype='catch start-end';
                            xlabelis = {['Time from catch onset [s]']};
                        else
                            ctype= 'onset, BPaligned';
                            xlabelis = {['Time from reporting'];['catch onset [s]']};
                            
                        end
                        bsrem = [-3 -2]; %seconds
                        col=[0 .5 0]; %dark green.
                    case 2
                        
                        useBP=storeacrossPpant_offsetBP;
                        useSNR=storeacrossPpant_offsetSNR;
                        useRTs=storeacrossPpant_offsetRTs;
                        useTOPO=storeacrossPpant_offsetTOPO;
                        
                        if alignment==1
                            ctype= 'offset';
                            xlabelis = {['Time from catch offset [s]']};
                            
                        else
                            ctype= 'offset, BPaligned';
                            xlabelis = {['Time from reporting'];['catch offset [s]']};
                            chtype='catch report';
                        end
                        
                        bsrem = [2 3]; %seconds
                        col='r';
                        
                        
                end
                
                
                %%
                figure(1)
                
                for ip=1:2
                    switch ip
                        case 1
                            used=useBP;
                            col='b';
                            lint='-';
                            ylimsare = [0 2.5];
                            basep=1;
                        case 2
                            used=useSNR;
                            basep=3;
                            if hzis==1 ||hzis==3
                                col='r';
                            elseif hzis==2 || hzis==4
                                col=[.5 .5 .5];
                            end
                            
                            
                            ylimsare= [ -.1 1.2];
                    end
                    placeme = basep+ (1*itimezero-1);
                    
                    subplot(2,2,placeme)
                    
                    
                    
                    
                    
                    
                    if rmvbase==1
                        used=used- nanmean(used(:));
                    end
                    
                    
                    
                    if interpCatch==1 && itimezero==1 && ip==2 %&& hzis==2
                        %may need to interpolate at trial level. see what this looks like.
                        %time vector
                        timeid=tgrm;
                        tvector=timeid-3;
                        points= 1:length(tvector);
                        if alignment==2
                        %remove bad section from trace.
                        %we used a one second sliding window, so:
                        badsec = dsearchn(tvector', [-1.35 0]');
                        else
                        badsec = dsearchn(tvector', [-.4 .4]');
                            
                        end
                            
                        tmp=points;
                        tmp(badsec(1):badsec(2))=[];
                        %input x vector with points missing.
                        xIN=tmp;
                        plotd=squeeze(mean(used,2));
%                       
                       %use in plots.
                        ppantMeanSNR=plotd;
                    else 
                         
                    ppantMeanSNR = squeeze(nanmean(used,2));
                    end
                    
                    
                    
                    
                    hold on
                   
                    %adjust standard error as per COusineau(2005)
                    %confidence interval for within subj designs.
                    % y = x - mXsub + mXGroup,
                    
                    x = ppantMeanSNR;
                    if ip~=1
                    mXppant =squeeze( nanmean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
                    mXgroup = mean(nanmean(mXppant)); %mean overall (remove b/w sub differences
                    
                    %for each observation, subjtract the subj average, add
                    %the group average.
                    NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
                    
                    else
                        NEWdata=ppantMeanSNR;
                    end
                    
                    %compute new stErr %which version?
                    stE = std(NEWdata)/sqrt(size(x,1));
                    
                    if ip==1
                        sh=shadedErrorBar(-3:1/60:3, mean(ppantMeanSNR,1),stE,[lint],[1]);
                        ylabel(['Buttons pressed'])
                    else
                        sh=shadedErrorBar(tgrm-3, mean(ppantMeanSNR,1),stE,[lint],[1]);
                        ylabel({['log(SNR)']})
                    end
                    sh.mainLine.LineWidth=3;
                    
                    sh.mainLine.Color = col;
                    sh.patch.FaceColor = col;
                    sh.edge(1).Color = col;
                    sh.edge(2).Color = col;
                    
                    %                 if hzis==2
                    %                     sh.mainLine.Color= 'k';
                    %                 end
                    set(gca, 'fontsize', 15)
                    hold on;
                    %                 %%
                    %
                    
                    
                    
                    
                    %             title({[num2str(usehz) ' Hz SSVEP']})
                    xlabel(xlabelis)
                    %             xlabel('Time from perceptual report')
                    set(gca, 'fontsize', 25)
                    xlim([-3 3])
                    ylim( [ylimsare])
                    
                    
                    
                    set(gcf, 'color', 'w')
                    
                    if ip==2
                        %             ttestdata(plcount,:,:) = ppantMeanSNR;
                        legendprint(lc)=sh.mainLine;
                        legendnames = [legendnames, {[num2str(usehz) 'Hz']}];
                        lc=lc+1;
                    end
                    
                    plcount=plcount+1;
                    
                    
                    %place legend
                    if (hzis==2||hzis==4) && ip==2
                        if itimezero==2
                            lg=legend([legendprint(2) legendprint(4)], {legendnames{1}, legendnames{3}});
set(lg, 'location', 'NorthEast')
                        end
                        
                        
                        
                    end
                    
                end
                
                        %plot patch:
            axis tight
            end
        
        %%
        if checksigON==1
        %check for sig
        pvals=zeros(1,size(ttestdata,3));
        tvals=zeros(1,size(ttestdata,3));
        for itime = 1:size(ttestdata,3)
            
            try [h,pvals(itime),~,stat]=ttest(ttestdata(1,:,itime), ttestdata(2,:,itime));
                shuffType=1;
            catch
                [h,pvals(itime),~,stat]=ttest(ttestdata(plcount-1,:,itime)); %compares to zero.
                shuffType=2; %whether or not to skip the non-parametric test for sig.
            end
            
            tvals(itime)= stat.tstat;
        end
        %%
        %             q=fdr(pvals,.05);
        
            sigs=find(pvals<.05);
            
            %perform cluster based correction.
            if length(sigs)>2
                % find biggest cluster:
                %finds adjacent time points
                vect1 = diff(sigs);
                v1 = (vect1(:)==1);
                d = diff(v1);
                clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
                %grab largest
                %                 ignore bad points.
                if hzis==2 &&itimezero==1
                    pvals(11:22)=nan;
                end
                
                % find biggest cluster:
                %finds adjacent time points
                sigs = find(pvals<.05);
                
                vect1 = diff(sigs);
                v1 = (vect1(:)==1);
                d = diff(v1);
                clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
                [~,maxClust] = max(clusterSTandEND(:,2)-clusterSTandEND(:,1));
                if itimezero==1 && hzis==2 && alignment==2
                    maxClust=3;
                end
                %%
                for icl=1%:size(clusterSTandEND,1)
                    
                    %start and end are now:
                    STC=sigs(clusterSTandEND(maxClust,1));
                    ENDC=sigs(clusterSTandEND(maxClust,2)+1);
                    checktimes =STC:ENDC;
                    observedCV = sum(abs(tvals(checktimes)));
                    % now shuffle condition labels to see if this cluster is
                    % sig (compared to chance).
                    sumTestStatsShuff = zeros(1,2000);
                    
                    for irand = 1:2000
                        %testing the null that it isn't mismatched - matched at time 2
                        % which creates a diff. so select from either!
                        shD= zeros(2,length(checktimes),size(ttestdata,2));
                        
                        %change shuffle parameters based on test of
                        %interest. (ie between conditions, or temporal
                        %null).
                        if shuffType==1 %null is that no condition differences.
                            for ipartition = 1:2
                                for ippant = 1:21
                                    for itime=1:length(checktimes)
                                        
                                        if mod(randi(100),2)==0 %if random even number
                                            pdata = ttestdata(1,randi(21), checktimes(itime)); %select both chans
                                        else
                                            pdata = ttestdata(2,randi(21), checktimes(itime));
                                        end
                                        
                                        shD(ipartition,itime,ippant) = pdata;
                                    end
                                end
                            end
                        else %null is that there are no temporal coincident sig values.
                            for ipartition = 1:2
                                for ippant = 1:21
                                    for itime=1:length(checktimes)
                                        
                                        %take random timepoint.
                                        pdata = ttestdata(1,ippant, randi(size(ttestdata,3)));
                                        
                                        
                                        shD(ipartition,itime,ippant) = pdata;
                                    end
                                end
                            end
                        end
                        %now compute difference between out hypothetical topoplots,
                        % and test for sig, checking the accumulated test statistic at our
                        % times of interest
                        tvalspertimepoint = zeros(1,length(checktimes));
                        
                        testdata = squeeze(shD(1,:,:)) - squeeze(shD(2,:,:));
                        
                        for itest = 1:length(checktimes) %test each time point
                            
                            [~, p, ~,stat]= ttest(testdata(itest,:));
                            
                            tvalspertimepoint(1,itest) = stat.tstat;
                        end
                        
                        sumTestStatsShuff(1,irand) = sum(abs(tvalspertimepoint));
                    end %repeat nshuff times
                    
                    
                    %is the observed greater than CV?
                    % plot histogram:
                    figure(2);
                    if plcount==2
                        clf
                    end
                    subplot(2,1, plcount-1)
                    H=histogram(abs(sort(sumTestStatsShuff)));
                    % fit CDF
                    cdf= cumsum(H.Data)/ sum(H.Data);
                    %the X values (actual CV) corresponding to .01
                    [~,cv05uncorr] = (min(abs(cdf-.95)));
                    [~,cv01uncorr] = (min(abs(cdf-.99)));
                    [~,cv001uncorr] = (min(abs(cdf-.999)));
                    hold on
                    pCV=plot([observedCV observedCV], ylim, ['r-']);
                    
                    p05=plot([H.Data(cv05uncorr) H.Data(cv05uncorr)], ylim, ['k:']);
                    plot([H.Data(cv01uncorr) H.Data(cv01uncorr)], ylim, ['k:']);
                    plot([H.Data(cv001uncorr) H.Data(cv001uncorr)], ylim, ['k:']);
                    legend([pCV p05], {['observed'] ['p01'] })
                    
                    %%
                    if observedCV>H.Data(cv05uncorr)
                        title(['sum tvals = ' num2str(observedCV)]);
                        %              title('Spatial Cluster  Significant!')
                        for itime=checktimes
                            figure(1);
                            hold on
                            plot(tgrm(itime)-3, sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', sh.mainLine.Color)
                        end
                    end
                    
                    %
                end
            end
        end
        %%
        
        
        hold on
%         plot(xlim, [0 0], ['k:']);
%         plot([0 0] , ylim, ['k:'], 'linewidth', 1)

        
        shg
        
        hold on
        
        %%
%         axis tight
    end
   %%
   cd(basefol)
    cd('Figures')
    cd('SNR-timecourse')
    %%
    print('-dpng', ['Catch trace Bground SSVEP summary, during ' num2str(chtype) ' 2ndharm.png'])
end
end% end