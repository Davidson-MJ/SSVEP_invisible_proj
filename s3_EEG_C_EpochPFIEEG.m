% s3_CD_EpochPFIEEG
% clear all

%% UPDATED for Wills data:
try cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA')
catch
end
basefol=pwd;
cd('EEG')
dirs = dir([pwd filesep '*_*' 'EEG']);

%Epoch by number of targets involved, and direction (disapp and reapp).
% also stores duration per event.
% also stores cumularive BP (single array) per event, for next job plots.

dbstop if error

job.epochperppant_PFI=0; %epochs raw EEG and labels according to number and direction of target PFI involved (still time-domain).

%% previous (archived, extra jobs and sanity checks.).




job.erpimagePFIusingSNR=0; %take SNR, in freq domain within participant, any number of targets, just increasing or decreasing PFI.

job.ppantPFI_topotime=0; %appends to above, for topos of PFI increase or decrease. % completed for 15,20,5 Hz.

job.concaterpdataacrossppants=0; %Concats all Increase/ Decrease  % change here for MIMINIM PFI duration required in plot..
job.erpimageacrossppants=0;
job.concaterpdataacrossppants_keepnumberSeparate=0; %concats for sep num PFI, based on accum period after.

job.concatTOPOTIMEacrossppants=0;

% plotting TOPO TIME in new script for Will's data.
% job.plotTOPOtimeacrossppants=0; %%%%%%%%%%% use this to plot spatial correlation over time also (see job2 within).

%%%%%%% plot:
job.BPandSSVEPtimecourseacrossppants=0; %this uses all/any PFI. (any duration, any number).
%%%%%%%
job.SSVEPtimecourseSelectchans=0; %combines across certain channels.
job.FirstSIGintimecourse=0; %this uses all/any PFI. (any duration, any number).
job.BPandSSVEPtimecourseacrossppants_numSeparate=0; %this uses all/any PFI. (separating by number).
job.BPandSSVEPtimecourseacrossppants_dursSeparate=0; %this uses all/any PFI. (separating by duration, any number).

%remaining participants after behavioral data analysis/exclusion
% allppants=[1,2,4,6,9:16,18]; %
allppants=[1,2,4,6,7,9:19]; % new


param_spcgrm.tapers = [1 1];
param_spcgrm.Fs= [250];
param_spcgrm.Fpass= [0 50];
param_spcgrm.trialave=0;
param_spcgrm.pad = 2;
movingwin=[2.5,.25]; %

peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; % don't change!


if job.epochperppant_PFI==1
    cd(basefol)
    cd('Behaviour')
    
    %% load data to determine physical catch timing.
    load('MD_AllBP_perppant')
    load('PFI_data')
    %%
    
    window=[-3 3];
    
    srate=250;
    
    
    epochdur = sum(abs(window))*srate;
    onsetc = ceil(epochdur)/2;
    %timing
    tt = [0:1/srate:60];
    
    
    
    for ifol =allppants
        cd(basefol)
        cd('EEG')
        cd(dirs(ifol).name)
        %load new rereferenced data.
        load(['P' num2str(ifol) '_autopreprocd'])
        
        %separate by number involved
        %disappearances:
        ppant_SNREEG_PFI_0_1 = [];
        counter0_1=1;
        ppant_SNREEG_PFI_1_2 = [];
        counter1_2=1;
        ppant_SNREEG_PFI_2_3 = [];
        counter2_3=1;
        ppant_SNREEG_PFI_3_4 = [];
        counter3_4=1;
        %Reappearances
        ppant_SNREEG_PFI_4_3 = [];
        counter4_3=1;
        ppant_SNREEG_PFI_3_2 = [];
        counter3_2=1;
        ppant_SNREEG_PFI_2_1 = [];
        counter2_1=1;
        ppant_SNREEG_PFI_1_0 = [];
        counter1_0=1;
        
        
        durs0_1=0;
        durs1_0=0;
        durs1_2=0;
        durs2_1=0;
        durs2_3=0;
        durs3_2=0;
        durs3_4=0;
        durs4_3=0;
        
        BPs0_1=[];
        BPs1_0=[];
        BPs1_2=[];
        BPs2_1=[];
        BPs2_3=[];
        BPs3_2=[];
        BPs3_4=[];
        BPs4_3=[];
        
        
        Locationwas=[]; %same for location of PFI events BPs.
        for itrial=1:48
            
            
            
            SNRdata = squeeze(EEG_preprocd(:,:,itrial)); %all channels, 20 and 40 hz.
            trialdata= PFI_only(ifol).Trial(itrial);
            %first establish the direction of changes (disap /reappear
            %targets).
            accumBP = nansum(trialdata.allBPs,1);
            
            
            % want to remove transients?
            excludeTransientlength=15; % frames = 1/60
            if excludeTransientlength>0;%; %ms
                % then adjust for short transient button presses:
                
                %                 figure(1);
                %                 plot(accumBP, 'k'); %plot original.
                
                % find all diffs
                sws = find(diff(accumBP));
                
                %how long between BPs?
                diff_sws = diff(sws);
                
                rmv = find(diff_sws<excludeTransientlength);
                % now work through and replace.
                for irm= rmv
                    
                    tremovesw =sws(irm);
                    %replace values for this period:
                    vals = accumBP(1,tremovesw);
                    %replace to = sws
                    trep = sws(irm+1);
                    %replace now:
                    accumBP(1,tremovesw:trep)=vals;
                    
                    
                end
                
                %                 hold on;
                %                 plot(accumBP, 'r'); %check method.
                
            end
            
            
            
            
            
            directionBP = diff(accumBP);
            directionBPtimes = find(directionBP~=0);
            directionBPtimes=directionBPtimes+1; %account for diff function
            %now remove zeros, leaving time of all disap/reap
            directionBP(directionBP==0)=[];
            
            
            
            % present visualization as sanity check.
            %             figure(1); clf;
            %             plot(accumBP,'k', 'linew', 3); title(['Trial ' num2str(itrial)]);
            %             ylim([-1 5])
            
            if trialdata.Goodtrial==1
                
                for iPFI = 1:length(directionBPtimes)
                    %Collect event info.
                    timeis = directionBPtimes(iPFI);
                    perceptis = accumBP(directionBPtimes(iPFI));
                    perceptwas = accumBP(directionBPtimes(iPFI)-1);
                    
                    
                    outgoingPFI=[];
                    %disap or reap.
                    if perceptwas<perceptis
                        PFIdir=1; %disappearing, an extra button is being pressed
                        col='b';
                    else
                        PFIdir=-1; %Reappearing
                        col='r';
                    end
                    
                    %gather EEG time index.
                    [~,timePFI]= min(abs(tt-timeis/60));
                   
                    
                    % ignore events at start/end of block.
                    
                    if  timePFI-onsetc>1 && timePFI+onsetc<size(SNRdata,2)
                        %collect data
                        outgoingPFI= SNRdata(:,timePFI-onsetc:timePFI+onsetc);
                        
                        
                        %store according to direction, and number of
                        %targets involved.
                        if PFIdir==1 %increasing targets absent
                            switch perceptis
                                case 1
                                    ppant_SNREEG_PFI_0_1(counter0_1,:,:) = outgoingPFI;
                                    %how long was this event?
                                    %find the correct framestamp and duration.
                                    %                                     [~,indx]=min(abs(trialdata.PFI_disap_1target_framestart - timeis));
                                    %                                     durs=trialdata.PFI_disap_1target_durs;
                                    %save duration of this event for averaging
                                    %'best' later.
                                    %                                     durs0_1(counter0_1) = durs(indx);
                                    %                                     %also store BP trace.
                                    BPs0_1(counter0_1,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                    counter0_1=counter0_1+1;
                                case 2
                                    ppant_SNREEG_PFI_1_2(counter1_2,:,:) = outgoingPFI;
                                    %                                     [~,indx]=min(abs(trialdata.PFI_disap_2target_framestart - timeis));
                                    %                                     durs=trialdata.PFI_disap_2target_durs;%
                                    
                                    %                                     durs1_2(counter1_2) = durs(indx);%
                                    %also store BP trace.
                                    BPs1_2(counter1_2,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                    counter1_2=counter1_2+1;
                                case 3
                                    ppant_SNREEG_PFI_2_3(counter2_3,:,:) = outgoingPFI;
                                    %                                     [~,indx]=min(abs(trialdata.PFI_disap_3target_framestart - timeis));
                                    %                                     durs=trialdata.PFI_disap_3target_durs;
                                    %                                     durs2_3(counter2_3) = durs(indx);
                                    BPs2_3(counter2_3,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                    counter2_3 = counter2_3+1;
                                    
                                case 4
                                    ppant_SNREEG_PFI_3_4(counter3_4,:,:) = outgoingPFI;
                                    
                                    %                                     [~,indx]=min(abs(trialdata.PFI_disap_4target_framestart - timeis));
                                    %                                     durs=trialdata.PFI_disap_4target_durs;
                                    %                                     %save duration of this event for averaging
                                    %                                     %'best' later.
                                    %                                     durs3_4(counter3_4) = durs(indx);
                                    BPs3_4(counter3_4,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                    counter3_4 = counter3_4+1;
                                    
                            end
                            
                            
                        else %decreasing number targets absent.
                            switch perceptwas
                                case 1
                                    ppant_SNREEG_PFI_1_0(counter1_0,:,:) = outgoingPFI;
                                    %how long was this event?
                                    %find the correct framestamp and duration.
                                    %                                         [~,indx]=min(abs(trialdata.PFI_disap_1target_framestart - timeis));
                                    %                                         durs=trialdata.PFI_disap_1target_durs;
                                    % store BP also
                                    BPs1_0(counter1_0,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                    
                                    
                                    %save duration of this event for averaging
                                    %'best' later.
                                    %                                         durs1_0(counter1_0) = durs(indx);
                                    
                                    counter1_0=counter1_0+1;
                                    
                                case 2
                                    ppant_SNREEG_PFI_2_1(counter2_1,:,:) = outgoingPFI;
                                    %how long was this event?
                                    %find the correct framestamp and duration.
                                    %                                         [~,indx]=min(abs(trialdata.PFI_disap_2target_framestart - timeis));
                                    %                                         durs=trialdata.PFI_disap_2target_durs;
                                    % store BP also
                                    BPs2_1(counter2_1,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                    
                                    %save duration of this event for averaging
                                    %'best' later.
                                    %                                         durs2_1(counter2_1) = durs(indx);
                                    
                                    counter2_1=counter2_1+1;
                                    
                                case 3
                                    ppant_SNREEG_PFI_3_2(counter3_2,:,:) = outgoingPFI;
                                    %                                         [~,indx]=min(abs(trialdata.PFI_disap_3target_framestart - timeis));
                                    %                                         durs=trialdata.PFI_disap_3target_durs;%
                                    %                                         %save duration of this event for averaging
                                    %                                         %'best' later.
                                    %                                         durs3_2(counter3_2) = durs(indx);
                                    BPs3_2(counter3_2,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                    
                                    counter3_2 = counter3_2+1;
                                case 4
                                    ppant_SNREEG_PFI_4_3(counter4_3,:,:) = outgoingPFI;
                                    %                                     [~,indx]=min(abs(trialdata.PFI_disap_4target_framestart - timeis));
                                    %                                     durs=trialdata.PFI_disap_4target_durs;
                                    
                                    %                                     %save duration of this event for averaging
                                    %                                     %'best' later.
                                    %                                     durs4_3(counter4_3) = durs(indx);
                                    BPs4_3(counter4_3,:)=accumBP(1,(timeis+window(1)*60):(timeis+window(2)*60));
                                    
                                    counter4_3 = counter4_3+1;
                                    
                            end %nPFI
                        end%perceptdir
                    end %if epochrange
                end %all events
            end %good trials
            
        end %all trials
        
        savename= ['ppant_PFI_Epoched'];
        
        
        save(savename, ...
            'ppant_SNREEG_PFI_0_1',...
            'ppant_SNREEG_PFI_1_0',...
            'ppant_SNREEG_PFI_1_2',...
            'ppant_SNREEG_PFI_2_1',...
            'ppant_SNREEG_PFI_2_3',...
            'ppant_SNREEG_PFI_3_2',...
            'ppant_SNREEG_PFI_3_4',...
            'ppant_SNREEG_PFI_4_3',...
            'BPs0_1',...
            'BPs1_0',...
            'BPs1_2',...
            'BPs2_1',...
            'BPs2_3',...
            'BPs3_2',...
            'BPs3_4',...
            'BPs4_3',...
            'window','tt', 'Locationwas')
        
        disp(['Finished ppant ' num2str(ifol)])
        
    end
end
%%


if job.erpimagePFIusingSNR==1
    getelocs
    window=[-3 3];
    
    srate=250;
    rmvbase=0;
    epochdur = sum(abs(window))*srate;
    
    timeid = [0:1/srate:epochdur];
    timeid= timeid-3;
    
    onsetc = ceil(epochdur)/2;
    
    %timing
    tt = 0:1/srate:60;
    
    snrmethod=2; %1 for MXC division, 2 for convolution. (2=much faster)
    %%
    param_spctrm.tapers = [1 1];
    param_spctrm.Fs= [250];
    param_spctrm.fpass= [0 50];
    param_spctrm.trialave=0;
    
    
    %%
    for ifol = allppants
        for hzis=[1,2,3,4,7,8]
            switch hzis
                case 1
                    usehz=15; %
                case 2
                    usehz=20;
                case 3
                    usehz=30;
                case 4
                    usehz=40;
                case 5
                    usehz=45;
                case 6
                    usehz=60;
                case 7
                    usehz= 5;
                case 8
                    usehz = 25;
                case 9
                    usehz= 35;
            end
            
            
            
            icounter=1;
            
            cd(basefol)
            cd('EEG')
            
            cd(dirs(ifol).name)
            %%
            load(['ppant_PFI_Epoched'])
            %%
            cd(basefol)
            cd('Figures')
            cd('Ppant rawSNR summary (PFI)')
            %%
            for itimezero = 1:2
                if itimezero==1
                    
                    %append them all. TARG-> Disappearing (more buttons
                    %pressed).
                    datatouse = cat(1, ppant_SNREEG_PFI_0_1,ppant_SNREEG_PFI_1_2,ppant_SNREEG_PFI_2_3,ppant_SNREEG_PFI_3_4);
                    %                     RTstouse = [durs0_1'; durs1_2'; durs2_3'];
                    BPstouse = cat(1, BPs0_1, BPs1_2, BPs2_3, BPs3_4);
                    ctype = 'PFI increase';
                    bsrem = [-3 -1]; %seconds
                    
                else
                    
                    datatouse = cat(1, ppant_SNREEG_PFI_1_0,ppant_SNREEG_PFI_2_1,ppant_SNREEG_PFI_3_2, ppant_SNREEG_PFI_4_3);
                    %                     RTstouse = [durs1_0'; durs2_1'; durs3_2'];
                    BPstouse = cat(1, BPs1_0, BPs2_1, BPs3_2, BPs4_3);
                    ctype = 'PFI decrease';
                    
                    
                    bsrem = [1 3]; %seconds
                    
                    
                end
                
                
                
                %plot the topo for sanity check: (pre disap)
                windcheck= [-3 -0.1];
                tidx=dsearchn(timeid', [windcheck]');
                allt=1:size(datatouse,1);
                
                %                 % Do later (future scripts) instead. for now store all.
                %                 if itimezero==2 % check the duration of disap was not a transient
                %                     %onset or offset window?
                %                     % find duration of pre zero bp. (and
                %                     % flip)
                %                     BPsshort= fliplr(BPstouse(:,1:180));
                %
                %                 else% disappearing target types
                %                     %
                %                     BPsshort = BPstouse(:,181:end);
                %                 end
                %                 %first difference in 2nd dimension.
                %                 tmp=(diff(BPsshort,1,2));
                %                 dursCheck=zeros(size(BPsshort,1),1);
                %                 for itrial=1:size(BPsshort,1)
                %                     try dursCheck(itrial,:)= find(squeeze(tmp(itrial,:)),1, 'first');
                %                     catch % the button press didn't end!
                %                         dursCheck(itrial,:) = 180; %ie at least epoch length.
                %                     end
                %                 end
                %
                %                 %retain which trials?
                %                 trialtypesperTargPresent = find(dursCheck>exclTransientBP);
                %
                %
                
                
                %reduce size.
                topod= squeeze(datatouse(:,:,tidx(1):tidx(2)));
                %calc statSNR per channel, over trials
                
                
                %%
                snr20=zeros(64,1);
                
                for ichan= 1:64
                    
                    [s,f]=mtspectrumc(squeeze(topod(:,ichan,:))', param_spctrm)    ;
                    
                    %comp SNR;
                    kernelw = [-.25 -.25 0 0 1 0 0 -.25 -.25];
                    
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
                usechan=62; %Oz for all
                %%
                figure(1)
                clf
                subplot(4,2,1:2)
                topoplot(snr20, elocs(1:64), 'emarker2', {usechan, 'o', 'w'});
                c=colorbar;
                caxis([0 2])
                ylabel(c, 'SNR')
                title({[num2str(usehz) 'Hz SNR '];['-3000:-100ms']})
                set(gca, 'fontsize', 15)
                %%
                
                %plot BP for comparison
                subplot(4,2, [3 5])
                %sort trial index by longest RT first.
                %                 [sortedRTs, cid] = sort(RTstouse(allt), 'descend');
                
                %sort by accum PFI in window after onset.
                if itimezero==1
                    accumPeriod = sum(BPstouse(:,180:end),2);
                else
                    accumPeriod = sum(BPstouse(:,1:180),2);
                    
                end
                
                [checkBP, cid] = sort(accumPeriod, 'descend');
                
                
                %% MOVE nan to bottom of order (not top)
                nanind= find(isnan(checkBP));
                %new start point
                cid=[cid(length(nanind)+1:end) ; cid(nanind)];
                %                     checkBP=[checkBP(length(nanind)+1:end); checkBP(nanind)]
                
                %rearrange.
                BPstouse=BPstouse(cid,:);
                datais=squeeze(datatouse(cid,usechan,:));
                
                %                 RTstouse=RTstouse(cid);
                %%
                
                
                imagesc(-3:1/60:3,1:length(cid),BPstouse)
                c=colorbar;
                ylabel(c, 'buttons pressed')
                %                 set(gca, 'ytick', 1:length(cid), 'yticklabel', round(sortedRTs./60,2))
                ylabel('PFI events')
                xlabel('Time [sec]')
                hold on
                plot([0 0] , ylim, ['w:'], 'linewidth', 3)
                title({['sorted BP data for ' ctype ','];['ppant' num2str(ifol)]})
                set(gca, 'fontsize', 15)
                
                %% show mean over time.
                subplot(4,2,7)
                plot([-3:1/60:3], nanmean(BPstouse,1),['k-'], 'linewidth', 3)
                %         shadedErrorBar([-3:1/60:3], nanmean(acrBP,1),nanstd(acrBP,1))
                set(gca, 'fontsize', 15)
                hold on;
                ylabel('nanmean BP ')
                %
                xlabel('Time secs')
                axis tight
                plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                
                
                set(gca, 'fontsize', 15)
                %% now plot spectrogram SNR.
                
                datais=squeeze(datatouse(cid,usechan,:));
                
                
                %
                %             %rmvbaseline from EEG.
                bsdata = zeros(size(datais));
                
                for itrial = 1:size(datais,1)
                    td = detrend(datais(itrial,:), 'linear');
                    %                     tdrm= mean(td(1,1:250));
                    %                     rmb= repmat(tdrm, [1 length(td)]);
                    bsdata(itrial,:) = td;%-rmb;
                end
                datais=bsdata;
                %%
                [sgrm ,tgrm, fgrm] = mtspecgramc(datais', movingwin, param_spcgrm);
                %%
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
                %%
                % smooth across trials
                sm_snrgrm20=zeros(size(snrgrm20));
                for itime=1:size(snrgrm20,1)
                    sm_snrgrm20(itime,:)= smooth(snrgrm20(itime,:),5);
                end
                
                subplot(4,2,[4 6]);
                imagesc(tgrm-3,  1:size(snr_sgrm,3), sm_snrgrm20');
                c=colorbar;
                ylabel(c, 'SNR')
                xlabel('Time secs')
                caxis([-1*max(max(sm_snrgrm20)) max(max(sm_snrgrm20))])
                caxis([0  max(max(sm_snrgrm20))])
                hold on
                plot([0 0] , ylim, ['w:'], 'linewidth', 3)
                title({['sorted ' num2str(usehz) 'Hz SNR data (smoothed) for '];[ ctype ', ppant' num2str(ifol)]})
                set(gca, 'fontsize', 15)
                
                
                %         ylim([0 20])
                %% show mean over time.
                subplot(4,2,8)
                plot(tbase, nanmean(snrgrm20,2),['k-'], 'linewidth', 3)
                %         shadedErrorBar([-3:1/60:3], nanmean(acrBP,1),nanstd(acrBP,1))
                set(gca, 'fontsize', 15)
                hold on;
                ylabel('nanmean SNR ')
                %%
                axis tight
                plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                %
                xlabel('Time secs')
                set(gcf,'color', 'w');
                
                %         ylim([0 20])
                
                %%
                
                print('-dpng', ['figure_PFI_' ctype '_summary_BPandSSVEP_' num2str(usehz) 'ppant ' num2str(ifol) '.png']);
                
                
                switch itimezero
                    case 1 %store for across ppant plots:
                        Ppant_onsetBP=BPstouse;
                        Ppant_onsetSNR=snrgrm20; %sorted.
                        %                         Ppant_onsetRTs=RTstouse;
                        Ppant_onsetTOPO=snr20;
                    case 2
                        
                        Ppant_offsetBP=BPstouse;
                        Ppant_offsetSNR=snrgrm20; %always sorted in descending order of PFI.
                        %                         Ppant_offsetRTs=RTstouse;
                        Ppant_offsetTOPO=snr20;
                end
            end
            
            
            %%
            savename=['PFIperformance_withSNR_' num2str(usehz)];
            
            cd(basefol)
            cd('EEG')
            cd(dirs(ifol).name)
            save(savename,...
                'Ppant_onsetBP','Ppant_offsetBP',...
                'Ppant_onsetSNR', 'Ppant_offsetSNR', ...
                'Ppant_onsetTOPO', 'Ppant_offsetTOPO', 'tgrm')
            
            
        end
    end
end





if job.ppantPFI_topotime==1
    getelocs
    window=[-3 3];
    
    srate=250;
    rmvbase=0;
    epochdur = sum(abs(window))*srate;
    
    timeid = [0:1/srate:epochdur];
    timeid= timeid-3;
    
    onsetc = ceil(epochdur)/2;
    % peakfreqsare=[20,40]; %hz
    %timing
    tt = 0:1/srate:60;
    snrmethod=2;
    %%
    param_spctrm.tapers = [1 1];
    param_spctrm.Fs= [250];
    param_spctrm.Fpass= [0 50];
    param_spctrm.trialave=0;
    
    %%
    for ifol = allppants
        
        for hzis=[1,2,7];
            usehz = peakfreqsare(hzis);
            
            
            
            icounter=1;
            cd(basefol)
            cd('EEG')
            cd(dirs(ifol).name)
            
            
            load(['ppant_PFI_Epoched'])
            
            for itimezero = 1:2
                if itimezero==1
                    
                    %append them all. TARG-> Disappearing (more buttons
                    %pressed).
                    datatouse = cat(1, ppant_SNREEG_PFI_0_1,ppant_SNREEG_PFI_1_2,ppant_SNREEG_PFI_2_3, ppant_SNREEG_PFI_3_4);
                    
                    ctype = 'PFI increase';
                    bsrem = [-3 -1]; %seconds
                    
                else
                    
                    datatouse = cat(1, ppant_SNREEG_PFI_1_0,ppant_SNREEG_PFI_2_1,ppant_SNREEG_PFI_3_2,ppant_SNREEG_PFI_4_3);
                    %                     RTstouse = [durs1_0'; durs2_1'; durs3_2'];
                    %                     BPstouse = cat(1, BPs1_0, BPs2_1, BPs3_2);
                    ctype = 'PFI decrease';
                    
                    
                    bsrem = [1 3]; %seconds
                    
                    
                end
                
                %             %rmvbaseline from EEG.
                datais=datatouse;
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
                %                 snr_grmout=zeros(64,14); % chans by time points.
                snr_grmout=zeros(64,size(datais,1),14);
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
                    
                    %                 snr_grmout(ichan,:)=squeeze(nanmean(snrgrm20,2));
                    snr_grmout(ichan,:,:)=snrgrm20';
                end
                %%
                
                switch itimezero
                    case 1 %store for across ppant plots:
                        
                        Ppant_onsetSNR_allchan=snr_grmout;
                        
                    case 2
                        
                        
                        Ppant_offsetSNR_allchan=snr_grmout;
                        
                end
            end
            
            
            %%
            
            savename=['PFIperformance_withSNR_' num2str(usehz)];
            
            
            
            save(savename,...
                'Ppant_onsetSNR_allchan','Ppant_offsetSNR_allchan','-append');
            
            disp([ 'fin all chan append for ' num2str(ifol)])
        end
    end
end


if job.concaterpdataacrossppants==1
    
    
    dursMINIMUM= 0; %1 second. % or change to  %this works.
    %     dursMINIMUM=123; % selects top third!.
    
    ppantsmoothing=1; % average across participants, after smoothing, or no.
    
    for ihz=[1:4,7]
        
        usehz=peakfreqsare(ihz);
        loadname=['PFIperformance_withSNR_' num2str(usehz)];
        
        storeacrossPpant_onsetBP=[];
        storeacrossPpant_onsetSNR=[];
        storeacrossPpant_onsetRTs=[];
        storeacrossPpant_onsetTOPO=[];
        storeacrossPpant_offsetBP=[];
        storeacrossPpant_offsetSNR=[];
        storeacrossPpant_offsetRTs=[];
        storeacrossPpant_offsetTOPO=[];
        
        icounter=1;
        
        for ippant = allppants
            cd(basefol)
            cd('EEG')
            cd(dirs(ippant).name);
            %onset types
            
            load(loadname)
            Ppant_onsetSNR=Ppant_onsetSNR';
            Ppant_offsetSNR=Ppant_offsetSNR';
            
            %restrict to only 'long' PFIs if needed.
            try
                if dursMINIMUM>0 && dursMINIMUM~=123
                    
                    %restrict to those trials of interest!
                    all_onset = 1:size(Ppant_onsetBP,1);
                    all_offset = 1:size(Ppant_offsetBP,1);
                    
                    shrton= find(Ppant_onsetRTs<dursMINIMUM); %note switch
                    all_onset(shrton)=[];
                    keepON=all_onset;
                    
                    %same for offsets.
                    shrtoff=find(Ppant_offsetRTs<dursMINIMUM);
                    all_offset(shrtoff)=[];
                    keepOFF=all_offset;
                    
                    Ppant_onsetBP=Ppant_onsetBP(keepON,:);
                    Ppant_offsetBP=Ppant_offsetBP(keepOFF,:);
                    
                elseif dursMINIMUM==123
                    
                    %select top third of disap durations.
                    %restrict to those trials of interest!
                    all_onset = 1:size(Ppant_onsetBP,1);
                    all_offset = 1:size(Ppant_offsetBP,1);
                    
                    
                    %find top third of distribution.
                    [ndurs,nid] = sort(Ppant_onsetRTs, 'descend');
                    %top third
                    keepON = nid(1:ceil(length(ndurs)/4));
                    
                    %same for offsets.
                    %find top third of distribution.
                    [ndurs,nid] = sort(Ppant_offsetRTs, 'descend');
                    %top third
                    keepOFF = nid(1:ceil(length(ndurs)/3));
                    
                    Ppant_onsetBP=Ppant_onsetBP(keepON,:);
                    Ppant_offsetBP=Ppant_offsetBP(keepOFF,:);
                    
                end
                
                
                
                
                if any(isnan(Ppant_onsetBP(:,1)))
                    lc= min(find(isnan(Ppant_onsetBP(:,1))));
                    lc=lc-1;
                else
                    lc=size(Ppant_onsetBP,1);
                end
                
                % we need to resample the BP and SNR data, to equate across
                % trial types    (ignoring nan)
                % currently at 100 fs.
                %%
                Ppant_onsetBP= resample(Ppant_onsetBP(1:lc,:),100,lc);
                Ppant_onsetSNR= resample(Ppant_onsetSNR(1:lc,:),100,lc);
                
                
                %offsets too
                if any(isnan(Ppant_offsetBP(:,1)));
                    lc= min(find(isnan(Ppant_offsetBP(:,1))));
                    lc=lc-1;
                else
                    lc=size(Ppant_offsetBP,1);
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
                        pdataN(1:16,:) = repmat(pData(1,:), [16,1]);
                        %and bottom;
                        pdataN(end-15:end,:) = repmat(pData(end,:), [16,1]);
                        
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
                
                
                
                
                storeacrossPpant_onsetBP(icounter,:,:)=Ppant_onsetBP;
                storeacrossPpant_onsetSNR(icounter,:,:)=Ppant_onsetSNR;
                %                 storeacrossPpant_onsetRTs(icounter,:)=Ppant_onsetRTs;
                storeacrossPpant_onsetTOPO(icounter,:)=Ppant_onsetTOPO';
                
                
                storeacrossPpant_offsetBP(icounter,:,:)=Ppant_offsetBP;
                storeacrossPpant_offsetSNR(icounter,:,:)=Ppant_offsetSNR;
                %                 storeacrossPpant_offsetRTs(icounter,:)=Ppant_offsetRTs;
                storeacrossPpant_offsetTOPO(icounter,:)=Ppant_offsetTOPO';
                
                
                icounter=icounter+1;
            catch SETUP
                if ippant~=27
                    rethrow(SETUP)
                end
            end
        end
        
        %save appropriately
        cd(basefol)
        cd('EEG')
        
        
        switch ihz
            case 1
                savename=['GFX_PFIperformance_withSNR_15_min' num2str(dursMINIMUM) ];
            case 2
                savename=['GFX_PFIperformance_withSNR_20_min' num2str(dursMINIMUM) ];
            case 3
                savename=['GFX_PFIperformance_withSNR_30_min' num2str(dursMINIMUM) ];
            case 4
                savename=['GFX_PFIperformance_withSNR_40_min' num2str(dursMINIMUM) ];
        end
        
        
        save(savename,...
            'storeacrossPpant_onsetBP','storeacrossPpant_offsetBP',...
            'storeacrossPpant_onsetSNR', 'storeacrossPpant_offsetSNR', ...
            'storeacrossPpant_onsetTOPO', 'storeacrossPpant_offsetTOPO', 'tgrm')
        
        
        
    end
    
    
    
end


if job.concaterpdataacrossppants_keepnumberSeparate==1
    % as the previous, except now we will store the outgoing MEAN
    %SSVEP SNR for accumperiods in thirds 0-1 1-2 2and greater
    
    
    dursMINIMUM= 0; %1 second. % or change to  %this works.
    %     dursMINIMUM=123; % selects top third!.
    
    ppantsmoothing=1; % average across participants, after smoothing, or no.
    
    for ihz=1:2
        
        if ihz==1
            loadname='PFIperformance_withSNR_20';
        else
            loadname='PFIperformance_withSNR_40';
        end
        storeacrossPpant_onsetSNR_PFI1=zeros(21,33);
        storeacrossPpant_onsetSNR_PFI2=zeros(21,33);
        storeacrossPpant_onsetSNR_PFI3=zeros(21,33);
        
        storeacrossPpant_offsetSNR_PFI1=zeros(21,33);
        storeacrossPpant_offsetSNR_PFI2=zeros(21,33);
        storeacrossPpant_offsetSNR_PFI3=zeros(21,33);
        
        icounter=1;
        
        for ippant = allppants
            cd(basefol)
            cd(num2str(ippant));
            %onset types
            
            load(loadname)
            Ppant_onsetSNR=Ppant_onsetSNR';
            Ppant_offsetSNR=Ppant_offsetSNR';
            
            %restrict to only 'long' PFIs if needed.
            
            if dursMINIMUM>0 && dursMINIMUM~=123
                
                %restrict to those trials of interest!
                all_onset = 1:size(Ppant_onsetBP,1);
                all_offset = 1:size(Ppant_offsetBP,1);
                
                shrton= find(Ppant_onsetRTs<dursMINIMUM); %note switch
                all_onset(shrton)=[];
                keepON=all_onset;
                
                %same for offsets.
                shrtoff=find(Ppant_offsetRTs<dursMINIMUM);
                all_offset(shrtoff)=[];
                keepOFF=all_offset;
                
                Ppant_onsetBP=Ppant_onsetBP(keepON,:);
                Ppant_offsetBP=Ppant_offsetBP(keepOFF,:);
                
            elseif dursMINIMUM==123
                
                %select top third of disap durations.
                %restrict to those trials of interest!
                all_onset = 1:size(Ppant_onsetBP,1);
                all_offset = 1:size(Ppant_offsetBP,1);
                
                
                %find top third of distribution.
                [ndurs,nid] = sort(Ppant_onsetRTs, 'descend');
                %top third
                keepON = nid(1:ceil(length(ndurs)/4));
                
                %same for offsets.
                %find top third of distribution.
                [ndurs,nid] = sort(Ppant_offsetRTs, 'descend');
                %top third
                keepOFF = nid(1:ceil(length(ndurs)/3));
                
                Ppant_onsetBP=Ppant_onsetBP(keepON,:);
                Ppant_offsetBP=Ppant_offsetBP(keepOFF,:);
                
            end
            
            
            
            
            if any(isnan(Ppant_onsetBP(:,1)))
                lc= min(find(isnan(Ppant_onsetBP(:,1))));
                lc=lc-1;
            else
                lc=size(Ppant_onsetBP,1);
            end
            
            % we need to resample the BP and SNR data, to equate across
            % trial types    (ignoring nan)
            % currently at 100 fs.
            %%
            Ppant_onsetBP= resample(Ppant_onsetBP(1:lc,:),100,lc);
            Ppant_onsetSNR= resample(Ppant_onsetSNR(1:lc,:),100,lc);
            
            
            %offsets too
            if any(isnan(Ppant_offsetBP(:,1)));
                lc= min(find(isnan(Ppant_offsetBP(:,1))));
                lc=lc-1;
            else
                lc=size(Ppant_offsetBP,1);
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
                    pdataN(1:16,:) = repmat(pData(1,:), [16,1]);
                    %and bottom;
                    pdataN(end-15:end,:) = repmat(pData(end,:), [16,1]);
                    
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
            
            % now sort into chunks - median split?
            
            for id=1:2 %onsets and offsets
                switch id
                    case 1
                        dtype=Ppant_onsetBP;
                        %                                 accumPeriod = (sum(dtype(:, 180:end),2)/60)/3;
                        
                    case 2
                        dtype=Ppant_offsetBP;
                        %                                 accumPeriod = (sum(dtype(:, 1:180),2)/60)/3;
                end
                
                %now we have accum in seconds, sort into 3 separate
                %traces.
                %                         p1t = find(accumPeriod>0 & accumPeriod<1);
                %                         p2t = find(accumPeriod>1 & accumPeriod<2);
                %                         p3t = find(accumPeriod>2);
                
                %or find max npressed in a trial:
                % maxN = max(dtype,[],2);
                %   p1t = find(maxN<1);
                %   p2t = find(maxN>1 & maxN<2);
                %   p3t = find(maxN>2 & maxN<3);
                %                         p2t = find(accumPeriod>1 & accumPeriod<2);
                %                         p3t = find(accumPeriod>2);
                
                
                switch id
                    case 1
                        p1D_on= mean(Ppant_onsetSNR(51:100,:),1);
                        p2D_on= mean(Ppant_onsetSNR(1:50,:),1);
                        %                                 p3D_on= mean(Ppant_onsetSNR(1:33,:),1);
                        %                                 if isempty(p3t)
                        %                                 p3D_on= zeros(1,33); %anovan can't handle nan.
                        %                                 else
                        %                                     p3D_on= mean(Ppant_onsetSNR(p2t,:),1);
                        %                                 end
                    case 2
                        p1D_off= mean(Ppant_offsetSNR(51:100,:),1);
                        p2D_off= mean(Ppant_offsetSNR(1:50,:),1);
                        %                                 p3D_off= mean(Ppant_offsetSNR(1:33,:),1);
                        %                                 if isempty(p3t)
                        %                                 p3D_off= zeros(1,33); %anovan can't handle nan.
                        %                                 else
                        %                                     p3D_off= mean(Ppant_offsetSNR(p3t,:),1);
                        %                                 end
                        
                        
                end
                
            end
            
            %% sanity check:
            %                     figure(1); clf
            %                     plot(1:size(p1D_off,2), p1D_off), hold on;
            %                     plot(1:size(p1D_off,2), p2D_off, 'k');
            %                     title(num2str(ippant))
            %                     shg
            %                     figure(2);imagesc(Ppant_offsetSNR); colorbar, caxis([-5 5]); hold on, plot([xlim], [50 50], ['-k'])
            %                     %%
            storeacrossPpant_onsetSNR_PFI1(icounter,:)= p1D_on;
            storeacrossPpant_onsetSNR_PFI2(icounter,:)= p2D_on;
            %                     storeacrossPpant_onsetSNR_PFI3(icounter,:)= p3D_on;
            
            storeacrossPpant_offsetSNR_PFI1(icounter,:)= p1D_off;
            storeacrossPpant_offsetSNR_PFI2(icounter,:)= p2D_off;
            %                     storeacrossPpant_offsetSNR_PFI3(icounter,:)= p3D_off;
            icounter=icounter+1;
        end
        %save appropriately
        cd(basefol)
        cd('newplots-MD')
        
        
        switch ihz
            case 1
                savename=['GFX_PFIperformance_withSNR_20_min' num2str(dursMINIMUM) ];
            case 2
                savename=['GFX_PFIperformance_withSNR_40_min' num2str(dursMINIMUM) ];
        end
        
        
        save(savename,...
            'storeacrossPpant_onsetSNR_PFI1',...
            'storeacrossPpant_onsetSNR_PFI2',...
            'storeacrossPpant_onsetSNR_PFI3',...
            'storeacrossPpant_offsetSNR_PFI1',...
            'storeacrossPpant_offsetSNR_PFI2',...
            'storeacrossPpant_offsetSNR_PFI3','-append')
        
        
        
    end
    
    
    
end

%%
if job.concatTOPOTIMEacrossppants==1
    %%
    for verData=2
        
        for ihz=[1,2,7]
            
            usehz=peakfreqsare(ihz);
            loadname=['PFIperformance_withSNR_' num2str(usehz) ];
            
            storeacrossPpant_onsetSNR_chans=zeros(length(allppants),64,14);
            storeacrossPpant_offsetSNR_chans=zeros(length(allppants),64,14);
            
            icounter=1;
            if verData==1
                for ippant = allppants
                    cd(basefol)
                    cd('EEG')
                    cd(pdirs(ippant).name);
                    %onset types
                    
                    load(loadname)
                    %             storeacrossPpant_onsetSNR_chans(icounter,:,:)=squeeze(nanmean(Ppant_onsetSNR_allchan,2));
                    %             storeacrossPpant_offsetSNR_chans(icounter,:,:)=squeeze(nanmean(Ppant_offsetSNR_allchan,2));
                    storeacrossPpant_onsetSNR_chans(icounter,:,:)=Ppant_onsetSNR_allchan;
                    storeacrossPpant_offsetSNR_chans(icounter,:,:)=Ppant_offsetSNR_allchan;
                    
                    icounter=icounter+1;
                end
                %%
                cd(basefol)
                cd('EEG')
                cd('GFX_Pre-RESS')
                %%
                
                savename=['GFX_PFIperformance_withSNR_' num2str(usehz) '_allchan'];
                
                save(savename, 'storeacrossPpant_onsetSNR_chans',...
                    'storeacrossPpant_offsetSNR_chans', 'tgrm');
                
                
            else % downsample PFI.
                %downsample to how many trials?
                
                dsnum=48;
                %how many times?
                nshuff=100;
                
                
                storeacrossPpant_onsetSNR_chans = zeros(length(allppants),64,nshuff,14); %retain across shuffles.
                storeacrossPpant_offsetSNR_chans=storeacrossPpant_onsetSNR_chans ;
                
                for ippant = allppants
                    
                    
                    cd(basefol)
                    cd('EEG')
                    cd(pdirs(ippant).name);
                    
                    %onset types
                    
                    load(loadname)
                    
                    % downsample with replacement,
                    dsdata_ons = zeros(64, nshuff,14); %
                    dsdata_offs= zeros(64, nshuff,14); %
                    
                    
                    for ids=1:nshuff
                        %random trials (since 24 catch trials (max).
                        dstrials= randi([1 size(Ppant_onsetSNR_allchan,2)],[1 dsnum]);
                        % % dstrials=randperm(size(Ppant_onsetSNR_allchan,2), dsnum);
                        %
                        dsdata_ons(:,ids,:) =squeeze(nanmean(Ppant_onsetSNR_allchan(:, dstrials,:),2));
                        
                        %can have different trial numbers for onsets and
                        %offsets. so reperform.
                        dstrials= randi([1 size(Ppant_offsetSNR_allchan,2)],[1 dsnum]);
                        %                         dstrials=randperm(size(Ppant_offsetSNR_allchan,2), dsnum);
                        dsdata_offs(:,ids,:) =squeeze(nanmean(Ppant_offsetSNR_allchan(:, dstrials,:),2));
                    end
                    
                    % mean over permutations.
                    storeacrossPpant_onsetSNR_chans(icounter,:,:,:)=dsdata_ons;
                    storeacrossPpant_offsetSNR_chans(icounter,:,:,:)=dsdata_offs;
                    
                    icounter=icounter+1;
                    
                    
                    
                end
                %%
                cd(basefol)
                cd('EEG')
                cd('GFX_Pre-RESS')
                %%
                %                 cd('newplots-MD')
                
                
                savename=['GFX_PFIperformance_withSNR_' num2str(usehz) '_allchan'];
                
                
                %adjust variable name
                storeacrossPpant_onsetSNR_chans_downsampled =storeacrossPpant_onsetSNR_chans;
                storeacrossPpant_offsetSNR_chans_downsampled =storeacrossPpant_offsetSNR_chans;
                
                save(savename, 'storeacrossPpant_onsetSNR_chans_downsampled',...
                    'storeacrossPpant_offsetSNR_chans_downsampled', 'tgrm', '-append');
                
            end
        end
        
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
    
    
    
    for hzis=3:4
        cd(basefol)
        cd('EEG')
        switch hzis
            case 1
                
                load('GFX_PFIperformance_withSNR_15_min0')
                
                usehz=15;
            case 2
                
                load('GFX_PFIperformance_withSNR_20_min0')
                
                usehz=20;
            case 3
                
                load('GFX_PFIperformance_withSNR_30_min0')
                
                usehz=30;
            case 4
                
                load('GFX_PFIperformance_withSNR_40_min0')
                
                usehz=40;
        end
        if hzis==3
            clf
        end
        tbase = tgrm-3;
        %%
        for itimezero=1
            
            
            switch itimezero
                case 1
                    useBP=storeacrossPpant_onsetBP;
                    useSNR=storeacrossPpant_onsetSNR;
                    %                     useRTs=storeacrossPpant_onsetRTs;
                    useTOPO=storeacrossPpant_onsetTOPO;
                    
                    ctype= 'PFI increase';
                    chtype='button press';
                    
                case 2
                    
                    useBP=storeacrossPpant_offsetBP;
                    useSNR=storeacrossPpant_offsetSNR;
                    %                     useRTs=storeacrossPpant_offsetRTs;
                    useTOPO=storeacrossPpant_offsetTOPO;
                    
                    ctype= 'PFI decrease';
                    chtype='button release';
                    
                    
            end
            
            
            
            %                 for ip= 1:length(allppants)
            %                     trialind = [1:24] + 24*(ip-1);
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
            %             clf
            %% not topoplotting anymore.
            %             usechan=62;
            %             subplot(2,2,1:2)
            %             topoplot(acrTOPO, elocs(1:64), 'emarker2', {usechan, 'o', 'w'});
            %             c=colorbar;
            %             caxis([0 15])
            %             ylabel(c, {['10log10(SNR)'];['-3000:-100ms']})
            %             hold on
            %             title({[num2str(usehz) ' Hz SSVEP']})
            %             set(gca, 'fontsize', 20)
            %%
            subplot(311)
            %         [sortedRTs, cid] = sort(acrRTs, 'descend');
            
            %
            imagesc([-3:1/60:3], 1:size(acrBP,1), acrBP);%(cid,:));
            %
            title({['Buttons Pressed'] }, 'fontsize', 25)
            c=colorbar;
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
            xlabel(['Time from ' ctype])
            
            %SNR
            try
                subplot(3,1, hzis+1);
            catch
                subplot(3,1, hzis-1);
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
            %             ylabel(c, {['SNR(dB/Hz)'];['baseline corrected']})
            ylabel(c, {['SNR(dB/Hz)']})
            %
            %             caxis([-1*max(max(sm_snrgrm20)) max(max(sm_snrgrm20))])
            %             caxis([0 1])
            %                     caxis([10 15])
            title([num2str(usehz) 'Hz dyn-SSVEP SNR'])
            %             set(gca, 'ytick', 1:24, 'yticklabel', round(sortedRTs./60,2))
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
            x = ppantMeanSNR;
            
            mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
            mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
            
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
            xlabel(['Time from ' ctype])
            ylabel({'Normalized'; 'trial count'})
            
            set(gca, 'fontsize', 25)
            cd(basefol)
            cd('Figures')
            set(gcf, 'color', 'w')
            shg
            %%
            colormap('jet')
            shg
        end
    end
    %%
    shg
    %%
    cd([basefol filesep 'Figures' filesep 'SNR trial-by-trial'])
    %%
    print('-dpng', ['PFI SNR summary ' ctype ' 2ndharm.png'])
end

%%
if job.BPandSSVEPtimecourseacrossppants==1
    %%
    getelocs
    
    rmvbase=0;
    checksigON=0;
    checkcluster=0;
    
    clf
    plcount=1;
    legendprint=[];
    legendnames=[];
    ttestdata=[];
    lc=1;
    
    
    for hzis=1:2%
        cd([basefol filesep 'EEG'])
        switch hzis
            case 1
                load('GFX_PFIperformance_withSNR_15_min0')
                
                lint=':';
                usehz=15;
                sigheight=2;
                if rmvbase==1
                    %                     sigheight=-3.7;
                    sigheight=1.9;
                end
                col='r';
                
            case 2
                load('GFX_PFIperformance_withSNR_20_min0')
                
                lint=':';
                usehz=20;
                sigheight=2;
                if rmvbase==1
                    %                     sigheight=-3.7;
                    sigheight=1.9;
                end
                col=[.5 .5 .5];
            case 3
                load('GFX_PFIperformance_withSNR_30_min0')
                
                lint=':';
                usehz=30;
                sigheight=2;
                if rmvbase==1
                    %                     sigheight=-3.7;
                    sigheight=1.9;
                end
                col='r';
            case 4
                
                load('GFX_PFIperformance_withSNR_40_min0')
                
                usehz=40;
                lint=':';
                sigheight=1.85;
                if rmvbase==1
                    %                     sigheight=-4.1;
                    sigheight=1.5;
                end
                
                col=[.5 .5 .5];
        end
        tbase = tgrm-3;
        %%
        
        
        %         legendprint=[];
        for itimezero=1:2
            
            
            switch itimezero
                case 1
                    useBP=storeacrossPpant_onsetBP;
                    useSNR=storeacrossPpant_onsetSNR;
                    
                    useTOPO=storeacrossPpant_onsetTOPO;
                    
                    chtype='target invisible';
                    %                     chtype='button press';
                    
                    %                     col=[0 .5 0];
                    
                case 2
                    
                    useBP=storeacrossPpant_offsetBP;
                    useSNR=storeacrossPpant_offsetSNR;
                    
                    useTOPO=storeacrossPpant_offsetTOPO;
                    
                    
                    chtype='target visible';
                    %                     chtype='button release';
                    
                    
                    %                     col='r';
                    bsrem=[-3 -1];
                    
                    
            end
            
            
            %%
            %
            %SNR
            figure(1);
            %
            
            for ip=1:2
                switch ip
                    case 1
                        used=useBP;
                        basep=1;
                        col='b';
                        lint='-';
                        ylimsare = [0 4];
                    case 2
                        used=useSNR;
                        basep=3;
                        if hzis==1 ||hzis==3
                            col='r';
                            lgc=1;
                        elseif hzis==2|| hzis==4
                            col=[0, .5, 0];
                            lgc=2;
                        end
                        lint=':';
                        %                         ylimsare= [ 1.2 1.9];
                end
                
                placeme = basep+ (1*itimezero-1);
                
                subplot(2,2,placeme)
                
                
                hold on
                
                %
                if rmvbase==1
                    useSNR=useSNR-mean(useSNR(:));
                end
                %plot across ppant trace
                ppantMeanSNR= squeeze(mean(used,2));
                
                %adjust standard error as per COusineau(2005)
                %confidence interval for within subj designs.
                % y = x - mXsub + mXGroup,
                if ip~=1
                    x = ppantMeanSNR;
                    
                    mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
                    mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
                    
                    %for each observation, subjtract the subj average, add
                    %the group average.
                    NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
                    
                    %             %compute new stErr %which version?
                    
                    
                else
                    NEWdata=ppantMeanSNR;
                end
                stE = std(NEWdata)/sqrt(size(NEWdata,1));
                
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
                xlabel(['Time from ' chtype ' [s]'])
                %             xlabel('Time from perceptual report')
                set(gca, 'fontsize', 25)
                xlim([-3 3])
                ylim( [ylimsare])
                
                
                
                
                
                set(gcf, 'color', 'w')
                
                if ip==2 %then SNR data
                    if itimezero==1 && hzis==1
                        ttestdata(1,:,:) = ppantMeanSNR;
                    elseif itimezero==1 && hzis==2
                        ttestdata(2,:,:) = ppantMeanSNR;
                        
                    elseif itimezero==2 && hzis==1
                        ttestdata(3,:,:) = ppantMeanSNR;
                    elseif itimezero==2 && hzis==2
                        ttestdata(4,:,:) = ppantMeanSNR;
                        
                    end
                    legendprint(lgc)=sh.mainLine;
                end
                
                plcount=plcount+1;
                
                
                %place legend
                if hzis==2 && ip==2
                    if itimezero==1
                        lg=legend([legendprint(1) legendprint(2)], {'20 Hz', '40 Hz'});
                    end
                    set(lg, 'location', 'NorthEast')
                    
                end
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
                        lg=legend([legendprint(1) legendprint(3)], {legendnames{1}, legendnames{3}});
                        set(lg, 'location', 'NorthEast')
                    end
                    
                    
                    
                end
            end
            
            
            
            
            if ip==2
                axis tight
            end
            %%
            if checksigON==1 && hzis==2
                %check for sig
                
                if itimezero==1 %accomodate for the indexing in ttest data.
                    testSNRs = [1,2];
                else
                    testSNRs = [3,4];
                end
                pvals=zeros(1,size(ttestdata,3));
                tvals=zeros(1,size(ttestdata,3));
                for itime = 1:size(ttestdata,3)
                    
                    try [h,pvals(itime),~,stat]=ttest(ttestdata(testSNRs(1),:,itime), ttestdata(testSNRs(2),:,itime));
                        shuffType=1;
                    catch
                        %                 [h,pvals(itime),~,stat]=ttest(ttestdata(plcount-1,:,itime)); %compares to zero.
                        %                 shuffType=2; %whether or not to skip the non-parametric test for sig.
                    end
                    
                    tvals(itime)= stat.tstat;
                end
                sigs=find(pvals<.05);
                %
                %             %perform cluster based correction.
                if length(sigs)>2 &&checkcluster==1
                    % find biggest cluster:
                    %finds adjacent time points
                    vect1 = diff(sigs);
                    v1 = (vect1(:)==1);
                    d = diff(v1);
                    clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
                    %grab largest
                    %                 ignore bad points.
                    
                    
                    % find biggest cluster:
                    %finds adjacent time points
                    sigs = find(pvals<.05);
                    
                    vect1 = diff(sigs);
                    v1 = (vect1(:)==1);
                    d = diff(v1);
                    clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
                    [~,maxClust] = max(clusterSTandEND(:,2)-clusterSTandEND(:,1));
                    
                    %%
                    for icl=1%:size(clusterSTandEND,1)
                        
                        %start and end are now:
                        % change icl to maxClust if we only want the largest
                        % cluster.
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
                                    for ippant = 1:size(ttestdata,2)
                                        for itime=1:length(checktimes)
                                            
                                            if mod(randi(100),2)==0 %if random even number
                                                pdata = ttestdata(1,randi(size(ttestdata,2)), checktimes(itime)); %select both chans
                                            else
                                                pdata = ttestdata(2,randi(size(ttestdata,2)), checktimes(itime));
                                            end
                                            
                                            shD(ipartition,itime,ippant) = pdata;
                                        end
                                    end
                                end
                            else %null is that there are no temporal coincident sig values.
                                for ipartition = 1:2
                                    for ippant = 1:size(ttestdata,2)
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
                        if itimezero==1
                            clf
                        end
                        subplot(2,1, itimezero)
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
                                figure(1); subplot(2,2,placeme)
                                hold on
                                
                                %                             plot(tgrm(itime)-3, sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', sh.mainLine.Color)
                                plot(tgrm(itime)-3, sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'k')
                            end
                        end
                        
                        %
                    end
                end
            end
            
        end
        %%
    end
    
    hold on
    %             plot(xlim, [0 0], ['k:']);
    %             plot([0 0] , ylim, ['k:'], 'linewidth', 1)
    %plot patch:
    
    
    shg
    
    hold on
    
    %%
    
    
    %%
    cd(basefol)
    cd('Figures')
    cd('SNR-timecourse')
    %%
    print('-dpng', ['PFI trace Bground SSVEP summary, during PFI 2nd harm- meansubtracted.png'])
    %         print('-dpng', ['PFI trace Bground SSVEP summary, during  PFI 2nd harm.png'])
    
    %%
    
    % print('-dpng', ['PFI trace Bground SSVEP summary, during both, 40Hz.png'])
end

if job.SSVEPtimecourseSelectchans==1
    
    %%
    getelocs
    
    rmvbase=0;
    checkcluster=0;
    uniqueto40hz=1;
    clf
    lateonly=1;
    plcount=1;
    for hzis=1:2
        cd(basefol)
        cd('newplots-MD')
        switch hzis
            case 1
                load('GFX_PFIperformance_withSNR_20_allchan')
                
                lint='-';
                usehz=20;
                sigheight=4.1;
                if rmvbase==1
                    sigheight=-3.7;
                end
                
                
                
                if lateonly~=1
                    checkchans1= [25,62,30,61,29,63,31]; % PFI 20Hz, early.
                    
                    
                else
                    checkchans1=[];
                end
                
                %combine above.
                checkchans2= [29,30,61:63]; % PFI 20Hz, late.
                usechans20 = unique([checkchans1, checkchans2]);
                
                usechans=usechans20;
                
                
                col='r';
            case 2
                
                load('GFX_PFIperformance_withSNR_40_allchan')
                
                usehz=40;
                %                 lint=':';
                sigheight=3.7;
                if rmvbase==1
                    sigheight=-4.1;
                end
                
                if lateonly~=1
                    checkchans1= [4 9 18 23 24 25 26 29 30 47 54 56 57 58 60 62 63];
                else
                    checkchans1=[];
                end
                checkchans2=[10,44,53,58,59,26,60:63,29:31];
                
                usechans40= unique([checkchans1,checkchans2]) ;
                
                
                
                if uniqueto40hz ==1
                    remv=find(ismember(usechans40,usechans20));
                    usechans40(remv)=[];
                end
                usechans=usechans40;
                
                col=[0 .5 0];
        end
        tbase = tgrm-3;
        %%
        %         usechans=30;
        ttestdata=[];
        for itimezero=2
            
            
            switch itimezero
                case 1
                    
                    useSNR=squeeze(nanmean(storeacrossPpant_onsetSNR_chans(:,usechans,:),2));
                    
                    
                    
                    chtype='target invisible';
                    %                     chtype='button press';
                    
                    %                     col=[0 .5 0];
                    
                case 2
                    
                    
                    useSNR=squeeze(nanmean(storeacrossPpant_offsetSNR_chans(:,usechans,:),2));
                    
                    chtype='target visible';
                    %                     chtype='button release';
                    
                    
                    %                     col='r';
                    bsrem=[-3 -1];
                    
                    
            end
            
            
            %%
            %
            %SNR
            figure(1);
            plis = 3;%+(itimezero-1);
            %             subplot(2,2,plis)
            hold on
            
            
            if rmvbase==1
                tmp = zeros(size(useSNR));
                tIND = dsearchn(tbase', [bsrem]');
                for ippant=1:size(useSNR,1)
                    for itrial = 1:size(useSNR,2)
                        tbase = squeeze(nanmean(useSNR(ippant,itrial,tIND(1):tIND(2)),3));
                        rbase = repmat(tbase, [1  size(useSNR,3)]);
                        
                        tmp(ippant,itrial,:) =squeeze(useSNR(ippant,itrial,:))' - rbase;
                    end
                end
                useSNR=tmp;
            end
            
            %plot across ppant trace
            ppantMeanSNR= useSNR;
            
            %adjust standard error as per COusineau(2005)
            %confidence interval for within subj designs.
            % y = x - mXsub + mXGroup,
            x = ppantMeanSNR;
            
            mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
            mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
            
            %             %compute new stErr %which version?
            stE = std(NEWdata)/sqrt(size(x,1));
            %
            sh=shadedErrorBar(tgrm-3, mean(ppantMeanSNR,1),stE,[lint],[1]);
            sh.mainLine.LineWidth=3;
            
            sh.mainLine.Color = col;
            sh.patch.FaceColor = col;
            sh.edge(1).Color = col;
            sh.edge(2).Color = col;
            
            
            set(gca, 'fontsize', 15)
            hold on;
            %                 %%
            %
            ylabel({['log(SNR)']})
            
            %             title({[num2str(usehz) ' Hz SSVEP']})
            xlabel(['Time from ' chtype])
            %             xlabel('Time from perceptual report')
            set(gca, 'fontsize', 25)
            cd(basefol)
            cd('newplots-MD')
            
            
            cd('Figures')
            set(gcf, 'color', 'w')
            ttestdata(plcount,:,:) = ppantMeanSNR;
            legendprint(plcount)=sh.mainLine;
            plcount=plcount+1;
        end
        %%
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
        sigs=find(pvals<.05);
        %
        %             %perform cluster based correction.
        if length(sigs)>2 &&checkcluster==1
            % find biggest cluster:
            %finds adjacent time points
            vect1 = diff(sigs);
            v1 = (vect1(:)==1);
            d = diff(v1);
            
            %grab largest
            %                 ignore bad points.
            
            
            % find biggest cluster:
            %finds adjacent time points
            sigs = find(pvals<.05);
            
            vect1 = diff(sigs);
            v1 = (vect1(:)==1);
            d = diff(v1);
            clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
            [~,maxClust] = max(clusterSTandEND(:,2)-clusterSTandEND(:,1));
            
            %%
            for icl=1%:size(clusterSTandEND,1)
                
                %start and end are now:
                % change icl to maxClust if we only want the largest
                % cluster.
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
                            for ippant = 1:size(ttestdata,2)
                                for itime=1:length(checktimes)
                                    
                                    if mod(randi(100),2)==0 %if random even number
                                        pdata = ttestdata(1,randi(size(ttestdata,2)), checktimes(itime)); %select both chans
                                    else
                                        pdata = ttestdata(2,randi(size(ttestdata,2)), checktimes(itime));
                                    end
                                    
                                    shD(ipartition,itime,ippant) = pdata;
                                end
                            end
                        end
                    else %null is that there are no temporal coincident sig values.
                        for ipartition = 1:2
                            for ippant = 1:size(ttestdata,2)
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
                if plcount==3
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
                        %                             plot(tgrm(itime)-3, sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'm')
                    end
                end
                
                %
            end
        end
        
        
        %%
    end
    axis tight
    if rmvbase~=1
        ylim([-2 4.5])
    else
        ylim([-4.5 2])
    end
    hold on
    plot(xlim, [0 0], ['k:']);
    plot([0 0] , ylim, ['k:'], 'linewidth', 1)
    %plot patch:
    
    
    shg
    
    hold on
    
    %%
    
    %           lg=legend(legendprint, {'20 Hz', '40 Hz'});
    lg=legend(legendprint, {'PFI increase', 'PFI decrease'});
    if itimezero==1
        set(lg, 'location', 'NorthEast')
    end
    %%
    %         print('-dpng', ['PFI trace Bground SSVEP summary, during ' num2str(chtype) '.png'])
    
    %%
    
    
    
    
end



if job.FirstSIGintimecourse==1 %this uses
    
    getelocs
    
    rmvbase=1;
    checkcluster=1;
    clf
    plcount=1;
    ppantFIRSTSIG=nan(2,2,21);
    ttestdata=[];
    uniqueto40hz=0;
    
    icounter=1;
    
    for ifol=allppants
        cd(basefol)
        
        cd(num2str(ifol))
        
        for hzis=1:2
            
            %%
            
            
            switch hzis
                case 1
                    load('PFIperformance_withSNR_20.mat', 'tgrm', 'Ppant_onsetSNR_allchan','Ppant_offsetSNR_allchan');
                    checkchans1= [25,62,30,61,29,63,31]; % PFI 20Hz, early.
                    checkchans2= [29,30,61:63]; % PFI 20Hz, late.
                    checkchans1= []; % PFI 20Hz, late.
                    %combine above
                    usechans20 = unique([checkchans1, checkchans2]);
                    
                    usechans=usechans20;
                    
                    %       usechans=[57, 60:64,29:31]; %PFI increase
                    usechans=[25 29 30 58 60 61 62 63];%PFI decrease (20Hz).
                    
                case 2
                    load('PFIperformance_withSNR_40.mat', 'tgrm', 'Ppant_onsetSNR_allchan','Ppant_offsetSNR_allchan');
                    checkchans1= [4 9 18 23 24 25 26 29 30 47 54 56 57 58 60 62 63];
                    checkchans2=[10,44,53,58,59,26,60:63,29:31];
                    checkchans1=[];
                    usechans40= unique([checkchans1,checkchans2]) ;
                    
                    if uniqueto40hz ==1
                        remv=find(ismember(usechans40,usechans20));
                        usechans40(remv)=[];
                    end
                    usechans=usechans40;
                    
                    % or use ONSET chans: from fig 12
                    %                 usechans=[19 21 23 25 26 27 29 30 31 53 56 57 58 59 60 61 62 63 64];% PFI increase
                    
                    usechans=[15 19 23 24 25 26 29 30 31 50 52 53 54 56 57 58 59 60 61 62 63];%PFI decrease
                    
            end
            
            
            for itimezero=1:2
                
                
                switch itimezero
                    case 1
                        
                        useSNR=squeeze(nanmean(Ppant_onsetSNR_allchan(usechans,:,:),1));
                        
                        chtype='PFI increase';
                        %                     chtype='button press';
                        rmvbase=1;
                        col=[0 .5 0];
                        bsrem=[-3 -2];
                        
                    case 2
                        
                        
                        useSNR=squeeze(nanmean(Ppant_offsetSNR_allchan(usechans,:,:),1));
                        chtype='PFI decrease';
                        %                     chtype='button release';
                        rmvbase=1;
                        
                        col='r';
                        bsrem=[-3 -2];
                        
                        
                end
                
                tbase=tgrm-3;
                %%
                %
                %SNR
                
                hold on
                
                
                if rmvbase==1
                    tmp = zeros(size(useSNR));
                    tIND = dsearchn(tbase', [bsrem]');
                    
                    for itrial = 1:size(useSNR,1)
                        presub=squeeze(useSNR(itrial,:));
                        ttmp = squeeze(nanmean(useSNR(itrial,tIND(1):tIND(2)),2));
                        rbase = repmat(ttmp, [1  size(useSNR,2)]);
                        bremt=presub - rbase;
                        tmp(itrial,:) =bremt;
                    end
                    
                    useSNR=tmp;
                end
                
                %Calc first time point per ppant above zero.
                
                
                
                
                pvals=nan(1,size(useSNR,2));
                h=0;
                
                for itime = 1:size(useSNR,2)
                    
                    
                    [h,pvals(itime),~,stat]=ttest(useSNR(:,itime)); %compares to zero.
                    
                    
                    
                    tvals(itime)= stat.tstat;
                    
                    
                    switch itimezero
                        case 1
                            if squeeze(mean(useSNR(:,itime),1))<0
                                pvals(itime)=nan;
                            end
                        case 2
                            if squeeze(mean(useSNR(:,itime),1))>0
                                pvals(itime)=nan;
                            end
                    end
                end
                
                q=fdr(pvals,.05);
                
                %                 if q~=0
                %                 sigs=find(pvals<q);
                %                 else
                sigs=find(pvals<.05);
                %                 end
                %             %perform cluster based correction.
                
                %                 find biggest cluster:
                %finds adjacent time points
                vect1 = diff(sigs);
                v1 = (vect1(:)==1);
                if sum(v1)>0
                    d = diff(v1);
                    clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
                    %grab largest
                    %                 ignore bad points.
                    %
                    [~,maxClust] = max(clusterSTandEND(:,2)-clusterSTandEND(:,1));
                    % [mr,maxClust] = find((clusterSTandEND(:,2)-clusterSTandEND(:,1))>1);
                    %                maxClust=min(mr);
                    %start and end are now:
                    STC=sigs(clusterSTandEND(maxClust,1));
                    
                    if STC<1
                        errror('check')
                    end
                    
                    ppantFIRSTSIG(hzis,itimezero,icounter)= STC;
                    
                end
                %  ppantFIRSTSIG(hzis,itimezero,ippant)= min(sigs);
                
            end
            %     disp(['chans used hz ' num2str(hzis) ',= ' num2str(usechans)])
        end
        icounter=icounter+1;
    end
    
    %%
    figure(1);
    clf
    plcounter=1;
    for itime=1:2
        figure(1); hold on
        for ihz=1:2
            subplot(2,2,itime)
            switch ihz
                case 1
                    if itime==1
                        colb='g';
                    else
                        colb='r';
                    end
                case 2
                    colb='k';
            end
            h=histogram(ppantFIRSTSIG(ihz,itime,:));
            h.FaceColor= colb;
            
            hold on;
            pl(plcounter)=plot([nanmean(ppantFIRSTSIG(ihz,itime,:),3) nanmean(ppantFIRSTSIG(ihz,itime,:),3)], [ylim], 'linew', 4, 'color', colb);
            plcounter=plcounter+1;
            
            
            
        end
        [h,p,~,stat]= ttest(squeeze(ppantFIRSTSIG(1,itime,:)), squeeze(ppantFIRSTSIG(2,itime,:)));
        title(['t=' num2str(stat.tstat), 'p=' num2str(p)]);
        
    end
    %%
    figure(2); clf ;hold on
    
    
    
    itimezero=2;
    
    
    
    plot(tbase, ones(1, length(tbase)), 'color', 'w'); hold on; ylim([0 1])
    m1=tbase(floor(nanmean(ppantFIRSTSIG(1,itimezero,:),3)));
    m2=tbase(floor(nanmean(ppantFIRSTSIG(2,itimezero,:),3)));
    
    %adjust error bars.
    m1x=squeeze(ppantFIRSTSIG(1,itimezero,:));
    m2x=squeeze(ppantFIRSTSIG(2,itimezero,:));
    
    
    mXppant =squeeze( nanmean(ppantFIRSTSIG(:,itimezero,:),1)); %mean across conditions we are comparing (within ppant ie. hztypes).
    mXgroup = nanmean(mXppant); %mean overall (remove b/w sub differences
    
    %for each observation, subjtract the subj average, add
    %the group average.
    NEWdata1 = m1x - mXppant + repmat(mXgroup, length(m1x),1);
    NEWdata2 = m2x - mXppant + repmat(mXgroup, length(m2x),1);
    
    %             %compute new stErr %which version?
    stE1 = nanstd(NEWdata1)/sqrt(size(m1x,1));
    stE2 = nanstd(NEWdata2)/sqrt(size(m2x,1));
    
    std1=nanstd(ppantFIRSTSIG(1,itimezero,:))/ sqrt(21);
    std2=nanstd(ppantFIRSTSIG(2,itimezero,:))/ sqrt(21);
    
    plot([m1 m1], [.4 .4], 'marker', '*', 'linew', 400, 'color', 'r');
    plot([m2 m2], [.4 .4], 'marker', '*', 'linew', 400, 'color', 'k');
    hold on; %place errbars
    plot([m1-stE1/2 m1+stE1/2], [.4 .4], 'marker', '.', 'linew', 1, 'color', 'r');
    plot([m2-stE2/2 m2+stE2/2], [.4 .4], 'marker', '.', 'linew', 1, 'color', ['k']);
    xlim([-2.5 2.5]);
    xlabel('Time from PFI decrease')
    set(gca, 'fontsize', 25, 'ytick', [])
    set(gcf, 'color', 'w')
    %%
    
    %     legend(pl, {'PFIinc 20', 'PFIinc 40', 'PFIdec 20', 'PFIdec 40'})
    %     histogram(ppantFIRSTSIG(2,1,:));
    shg
    %%
end

if job.BPandSSVEPtimecourseacrossppants_numSeparate==1
    %%
    %%
    getelocs
    rmvbase=0;
    
    checkcluster=1;
    clf
    plcount=1;
    %
    for hzis=1%:2
        
        cd(basefol)
        cd('newplots-MD')
        switch hzis
            case 1
                loadname='GFX_PFIperformance_withSNR_20';
                
                lint='-';
                usehz=20;
                sigheight=4.1;
            case 2
                
                loadname='GFX_PFIperformance_withSNR_40';
                
                usehz=40;
                lint=':';
                sigheight=3.9;
        end
        
        %%
        for itimezero=1%2
            anovatestdata=[];
            
            for idur=1%:2
                cd(basefol)
                cd('newplots-MD')
                switch  idur
                    case 1
                        load([loadname '_min0'])
                        
                    case 2
                        load([loadname '_min30'])
                    case 3
                        load([loadname '_min60'])
                end
                
                tbase = tgrm-3;
                useSNRNUM=[];
                switch itimezero
                    case 1
                        
                        useSNRNUM(1,:,:)=squeeze(storeacrossPpant_onsetSNR_PFI1);
                        useSNRNUM(2,:,:)=squeeze(storeacrossPpant_onsetSNR_PFI2);
                        useSNRNUM(3,:,:)=squeeze(storeacrossPpant_onsetSNR_PFI3);
                        
                        chtype='PFI increase';
                        %                     chtype='button press';
                        
                        col=[0 .5 0];
                        
                    case 2
                        
                        
                        useSNRNUM(1,:,:)=squeeze(storeacrossPpant_offsetSNR_PFI1);
                        useSNRNUM(2,:,:)=squeeze(storeacrossPpant_offsetSNR_PFI2);
                        useSNRNUM(3,:,:)=squeeze(storeacrossPpant_offsetSNR_PFI3);
                        
                        chtype='PFI decrease';
                        %                     chtype='button release';
                        
                        
                        col='r';
                        
                        
                end
                
                
                
                %
                %% SNR
                figure(1);
                for inum=1:2%:size(useSNRNUM,1)
                    ppantMeanSNR=squeeze(useSNRNUM(inum,:,:));
                    
                    plis = 3;%+(itimezero-1);
                    %             subplot(2,2,plis)
                    hold on
                    %plot across ppant trace
                    
                    
                    %adjust standard error as per COusineau(2005)
                    %confidence interval for within subj designs.
                    % y = x - mXsub + mXGroup,
                    x = ppantMeanSNR;
                    
                    mXppant =squeeze( nanmean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
                    mXgroup = nanmean(nanmean(mXppant)); %mean overall (remove b/w sub differences
                    
                    %for each observation, subjtract the subj average, add
                    %the group average.
                    NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
                    
                    %             %compute new stErr %which version?
                    stE = nanstd(NEWdata)/sqrt(size(x,1));
                    %
                    sh=shadedErrorBar(tgrm-3, nanmean(ppantMeanSNR,1),stE,[lint],[1]);
                    
                    sh.mainLine.Color = col;
                    sh.mainLine.LineWidth = inum*2;
                    sh.patch.FaceColor = col;
                    sh.edge(1).Color = col;
                    sh.edge(2).Color = col;
                    
                    if hzis==2
                        sh.mainLine.Color= 'k';
                    end
                    set(gca, 'fontsize', 15)
                    hold on;
                    %                 %%
                    %
                    ylabel({['SNR(dB/Hz)'];['baseline corrected']})
                    
                    %             title({[num2str(usehz) ' Hz SSVEP']})
                    xlabel(['Time from ' chtype])
                    set(gca, 'fontsize', 25)
                    
                    set(gcf, 'color', 'w')
                    anovatestdata(inum,1:size(ppantMeanSNR,1),:) = ppantMeanSNR;
                    legendprint(inum)=sh.mainLine;
                    plcount=plcount+1;
                end
            end
        end
        %%
        %check for sig
        %use RMANOVA.
        %also calculate and plot significance, using ANOVA across
        %conditions.
        factorarray = [ones(21,1);repmat(2,21,1);repmat(3,21,1)];
        ss=[1:21]';
        subjects=[ss;ss;ss];
        
        %store ANOVA results
        pvals=nan(1,size(anovatestdata,3));
        Fvals=pvals;
        %store planned comparisons; for AttL vs nAttL
        %                 presultPC=nan(1,size(anovatestdata,3));
        
        for itime=2:size(anovatestdata,3) %at each timeponit
            
            %arrange for rmanova. single cols.
            %                     rmanova(data,factorarray,subjects,varnames, btw_ss_col)
            %                     datat= [squeeze(anovatestdata(1,:,itime))'; squeeze(anovatestdata(2,:,itime))';...
            %                         squeeze(anovatestdata(3,:,itime))'];
            %
            %                     anvret=rmanova(datat,factorarray,subjects);
            %         pvals(itime)=anvret.p(1);
            %         Fvals(itime)=anvret.table{2,6};
            
            
            
            [~, pvals(itime),~,stat]=ttest(squeeze(anovatestdata(1,:,itime)),squeeze(anovatestdata(2,:,itime)), .05);
            Fvals(itime)=stat.tstat;
        end
        %         %%
        
        sigs=find(pvals<.05);
        %
        %             %perform cluster based correction.
        if length(sigs)>1 && checkcluster==1
            % find biggest cluster:
            %finds adjacent time points
            vect1 = diff(sigs);
            v1 = (vect1(:)==1);
            d = diff(v1);
            clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
            %grab largest
            %                 ignore bad points.
            
            
            % find biggest cluster:
            %finds adjacent time points
            sigs = find(pvals<.05);
            
            vect1 = diff(sigs);
            v1 = (vect1(:)==1);
            d = diff(v1);
            clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
            [~,maxClust] = max(clusterSTandEND(:,2)-clusterSTandEND(:,1));
            
            %%
            for icl=1%:size(clusterSTandEND,1)
                
                %start and end are now:
                STC=sigs(clusterSTandEND(icl,1));
                ENDC=sigs(clusterSTandEND(icl,2)+1);
                checktimes =STC:ENDC;
                observedCV = sum(abs(Fvals(checktimes)));
                % now shuffle condition labels to see if this cluster is
                % sig (compared to chance).
                sumTestStatsShuff = zeros(1,2000);
                
                for irand = 1:2000
                    %testing the null that it isn't mismatched - matched at time 2
                    % which creates a diff. so select from either!
                    shD= zeros(size(anovatestdata,1),length(checktimes),size(anovatestdata,2));
                    
                    %change shuffle parameters based on test of
                    %interest. (ie between conditions, or temporal
                    %null).
                    
                    for ipartition = 1:size(anovatestdata,1)
                        for ippant = 1:size(anovatestdata,2)
                            for itime=1:length(checktimes)
                                
                                if mod(randi(100),2)==0 %if random even number
                                    pdata = anovatestdata(1,randi(size(anovatestdata,2)), checktimes(itime)); %select both chans
                                else
                                    pdata = anovatestdata(2,randi(size(anovatestdata,2)), checktimes(itime));
                                end
                                
                                %                                         pdata = anovatestdata(randi(3),randi(size(anovatestdata,2)), checktimes(itime)); %select both chans
                                
                                shD(ipartition,itime,ippant) = pdata;
                            end
                        end
                    end
                    
                    %now compute difference between out hypothetical topoplots,
                    % and test for sig, checking the accumulated test statistic at our
                    % times of interest
                    tvalspertimepoint = zeros(1,length(checktimes));
                    
                    testdata = squeeze(shD(1,:,:)) - squeeze(shD(2,:,:));
                    pvalsNew=[];
                    FvalsNew=[];
                    for itest = 1:length(checktimes) %test each time point
                        
                        %
                        %                             datat= [squeeze(shD(1,:,itest))'; squeeze(shD(2,:,itest))';...
                        %                         squeeze(shD(3,:,itest))'];
                        %
                        %                     anvret=rmanova(datat,factorarray,subjects);
                        % %                     pvalsNew(itime)=anvret.p(1);
                        %                     FvalsNew(itime)=anvret.table{2,6};
                        
                        %
                        %                             tvalspertimepoint(1,itest) = anvret.table{2,6};
                        
                        
                        [~, p, ~,stat]= ttest(testdata(itest,:));
                        
                        tvalspertimepoint(1,itest) = stat.tstat;
                        
                    end
                    
                    sumTestStatsShuff(1,irand) = sum(abs(tvalspertimepoint));
                end %repeat nshuff times
                
                
                %is the observed greater than CV?
                % plot histogram:
                figure(2);
                
                clf
                
                %                     subplot(2,1, plcount-1)
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
                    
                    %observed pvalue in distribution=
                    [~, c2] = min(abs(H.Data-observedCV)); %closest value.
                    pvalis= 1-cdf(c2);
                    title(['sum tvals = ' num2str(observedCV), 'p=' num2str(pvalis)]);
                    
                    %              title('Spatial Cluster  Significant!')
                    for itime=checktimes
                        figure(1);
                        hold on
                        plot(tgrm(itime)-3, sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'm')
                    end
                end
                
                %
            end
        else
            try %just plotem
                for itime=sigs
                    figure(1);
                    hold on
                    %                             plot(tgrm(itime)-3, sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'k')
                end
            catch
            end
            
        end
        %
        %
        %%
    end
    %
    figure(1);
    hold on;
    
    axis tight
    ylim([-4 4.5])
    hold on
    plot(xlim, [0 0], ['k:']);
    plot([0 0] , ylim, ['k:'], 'linewidth', 1)
    %plot patch:
    title([num2str(usehz) ' Hz SSVEP'])
    
    shg
    
    hold on
    
    %
    
    
    %
    
    lg=legend(legendprint, {'PFI<median', 'PFI>median', '3 Targets PFI'});
    if itimezero==1
        
        set(lg, 'location', 'NorthWest')
        
    else
        
        set(lg, 'location', 'SouthWest')
    end
    set(lg, 'location', 'SouthWest')
    %
    cd(basefol)
    cd('newplots-MD')
    
    %
    cd('Figures')
    print('-dpng', ['PFI trace ' num2str(usehz) 'Hz SSVEP summary, during ' num2str(chtype) '_sepNumberPFI.png'])
end

if job.BPandSSVEPtimecourseacrossppants_dursSeparate==1
    %%
    getelocs
    rmvbase=0;
    
    
    
    clf
    plcount=1;
    for hzis=2%1:2
        
        cd(basefol)
        cd('newplots-MD')
        switch hzis
            case 1
                loadname='GFX_PFIperformance_withSNR_20';
                
                lint='-';
                usehz=20;
                sigheight=4.1;
            case 2
                
                loadname='GFX_PFIperformance_withSNR_40';
                
                usehz=40;
                lint=':';
                sigheight=3.9;
        end
        
        %%
        for itimezero=1%2
            ttestdata=[];
            for idur=1:2
                cd(basefol)
                cd('newplots-MD')
                switch  idur
                    case 1
                        load([loadname '_min60under'])
                        
                    case 2
                        load([loadname '_min60'])
                    case 3
                        load([loadname '_min120'])
                end
                
                tbase = tgrm-3;
                switch itimezero
                    case 1
                        useBP=storeacrossPpant_onsetBP;
                        useSNR=storeacrossPpant_onsetSNR;
                        
                        useTOPO=storeacrossPpant_onsetTOPO;
                        
                        chtype='PFI increase';
                        %                     chtype='button press';
                        
                        col=[0 .5 0];
                        
                    case 2
                        
                        useBP=storeacrossPpant_offsetBP;
                        useSNR=storeacrossPpant_offsetSNR;
                        
                        useTOPO=storeacrossPpant_offsetTOPO;
                        
                        
                        chtype='PFI decrease';
                        %                     chtype='button release';
                        
                        
                        col='r';
                        
                        
                end
                
                
                %%
                %             plis = 1;%+(itimezero-1);
                %             subplot(2,2,plis)
                %             ppantMeanBP= squeeze(mean(useBP,1));
                %             %adjust standard error as per COusineau(2005)
                %             %confidence interval for within subj designs.
                %             % y = x - mXsub + mXGroup,
                %             x = squeeze(mean(useBP,2));
                %
                %             mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
                %             mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
                %
                %             %for each observation, subjtract the subj average, add
                %             %the group average.
                %             NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
                %
                %             %compute new stErr %which version?
                %             stE = std(NEWdata)/sqrt(size(x,1));
                %             sh=shadedErrorBar([-3:1/60:3], mean(ppantMeanBP,1),stE, [col '-'], [1]);
                %             sh.mainLine.LineWidth=3;
                %
                %             set(gca, 'fontsize', 25)
                %             hold on;
                %             ylabel('Buttons pressed ')
                %             axis tight
                % %             ylim([0 2])
                %             plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                %             %
                %             xlabel(['Time from ' chtype])
                %
                %SNR
                figure(1);
                plis = 3;%+(itimezero-1);
                %             subplot(2,2,plis)
                hold on
                %plot across ppant trace
                ppantMeanSNR= squeeze(mean(useSNR,2));
                
                %adjust standard error as per COusineau(2005)
                %confidence interval for within subj designs.
                % y = x - mXsub + mXGroup,
                x = ppantMeanSNR;
                
                mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
                mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
                
                %for each observation, subjtract the subj average, add
                %the group average.
                NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
                
                %             %compute new stErr %which version?
                stE = std(NEWdata)/sqrt(size(x,1));
                %
                sh=shadedErrorBar(tgrm-3, mean(ppantMeanSNR,1),stE,[lint],[1]);
                
                sh.mainLine.Color = col;
                sh.mainLine.LineWidth = idur*2;
                sh.patch.FaceColor = col;
                sh.edge(1).Color = col;
                sh.edge(2).Color = col;
                
                if hzis==2
                    sh.mainLine.Color= 'k';
                end
                set(gca, 'fontsize', 15)
                hold on;
                %                 %%
                %
                ylabel({['SNR(dB/Hz)'];['baseline corrected']})
                
                %             title({[num2str(usehz) ' Hz SSVEP']})
                xlabel(['Time from ' chtype])
                set(gca, 'fontsize', 25)
                
                set(gcf, 'color', 'w')
                ttestdata(idur,1:size(ppantMeanSNR,1),:) = ppantMeanSNR;
                legendprint(idur)=sh.mainLine;
                plcount=plcount+1;
            end
        end
        %%
        %check for sig
        %use RMANOVA.
        
        
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
            
            
            % find biggest cluster:
            %finds adjacent time points
            sigs = find(pvals<.05);
            
            vect1 = diff(sigs);
            v1 = (vect1(:)==1);
            d = diff(v1);
            clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
            [~,maxClust] = max(clusterSTandEND(:,2)-clusterSTandEND(:,1));
            
            %%
            for icl=1:size(clusterSTandEND,1)
                
                %start and end are now:
                STC=sigs(clusterSTandEND(icl,1));
                ENDC=sigs(clusterSTandEND(icl,2)+1);
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
                            for ippant = 1:size(ttestdata,2)
                                for itime=1:length(checktimes)
                                    
                                    if mod(randi(100),2)==0 %if random even number
                                        pdata = ttestdata(1,randi(size(ttestdata,2)), checktimes(itime)); %select both chans
                                    else
                                        pdata = ttestdata(2,randi(size(ttestdata,2)), checktimes(itime));
                                    end
                                    
                                    shD(ipartition,itime,ippant) = pdata;
                                end
                            end
                        end
                    else %null is that there are no temporal coincident sig values.
                        for ipartition = 1:2
                            for ippant = 1:size(ttestdata,2)
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
                
                clf
                
                %                     subplot(2,1, plcount-1)
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
                        plot(tgrm(itime)-3, sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'm')
                    end
                end
                
                %
            end
        end
        
        
        %%
    end
    %%
    axis tight
    ylim([-4 4.5])
    hold on
    plot(xlim, [0 0], ['k:']);
    plot([0 0] , ylim, ['k:'], 'linewidth', 1)
    %plot patch:
    title([num2str(usehz) ' Hz SSVEP'])
    
    shg
    
    hold on
    
    %%
    
    
    %%
    
    lg=legend(legendprint, {'all PFI', 'PFI>1 second'});
    if itimezero==1
        set(lg, 'location', 'NorthWest')
    else
        set(lg, 'location', 'SouthWest')
    end
    cd(basefol)
    cd('newplots-MD')
    
    
    cd('Figures')
    print('-dpng', ['Catch trace ' num2str(usehz) 'Hz SSVEP summary, during ' num2str(chtype) '_sepDurations.png'])
end
