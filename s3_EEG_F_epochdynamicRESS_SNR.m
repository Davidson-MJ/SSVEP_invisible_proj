% Follows the format of s3_CE... only now applying to RESS timeseries

%having epoched tgs and catch, POST RESS. apply FFT and SNR to windows when
%all targets are present/ away from Catch stimulus onset.

clear all
cd('/Users/mDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA')

basefol=pwd;
clearvars -except basefol allppants
dbstop if error
%%

%%
%separate analysis at image level (trial dimension on Y axis).
% for PFI. First need to check individual variation in nPFI:
%this job sorts all participant PFI, by number of buttons pressed.
% Also prints the behaviour (button-press), and SNR side by side.
job.erpimagePpantlevel=0;

%Resamples along Y-dimension, to average across participants. 
%Concatenates across participants, the results of the above.
job.concaterpdataacrossppants=0; % concat above


%Same for catch/PMD
job.erpimage_catchPpantlevel=0; %stores responses as 100 ydimension trial x trial
job.concaterpdataforcatchacrossppants=0;


%%%%% %%%%%%% %%%%%% %%%%%% %%%%%%
%%%%% %%%%%%% %%%%%% %%%%%% %%%%%%
% PARAMETERS FOR ALL SPECTROGRAM! KEEP CONSISTENT
%%%%% %%%%%%% %%%%%% %%%%%% %%%%%%
%%%%% %%%%%%% %%%%%% %%%%%% %%%%%%
param_spcgrm.tapers = [1 1];
param_spcgrm.Fs= [250];
param_spcgrm.Fpass= [0 50];
param_spcgrm.trialave=0;
param_spcgrm.pad = 2;

movingwin=[2.5,.25]; %

%%%%% %%%%%%% %%%%%% %%%%%% %%%%%%
%%%%% %%%%%%% %%%%%% %%%%%% %%%%%%
% what size will the output be? (taken from inside mtspecgramc)
%frequency:
Nwin=round(250*movingwin(1)); % number of samples in window
Nstep=round(movingwin(2)*250); % number of samples to step through
nfft=max(2^(nextpow2(Nwin)+param_spcgrm.pad),Nwin);
f=getfgrid(250,nfft,[param_spcgrm.Fpass]);
% use to automatically preallocate variable size.
Nf=length(f);
%timesteps:
Nwin=round(250*movingwin(1)); % number of samples in window
Nstep=round(movingwin(2)*250); % number of samples to step through
N=1501; %size EEG epochs
winstart=1:Nstep:N-Nwin+1;
winmid=winstart+round(Nwin/2);
% use to automatically preallocate variable size.
Nt=winmid/250;


% begin.
getelocs;

cd(basefol)
cd('EEG')

pdirs = dir([pwd filesep '*EEG']);

%% load data to determine physical catch timing.

allppants=[1,2,4,6,7,9:19]; % NEW ppants.
window=[-3 3];

srate=250;

epochdur = sum(abs(window));

timeidDYN = [0:1/srate:epochdur];
timeidDYN= timeidDYN-3;
%%
onsetc = ceil(epochdur)/2;
% peakfreqsare=[20,40]; %hz
%timing
tt = 0:1/srate:60;

% which frequencies to analyze?/apply?
peakfreqsare=[15,20,30, 40, 45, 60, 5, 25, 35 ]; % don't change!

%% begin

if job.erpimagePpantlevel==1
    %%
    rmvbase=0;
    useSNRorhilbert=1;
    
    for ifol = allppants        
        for ihz=[1,2] % in peakfreqsare.
            
            usehz= peakfreqsare(ihz);
            
            icounter=1;
            
            cd(basefol)
            cd('EEG')
            cd(pdirs(ifol).name)
            load(['ppant_PFI_Epoched_RESS'])
            
            for itimezero = 1:2
                if itimezero==1
                    
                    %append them all. TARG-> Disappearing (more buttons
                    %pressed).
                    
                    datatouse = cat(1, ress_PFI_0_1_Hz(ihz).ressTS,ress_PFI_1_2_Hz(ihz).ressTS,ress_PFI_2_3_Hz(ihz).ressTS, ress_PFI_3_4_Hz(ihz).ressTS);
                    %                     RTstouse = [durs0_1'; durs1_2'; durs2_3'];
                    BPstouse = cat(1, ress_PFI_0_1_Hz(ihz).BPs,ress_PFI_1_2_Hz(ihz).BPs,ress_PFI_2_3_Hz(ihz).BPs,ress_PFI_3_4_Hz(ihz).BPs);
                    ctype = 'PFI increase';
                    bsrem = [-3 -1]; %seconds
                    
                else
                    
                    datatouse = cat(1, ress_PFI_1_0_Hz(ihz).ressTS,ress_PFI_2_1_Hz(ihz).ressTS,ress_PFI_3_2_Hz(ihz).ressTS,ress_PFI_4_3_Hz(ihz).ressTS);
                    %                     RTstouse = [durs1_0'; durs2_1'; durs3_2'];
                    BPstouse= cat(1, ress_PFI_1_0_Hz(ihz).BPs,ress_PFI_2_1_Hz(ihz).BPs,ress_PFI_3_2_Hz(ihz).BPs,ress_PFI_4_3_Hz(ihz).BPs);
                    ctype = 'PFI decrease';
                    
                    
                    bsrem = [1 3]; %seconds
                    
                    
                end

                %plot the topo for sanity check: (pre disap)
                windcheck= [-3 -0.1];
                tidx=dsearchn(timeidDYN', [windcheck]');
                               
                allt= 1:size(datatouse,1);
               
                clf
                colormap('viridis')
                %plot BP for comparison
                subplot(3,2, [1 3])
           
                %sort by accum PFI in window after onset.
                if itimezero==1
                    accumPeriod = sum(BPstouse(:,180:end),2);
                else % or before
                    accumPeriod = sum(BPstouse(:,1:180),2);
                    
                end
                
                [checkBP, cid] = sort(accumPeriod, 'descend');
                
                
                % MOVE nan to bottom of order (not top)
                nanind= find(isnan(checkBP));
                %new start point
                cid=[cid(length(nanind)+1:end) ; cid(nanind)];
                %                     checkBP=[checkBP(length(nanind)+1:end); checkBP(nanind)]
                
                %rearrange.
                BPstouse=BPstouse(cid,:);
                datais=squeeze(datatouse(cid,:));
                
                %%
                
                tBP=-3:1/60:3;
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
                %%
                % show mean over time.
                subplot(3,2,5)
                plot(tBP, nanmean(BPstouse,1),['k-'], 'linewidth', 3)
                %         shadedErrorBar([-3:1/60:3], nanmean(acrBP,1),nanstd(acrBP,1))
                set(gca, 'fontsize', 15)
                hold on;
                ylabel('nanmean BP ')
                %
                xlabel('Time secs')
                axis tight
                plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                
                
                set(gca, 'fontsize', 15)
                %%
                  %take moving window spectrogram
                    [sgrm ,tgrm, fgrm] = mtspecgramc(datais', movingwin, param_spcgrm);
                    
                    
                    %%
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
                        %                             hold on; plot(kernelw)
                    
                    %%
                    
                    %store just stim freq.
                    [~, hzid]= min(abs(fgrm-usehz));
                    %
                    % %                             %reduce size.
                    dynSSEP=squeeze(snr_sgrm(:,hzid,:))';
                    % %
                    timeidDYN=tgrm-3;
                   
                snrgrm20=dynSSEP;
                % smooth across trials ( ydimension)
                sm_snrgrm20=zeros(size(snrgrm20));
                for itime=1:size(snrgrm20,2)
                    sm_snrgrm20(:,itime)= smooth(snrgrm20(:,itime),5);
                end
                %%
                hold on
                subplot(3,2,[2 4]);
                hold on
                imagesc(timeidDYN,  1:size(dynSSEP,1), sm_snrgrm20);
                c=colorbar;               
                    ylabel(c, 'log(SNR)')
               
                xlabel('Time secs')
                caxis([-1*max(max(sm_snrgrm20)) max(max(sm_snrgrm20))])
                caxis([0  max(max(sm_snrgrm20))-1])%                 set(gca, 'ytick', 1:24, 'yticklabel', round(sortedRTs./60,2))
                % caxis([ 1 5])
                axis tight
                hold on
                plot([0 0] , ylim, ['w:'], 'linewidth', 3)
                title({['sorted ' num2str(usehz) 'Hz RESS data (smoothed) for  ' ctype ','];['ppant' num2str(ifol)]})
                set(gca, 'fontsize', 15, 'yticklabel', [])
                
                
                %% show mean over time.
                subplot(3,2,6)
                plot(timeidDYN, nanmean(snrgrm20,1),['k-'], 'linewidth', 3)
                %         shadedErrorBar([-3:1/60:3], nanmean(acrBP,1),nanstd(acrBP,1))
                set(gca, 'fontsize', 15)
                %                 ylim([1 3])
                hold on;
                ylabel('RESS logSNR')
                %%
                axis tight
                plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                shg
                xlabel('Time secs')
                set(gcf,'color', 'w');
                
                %%
                shg
                cd([basefol filesep 'Figures' filesep 'Participant summaries'])
                %%
                %%
                print('-dpng', ['figure_PFI_' ctype '_summary_BPandSSVEP_' num2str(usehz) '_RESS_subj' num2str(ifol) '_longerwindow.png']);
                
                
                switch itimezero
                    case 1 %store for across ppant plots:
                        Ppant_onsetBP=BPstouse;
                        Ppant_onsetSNR=snrgrm20; %sorted.
                       
                    case 2                        
                        Ppant_offsetBP=BPstouse;
                        Ppant_offsetSNR=snrgrm20; %always sorted in descending order of PFI.
                        
                end
            end
            
            %%
%             
            savename=['PFIperformance_withSNR_' num2str(usehz) '_RESS'];
            
            cd(basefol)
            cd('EEG')
            cd(pdirs(ifol).name)
            %%
            save(savename,...
                'Ppant_onsetBP','Ppant_offsetBP',...
                'Ppant_onsetSNR', 'Ppant_offsetSNR', 'timeidDYN', 'param_spcgrm', 'kernelw')
%             
            
        end
    end
end


%% This is the resampling along y-dimension, for comparion of nTargets by PFI/PMD

if job.concaterpdataacrossppants==1
    
    dursMINIMUM =0; % all PFI included.
    
    for ihz=1:2
        usehz=peakfreqsare(ihz);
        
        loadname=['PFIperformance_withSNR_' num2str(usehz) '_RESS'];
        
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
            cd([pdirs(ippant).name])
            
            %onset types
            load(loadname)
            
            % we need to resample the BP and SNR data, to equate across
            % trial types    (ignoring nan)
            
            % currently at 100 fs.
            
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
                
                % first find out if there are nans, and remove.
                if any(isnan(mean(pData,2)))
                    nantr= find(isnan(mean(pData,2)));
                    % for the trials with a nan. we need to
                    % remove.
                    pData(nantr,:)=[];
                end
                lc=size(pData,1);
                
                %discount the nans in data. (since these are
                %sorted, should be at bottom).
                %%
                pData= pData(1:lc,:);
                
                
                
                
                %then resample to 100 trials on the y dimension:
                
                pdata_final= resample(pData,100,size(pData,1));
                
                
                switch ionoff
                    case 1
                        Ppant_onsetSNR = pdata_final;
                    case 2
                        Ppant_offsetSNR = pdata_final;
                    case 3
                        Ppant_onsetBP = pdata_final;
                    case 4
                        Ppant_offsetBP = pdata_final;
                        
                end
                
                
                
            end
            
            storeacrossPpant_onsetBP(icounter,:,:)=Ppant_onsetBP;
            storeacrossPpant_onsetSNR(icounter,:,:)=Ppant_onsetSNR;
            
            storeacrossPpant_offsetBP(icounter,:,:)=Ppant_offsetBP;
            storeacrossPpant_offsetSNR(icounter,:,:)=Ppant_offsetSNR;
            
            
            icounter=icounter+1;
            
        end
       
        %% save appropriately
        cd(basefol)
        cd('EEG')
        cd('GFX_EEG-RESS')
        
        savename=['GFX_PFIperformance_withSNR_' num2str(usehz) '_min' num2str(dursMINIMUM) '_RESS'];
        %%
        save(savename,...
            'storeacrossPpant_onsetBP','storeacrossPpant_offsetBP',...
            'storeacrossPpant_onsetSNR', 'storeacrossPpant_offsetSNR', 'timeidDYN')
        
        
        
    end
    
end
%%

if job.erpimage_catchPpantlevel==1
    %% as above, but adapted for PMD data.
    
    for ifol = allppants
        for ihz=[1:2]
            usehz= peakfreqsare(ihz);
            
            icounter=1;
            
            cd(basefol)
            cd('EEG')
            cd(pdirs(ifol).name)
            load(['ppant_Catch_Epoched_RESS'])
            load('ppant_Catch_Epoched'); %for BP
            clearvars ppant_SNREEG*
            
            for itype=1:2
             switch itype
                    case 1 % BP to catch
                        
                        switch ihz
                            case {1,3,5} % TGs
                                %adjust to 1,2,3
                                ind = ceil(ihz/2);
                                datatouse = ress_BPcatchonsetTGs(ind,:,:);
                            case {2,4,6} % BGs
                                ind=ihz/2;
                                datatouse = ress_BPcatchonsetBGs(ind,:,:);
                            case 7
                                datatouse = ress_BPcatchonsetIMs(1,:,:);
                        end
                        
                        BPstouse = withincatchonsetBPs;
                        ctype = 'within catch onset BP';
                        
                    case 2
                        
                        switch ihz
                            case {1,3,5} % TGs
                                %adjust to 1,2,3
                                ind = ceil(ihz/2);
                                datatouse = ress_BPcatchoffsetTGs(ind,:,:);
                            case {2,4,6} % BGs
                                ind=ihz/2;
                                datatouse = ress_BPcatchoffsetBGs(ind,:,:);
                            case 7
                                datatouse = ress_BPcatchoffsetIMs(1,:,:);
                        end
                        
                        BPstouse = postcatchoffsetBPs;
                        ctype = 'post offset BP';
                       
                end
                datatouse=squeeze(datatouse);
                
                
                %plot the topo for sanity check: (pre disap)
                windcheck= [-3 -0.1];
                tidx=dsearchn(timeidDYN', [windcheck]');
               
                
                %restrict to certain PFI duration?
                allt= 1:size(datatouse,1);
               
              
                %%
              
                
                %sort by accum PFI in window after onset.
                switch itype
                    case {1,3,5}
                        accumPeriod = sum(BPstouse(:,180:end),2);
                    case {2,4}
                        accumPeriod = sum(BPstouse(:,1:180),2);
                        
                end
                
                [checkBP, cid] = sort(accumPeriod, 'descend');
                
                
                % MOVE nan to bottom of order (not top)
                nanind= find(isnan(checkBP));
                %new start point
                cid=[cid(length(nanind)+1:end) ];%  ; cid(nanind)];
                
                
                %rearrange.
                BPstouse=BPstouse(cid,:);
               
                datais=squeeze(datatouse(cid,:));
                %                 RTstouse=RTstouse(cid);
              
                %% moving window spectrogram.
                
                    [sgrm ,tgrm, fgrm] = mtspecgramc(datais', movingwin, param_spcgrm);
                    %% %% Two SNR methods:
                    
                    %%
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
                        
                        %                              %% sanity check:
                        %                               figure(2);
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
                        %                             hold on; plot(kernelw)
                    %%
                    
                    %store just stim freq.
                    [~, hzid]= min(abs(fgrm-usehz));
                    %
                    % %                             %reduce size.
                    dynSSEP=squeeze(snr_sgrm(:,hzid,:))';
                    % %
                    timeidDYN=tgrm-3;
                
                %%\
                snrgrm20=dynSSEP;
                % smooth across trials
                sm_snrgrm20=zeros(size(snrgrm20));
                for itime=1:size(snrgrm20,2)
                    sm_snrgrm20(:,itime)= smooth(snrgrm20(:,itime),5);
                end
                %% plot both
                subplot(3,2, [1 3])
           
                %%
                
                tBP=-3:1/60:3;
                imagesc(tBP,1:length(cid),BPstouse)
                c=colorbar;
                ylabel(c, 'buttons pressed')
                %                 set(gca, 'ytick', 1:length(cid), 'yticklabel', round(sortedRTs./60,2))
%                 ylabel('PFI events')
                xlabel('Time [sec]')
                hold on
                plot([0 0] , ylim, ['w:'], 'linewidth', 3)
                title({['sorted BP data for ' ctype ','];['ppant' num2str(ifol)]})
                set(gca, 'fontsize', 15)
                
                hold on
                subplot(3,2,[2 4]);
                hold on
                imagesc(timeidDYN,  1:size(dynSSEP,1), sm_snrgrm20);
                c=colorbar;                
                ylabel(c, 'log(SNR)')                
                xlabel('Time secs')
                
                try caxis([0  max(max(sm_snrgrm20))-1])%
                catch
                    caxis([0  max(max(sm_snrgrm20))])%
                end
               
                axis tight
                hold on
                plot([0 0] , ylim, ['w:'], 'linewidth', 3)
                title({['sorted ' num2str(usehz) 'Hz RESS data (smoothed) for  ' ctype ','];['ppant' num2str(ifol)]})
                set(gca, 'fontsize', 15, 'yticklabel', [])
               %%
                %% show mean over time.
                 subplot(3,2,5)
                plot(tBP, nanmean(BPstouse,1),['k-'], 'linewidth', 3)
                set(gca, 'fontsize', 15)               
                
                subplot(3,2,6)
                plot(timeidDYN, nanmean(snrgrm20,1),['k-'], 'linewidth', 3)                
                set(gca, 'fontsize', 15)               
                hold on;
                ylabel('RESS logSNR')
                %%
                axis tight
                plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                shg
                xlabel('Time secs')
                set(gcf,'color', 'w');
               
                %%
                shg
                cd([basefol filesep 'Figures' filesep 'Participant Catch summaries'])
                %%
                %%
                print('-dpng', ['figure_' ctype '_summary_BPandSSVEP_' num2str(usehz) '_RESS_subj' num2str(ifol) ' snr method ' num2str(snrmethod) '.png']);
                
                
                switch itype
                    case 1
                        BPcatch_onsetSNR_RESS_BPs=BPstouse;
                        BPcatch_onsetSNR_RESS=snrgrm20; %sorted (un smoothed)
                    case 2
                        BPcatch_offsetSNR_RESS_BPs=BPstouse;
                        BPcatch_offsetSNR_RESS=snrgrm20; %sorted (un smoothed)
                    
                end
            end
            
            
            %%
            
            savename=['Catch_performance_withSNR_' num2str(usehz) '_RESS'];
            
            cd(basefol)
            cd('EEG')
            cd(pdirs(ifol).name)
            %%
            save(savename,...
                'BPcatch_onsetSNR_RESS_BPs','BPcatch_onsetSNR_RESS',...
                'BPcatch_offsetSNR_RESS_BPs','BPcatch_offsetSNR_RESS',...                
                'timeidDYN', 'param_spcgrm', 'kernelw')
            
            
        end
    end
end
%%

if job.concaterpdataforcatchacrossppants==1
    %%
    
    for ihz=1:2
        usehz=peakfreqsare(ihz);
        
        loadname=['Catch_performance_withSNR_' num2str(usehz) '_RESS'];
        
        storeacrossPpant_catchEVENTS_BPs = zeros(length(allppants), 5, 132, 361);
        storeacrossPpant_catchEVENTS_SNR = zeros(length(allppants), 5, 132, length(Nt));
        icounter=1;
        
        for ippant = allppants
            cd(basefol)
            cd('EEG')
            cd([pdirs(ippant).name])
            
            %onset types
            
            load(loadname)
            
            % we need to resample the BP and SNR data, to equate across
            % trial types    (ignoring nan)
            
            % currently at 100 fs.
            
            % beware edge artefacts!
            
            % also apply participant level smoothing.
           
                
                for ionoff=1:4 %all SNR and BP combos.
                    
                    switch ionoff
                     
                        case 1
                            pData = BPcatch_onsetSNR_RESS;
                            storeid=3;
                        case 2
                            pData = BPcatch_onsetSNR_RESS_BPs;
                        case 3
                            pData = BPcatch_offsetSNR_RESS;                            
                            storeid=4;
                        case 4
                            pData = BPcatch_offsetSNR_RESS_BPs;
                       
                    end
                    
                    % first find out if there are nans, and remove.
                    if any(isnan(mean(pData,2)))
                        nantr= find(isnan(mean(pData,2)));
                        % for the trials with a nan. we need to
                        % remove.
                        pData(nantr,:)=[];
                    end
                    lc=size(pData,1);
                    %discount the nans in data. (since these are
                    %sorted, should be at bottom).
                    %%
                    %%
                    pData= pData(1:lc,:);
                     
                
                %then resample to 100 trials on the y dimension:
                
                pdata_final= resample(pData,100,size(pData,1));
                
                
                switch ionoff
                    case 1
                        Ppant_onsetSNR = pdata_final;
                    case 2
                        Ppant_offsetSNR = pdata_final;
                    case 3
                        Ppant_onsetBP = pdata_final;
                    case 4
                        Ppant_offsetBP = pdata_final;
                        
                end
                
                    if mod(ionoff,2)==0 %even numbers are BP
                        storeacrossPpant_catchEVENTS_BPs(icounter,storeid,:,:)=pdataOUT;
                    else
                        storeacrossPpant_catchEVENTS_SNR(icounter,storeid,:,:)=pdataOUT;
                    end
                end
                
            
            
            
            icounter=icounter+1;
        end
      
        %%
        %save appropriately
        cd(basefol)
        cd('EEG')
        cd('GFX_EEG-RESS')
        
        DIMSare = {'catch onset', 'catch offset', 'BPcatchonset', 'BPcatchoffset', 'invisibleonset'};
        
        savename=['GFX_Catch_performance_withSNR_' num2str(usehz) '_min' num2str(dursMINIMUM) '_RESS'];
        
        
        
        
        
        save(savename,...
            'storeacrossPpant_catchEVENTS_SNR','storeacrossPpant_catchEVENTS_BPs',...
            'DIMSare', 'timeidDYN')
        
        
        
    end
    
end
