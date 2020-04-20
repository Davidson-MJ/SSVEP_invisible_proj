%script 1
% Analyzing  BEH data, modelled on the scripts from previous analysis
% (SSVEP PFI project; Davidson et al., Neuroscience of Consc. 2020).

%MD July 18.
%updated MD 2020

dbstop if error
cd('/Users/MDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA/Behaviour')
basefol=pwd;
%


job.epocharoundCatchonset=0;

job.sortCatchtracesbynumTargets=0; %onset and offsets. appended to 'catch performance' % used in plotPpantProbandSig below.

%The below is labour intensive!, avoid rerunning. This job creates the
%shuffled likelihood button press, per participant, at catch onset and
%offset.
job.calcNEWShuffledCatchOnsetBPProb=0; %!  ! ! ! new ver MD # 04/18.

job.plotCatchtraces_ObsvsShuff_eachppant=0;  % also plots overall RTs to catch onset/offset.
job.plotCatchtraces_ObsvsShuff=1; % plots vs shuffled.

job.calcMissedcatchesbyTargperppant=0; % this appends to Catch Struct -> highlights those potentially missed catches in the all trials plots.


% Note based on all participant trials figures:
% % %
% ppant 3 was very bad, potentially switched button press with release
% % %
% ppant 5 failed many catch responses. Possibly pressed buttons for Targets
% PRESENT instead of disappearaed
% % %
% p 7 had many invisible catches, so do not remove based on catch
% performance.
% ppant 8 failed many catch responses. Not holding button for entire dur.
% % % - rest seem OK at first pass. continue catch analysis below:
%THUS: goodppants=[1,2,4,6,7,9:19]; %



%set up hardcoded dimensions, e.g. xaxis based on epoch length.
%start with a larger window, but trim the output when plotting.
%note that the 'catch start' time recorded, actually begins the ramp, so
%total target absence is catchstart+1.5seconds.
window = [-2 8].*60; %10 second window.
onset=abs(window(1));

%xaxis spacing for plots using this size window.
xtONSET=(window(1):1:window(2))/60;
xtOFFSET=(-1*(window(2)):1:abs(window(1)))/60;

%%
nppants=19;
fontsize=15;


%% BEGIN

%create catch Epoch for BPs
if job.epocharoundCatchonset==1
    load('MD_AllBP_perppant.mat');
    
    for ippant = 1:nppants
        for itrial=1:48
            Details= ppantTrialDatawithDetails(ippant).TrialDetails(itrial);
            tmpData= Details.BPData;
            Catch_any=[]; % can't prealllocate due to differing sizes per participant.
            catchstartTrial =Details.Catchstart_frames;
            if Details.TotalCatchTargetsRemoved>0
                %collect catch BP if it exists.
                
                if Details.TL_Catch==1 %epoch around catch
                    Catch_any= [Catch_any; tmpData(1,catchstartTrial-(abs(window(1))):catchstartTrial+window(2))];
                end
                if Details.TR_Catch==1 %epoch around catch
                    Catch_any= [Catch_any; tmpData(2,catchstartTrial-(abs(window(1))):catchstartTrial+window(2))];
                end
                if Details.BL_Catch==1 %epoch around catch
                    Catch_any= [Catch_any; tmpData(3,catchstartTrial-(abs(window(1))):catchstartTrial+window(2))];
                end
                if Details.BR_Catch==1 %epoch around catch
                    Catch_any= [Catch_any; tmpData(4,catchstartTrial-(abs(window(1))):catchstartTrial+window(2))];
                end
                
                
                %save the long form Epochs, for catch trace BPs
                ppantTrialDatawithDetails(ippant).TrialDetails(itrial).CatchBPs_longepoch = Catch_any; %may not contain 4 rows, and different catch lengths.
                
                %also save JUST the BPs during physical removal (removing pre and post windows), to analyse exact
                %catch accuracy (duration of physical absence recorded).
                ntrain_slim= zeros(size(Catch_any));
                
                physicalabsence = Details.totalCatchdur; %exclude off and on ramp, as not completely absent target.
                start_physicalabsence = abs(window(1)) ; %start in frames
                %record button press for the catch period, excluding the pre
                %catch epoched window,
                %ramp.
                catchexactonly = Catch_any(:,[start_physicalabsence:start_physicalabsence+ physicalabsence]);
                
                %store
                ppantTrialDatawithDetails(ippant).TrialDetails(itrial).CatchBPs_exact= catchexactonly;
                ppantTrialDatawithDetails(ippant).TrialDetails(itrial).truephysicalabsence_duration= physicalabsence;
            end
        end
    end
    %%
    
    save('MD_AllBP_perppant', 'ppantTrialDatawithDetails', '-append')
    %% now stack all in plots, and comparison with previous experiments.
    stackALL=[];
    stackALLinvis=[];
    %%
    for ippant = 1:nppants
        catchProb=[];
        invis_catchProb=[];
        for itrial=1:48
            tmpCatchBPs=ppantTrialDatawithDetails(ippant).TrialDetails(itrial).CatchBPs_longepoch;
            
            %only assess catch timing/Prob, if button wasn't already
            %pressed at catch onset.
            invislocs =[];
            onset=abs(window(1));
            keeptrials=[];
            invis_trials=[];
            for iBP=1:size(tmpCatchBPs,1)
                %for each trial, only retain for analysis of prob if BP
                %not existing
                tmpt=tmpCatchBPs(iBP,:);
                if sum(tmpt(1,onset-60:onset))<60 %60 frames is 1 second.
                    keeptrials=[keeptrials; tmpt];
                else %this catch was for an already invisible object!
                    invis_trials=[invis_trials;tmpt];
                    
                    %also make note of the invisible trials.
                    invislocs = [invislocs, iBP];
                end
                
            end
            if ~isempty(keeptrials)
                catchProb=[catchProb; squeeze(mean(keeptrials,1))]; %average
            end
            if ~isempty(invis_trials)
                invis_catchProb=[invis_catchProb; squeeze(mean(invis_trials,1))];
            end
            
            ppantTrialDatawithDetails(ippant).TrialDetails(itrial).invisiblecatchlocations= invislocs;
        end
        %for save
        ppantTrialDatawithDetails(ippant).eachTrialCatchProb_combOnetoThreetargets = catchProb;
        ppantTrialDatawithDetails(ippant).mTrialCatchProb_combOnetoThreetargets = mean(catchProb,1);
        
        ppantTrialDatawithDetails(ippant).mTrial_invisibleCatchProb_combOnetoThreetargets = mean(invis_catchProb,1);
        
        
        %separate trials for significance calculations
        stackALL(ippant).data= catchProb;
        stackALL_m(ippant,:)= squeeze(mean(catchProb,1));
        stackALL_invis(ippant).data= invis_catchProb;
        stackALL_invis_m(ippant,:)= squeeze(mean(invis_catchProb,1));
        
        
    end
    catchBP_ProbtraceperPpant_onetothreetargets = stackALL;
    catchBP_ProbtraceperPpant_onetothreetargets_invis = stackALL_invis;
    
    catchBP_ProbtraceperPpant_onetothreetargets_mean = stackALL_m;
    catchBP_ProbtraceperPpant_onetothreetargets_invis_mean = stackALL_invis_m;
    
    
    
    save('Catchperformance', 'catchBP_ProbtraceperPpant_onetothreetargets', ...
        'catchBP_ProbtraceperPpant_onetothreetargets_invis',...
        'catchBP_ProbtraceperPpant_onetothreetargets_mean', ...
        'catchBP_ProbtraceperPpant_onetothreetargets_invis_mean')
    
    %also resave the ppantTrialdata since we stored the invisible catch
    %index.
    save('MD_AllBP_perppant.mat', 'ppantTrialDatawithDetails', '-append')
    
end

if job.sortCatchtracesbynumTargets==1
    cd(basefol)
    
    load('MD_AllBP_perppant.mat')
    
    for onsetoroffset = 1:2
        dataOUT=[];
        for ippant = 1:nppants
            ppantCatch1=[];
            ppantCatch2=[];
            ppantCatch3=[];
            ppantCatch4=[];
            
            for itrial=1:48
                
                Details= ppantTrialDatawithDetails(ippant).TrialDetails(itrial);
                
                Catch_any=[]; %reset variable per trial. how many to keep will change.
                
                numtargetsneeded=Details.TotalCatchTargetsRemoved;
                
                
                
                
                if onsetoroffset==1
                    Catch_any = Details.CatchBPs_longepoch; %10 second epoch of all required channels.
                else
                    Catch_any = Details.CatchBPsOFFSET_longepoch; %10 second epoch of all required channels.
                end
                
                %filter out those with/without button press at t=0.
                if onsetoroffset==1
                    onset=abs(window(1));
                    keeptrials=[];
                    for iBP=1:numtargetsneeded
                        %for each trial, only retain for analysis of prob if BP
                        %   not existing
                        tmpt=Catch_any(iBP,:);
                        %                         if sum(tmpt(1,onset-12:onset))<1 %if no buttons already pressed.
                        keeptrials=[keeptrials; tmpt];
                        %                         end
                        
                    end
                else
                    %                     keeptrials=Catch_any;
                    %at the moment no filtering.
                    onset=abs(window(2));
                    keeptrials=[];
                    for iBP=1:numtargetsneeded
                        %                         %for each trial, only retain for analysis of prob if BP
                        %   not existing
                        tmpt=Catch_any(iBP,:);
                        %                         if mean(mean(tmpt(:,onset-12:onset,1))) ==1 %12 frames is 0.2seconds
                        keeptrials=[keeptrials; tmpt];
                        %                         end
                        
                    end
                end
                %take mean over retained trials, having filtered.
                Catch_any=mean(keeptrials,1);
                
                switch numtargetsneeded
                    case 1
                        ppantCatch1=[ppantCatch1; Catch_any];
                    case 2
                        
                        ppantCatch2=[ppantCatch2; Catch_any];
                    case 3
                        ppantCatch3=[ppantCatch3; Catch_any];
                    case 4
                        ppantCatch4=[ppantCatch4; Catch_any];
                end
                
                
            end
            
            
            
            %save per ppant.
            
            dataOUT(ippant).ProbTraceforOnetargetcatch=ppantCatch1;
            dataOUT(ippant).mProbTraceforOnetargetcatch=mean(ppantCatch1,1);
            
            dataOUT(ippant).ProbTraceforTwotargetcatch=ppantCatch2;
            dataOUT(ippant).mProbTraceforTwotargetcatch=mean(ppantCatch2,1);
            
            
            dataOUT(ippant).ProbTraceforThreetargetcatch=ppantCatch3;
            dataOUT(ippant).mProbTraceforThreetargetcatch=mean(ppantCatch3,1);
            
            dataOUT(ippant).ProbTraceforFourtargetcatch=ppantCatch4;
            dataOUT(ippant).mProbTraceforFourtargetcatch=mean(ppantCatch4,1);
            
        end
        
        if onsetoroffset==1
            catchBP_ProbtraceperPpant_by_numtargetsabsent = dataOUT;
        else
            catchBP_ProbtraceperPpantOFFSET_by_numtargetsabsent = dataOUT;
        end
        
    end
    
    
    
    %%
    save('Catchperformance', 'catchBP_ProbtraceperPpant_by_numtargetsabsent',...
        'catchBP_ProbtraceperPpantOFFSET_by_numtargetsabsent',...
        '-append')
    
end

if job.calcNEWShuffledCatchOnsetBPProb==1
    %for each participant, calculate likelihood of BP around onset/offset of catch
    % in adjacent trials, to estimate RT.
    
    %%
    cd(basefol)
    
    load('MD_AllBP_perppant.mat')
    
    nshuff=2000;
    ntrials= 48;
    acrossall_mShuffledCatch_MD=zeros(length(ppantTrialDatawithDetails), 2, nshuff, 601);
    for ippant = 1:length(ppantTrialDatawithDetails)
        % calc shuffled trials SETs, each composed of a single
        % selection from any (1:48) of the trials. present trial included.
        
        pDetails =ppantTrialDatawithDetails(ippant);
        
        %List trial index by number of targets removed.
        numCatchppant=[];
        for itrial=1:ntrials
            numCatchppant(itrial,:)=pDetails.TrialDetails(itrial).TotalCatchTargetsRemoved;
        end
        
        %for each type of catch trial (loc targets removed). extract a random
        %BP for the same time window, same location(s), any trial.
        
        ppantShuffled=nan(nshuff,2,601); %[nshuffled, disap/reap, nsamps]
        
        for ishuff=1:nshuff
            
            shuffdata=[];
            for itrial = 1:ntrials
                %how many PMD targets were removed on this trial?
                catchnum = numCatchppant(itrial);
                %find locations to compare in random trials.
                locs = [pDetails.TrialDetails(itrial).TL_Catch, pDetails.TrialDetails(itrial).TR_Catch, pDetails.TrialDetails(itrial).BL_Catch, pDetails.TrialDetails(itrial).BR_Catch];
                locs=find(locs);
                
                
                %find trials to compare,
                availTrials = 1:48;
                
                %find real onset and offset of this catch
                catchtimes = [pDetails.TrialDetails(itrial).Catchstart_frames, pDetails.TrialDetails(itrial).Catchend_frames];
                
                % for each trial onset and offset, select the same time
                % window at random.
                for icatch=1:2 %for disap and reap
                    
                    tmp_time=catchtimes(icatch); %this is the catch start to collect about.
                    
                    
                    foundtrial=0;
                    iteration=1; % lets not get stuck!
                    checkedtrials=[];
                    
                    while foundtrial~=1
                        %select random trial
                        availtrials = setdiff(availTrials, checkedtrials);
                        trialis=availTrials(randi(length(availTrials)));
                        
                        %we'll break the while loop if checked all.
                        checkedtrials = [checkedtrials, trialis];
                        if length(checkedtrials)==48
                            break
                        end
                        
                        BPis = squeeze(pDetails.AllBPData(trialis,locs,:));
                        
                        %window about same catch time
                        if size(BPis,1)>size(BPis,2)
                            BPis=BPis';
                        end
                        
                        
                        if icatch==1 %catch onset
                            onset=abs(window(1));
                            storeBP = BPis(:, tmp_time+(window(1)):tmp_time+window(2));
                            
                            % check if BP is down!
                            for iBP=1:size(BPis,1)
                                tmpBP = storeBP(iBP,:);
                                % check if buttons down.
                                if sum(tmpBP(1,onset-60:onset))<60 %60 frames is 1 second.
                                    foundtrial=1;
                                    %we can continue with this data.
                                    BPout = tmpBP;
                                    break % breaks for loop
                                end
                            end
                            
                        else %catch offset, so move inspection window.
                            onset=window(2);
                            storeBP = BPis(:, tmp_time-(window(2)):tmp_time+abs(window(1)));
                        end
                    end
                    
                    %store BP of this window
                    
                    shuffData(itrial,icatch,:) = BPout;
                    
                end
                
            end
            
            %store mean across trials for this participant, this
            %shuffle.
            ppantShuffled(ishuff,:,:)=squeeze(nanmean(shuffData,1));
            
            
            
        end
        
        % take mean over trials.
        %will plot 99%CI of the shuffled data across trials.
        acrossall_mShuffledCatch_MD(ippant,:,:,:,:) = permute(ppantShuffled,[2 1 3]);
        disp(['fin shuff for ppant ' num2str(ippant)])
    end
    dimsare = [{'ppants'}, {'Disap/Reap'},{'nshufftrials'},{'nsamps'}];
    cd(basefol)
    
    acrossall_mShuffledCatch_MD_filtered=acrossall_mShuffledCatch_MD;
    %%
    save('ShuffledCatch', 'acrossall_mShuffledCatch_MD_filtered', 'dimsare', '-append')
end

if job.plotCatchtraces_ObsvsShuff_eachppant==1;
    cd(basefol)
    
    load('Catchperformance.mat')
    load('ShuffledCatch.mat')
    %collect onset trace across ppants for each type
    %%
    clf
    %collect RTs for later.
    onsetRTs=[];offsetRTs=[];
    for onsetoroffset=1%:2
        %reset outgoing variable
        
        figure(1);
        %         clf
        
        if onsetoroffset==1
            %             datais=catchBP_ProbtraceperPpant_by_numtargetsabsent;
            datais=catchBP_ProbtraceperPpant_onetothreetargets_mean;
            col = ['r'];
            xt=xtONSET;
            eventis= 'onset' ;
            
        else
            datais=catchBP_ProbtraceperPpantOFFSET_by_numtargetsabsent;
            col = ['r'];
            xt=xtOFFSET;
            eventis= 'offset' ;
            
        end
        %%
        %also need shuffled data.
        
        acrossall_mShuffledCatch = acrossall_mShuffledCatch_MD_filtered;
        %(filtered seems better due to the large amount of
        %disappearances, and hence BP at catch onset)
        
        for ippant=1:nppants
            
            RT=[];
            %% plot results of shuffled
            
            shuffdata= squeeze(nanmean(acrossall_mShuffledCatch(ippant,onsetoroffset,:,:),2));
            Shuffupper=[];
            Shufflower=[];
            
            
            % for each time point, calculate width of distribution:
            for ip = 1:length(shuffdata)
                
                shd= shuffdata(:,ip)';
                shd(shd==0)=.0001;
                %note not normally distributed!
                %try logit transform
                lp = log(shd)-log(1-shd);
                
                %now rearrange so that we can place correct scores.
                [lpn,id]=sort(lp);%convert to zscore
                %sort first.
                tz= zscore((lpn));
                
                %find closest to zscores to 99% CI (2.57), 95% CI= 1.96;
                
                AB=dsearchn(tz', [-1.96 1.96]');
                %this is upper and lower bound.
                CI= [lpn(AB(1)), lpn(AB(2))];
                
                %now convert back from logit.
                %CI upper bound (95%)
                Shuffupperbound(1,ip) = exp(CI(2)) ./ (exp(CI(2))+1);
                Shufflowerbound(1,ip) = exp(CI(1)) ./ (exp(CI(1))+1);
                
                
            end
            %%
            %plot median for shuffle
            shf=plot(xt, squeeze(nanmedian(shuffdata,1)), 'color', [.2 .2 .2]);
            hold on
            %add patch for shuffled data
            %Calculate the error bars
            
            uE=Shuffupperbound;
            lE=Shufflowerbound;
            
            %Make the patch
            yP=[lE,fliplr(uE)];
            xP=[xt,fliplr(xt)];
            
            %remove nans otherwise patch won't work
            xP(isnan(yP))=[];
            yP(isnan(yP))=[];
            
            colsh=[.2 .2 .2];
            H.patch=patch(xP,yP,1,'facecolor',colsh,...
                'edgecolor','none',...
                'facealpha',.15);
            
            
            %Make pretty edges around the patch.
            H.edge(1)=plot(xt,lE,'-','color',colsh);
            H.edge(2)=plot(xt,uE,'-','color',colsh);
            xlim([-1 4])
          
            %% plot  observed now
            
            hold on
            dt= datais(ippant,:);
            pl=plot(xt, dt, 'color', [col], 'linew', 3);
            %%
            
            %store real catch trace across avail trials for comparison of
            %significant departure from shuffled.
            
            p=[];
            
          
            if onsetoroffset ==1
                %find first point data lower is greater than shuff upper
                %             sigpoint=find((dt(120:end))>Shuffupper(120:end), 1 ); %finds first.
                sigpoints=find(dt>Shuffupperbound); %finds first.
                sigpoint = min(sigpoints(sigpoints>120));
                placement = (sigpoint)/60 -2;
                
                try onsetRTs(ippant)=placement;
                catch
                    onsetRTs(ippant)=nan;
                end
                placeY=0.95;
                colis= [0 .5 0];
                
            else %for offsets, we want the first ns point.
                %             sigpoint=find(dt(480:end)<Shuffupper(480:end), 1 ); %finds first.
                sigpoints=find(dt<Shuffupperbound);
                
                %first after offset
                sigpoint = min(sigpoints(sigpoints>480));
                
                placement = (sigpoint)/60 -8 ;
                try offsetRTs(ippant)=placement;
                    
                    if placement<.02;
                        
                        offsetRTs(ippant) =nan;
                    end
                catch
                    offsetRTs(ippant)=nan;
                end
                placeY=0.95;
                colis = 'r';
            end
            
            try plot([placement placement], [0 1], '--','color', 'b', 'linewidth', 1);
                placementT=1;
            catch
                placementT=0;
            end
            
            title(['Participant ' num2str(ippant)]);
            set(gca, 'fontsize', 15)
            ylim([ 0 1])
            
            %%
            hold on
            plot([ 0 0], ylim, ['k:'])
            ylabel('P(press)')
            shg
            
        end
        %%
        %     st=suptitle('Reaction Time to Catch events');
        %         st.FontSize=30;
        %         st.FontWeight='bold';
        
        set(gcf, 'color', 'w')
        cd(basefol)
        cd ../
        cd('Figures')
        cd('Catch analysis')
        %%
        
    end
    %%
    %% if plotting single ppant:
    xlabel(['Time from catch ' eventis ' [sec]'], 'fontsize', 15)
    
    ylabel('P(Button Press)', 'fontsize', 15)
    
    set(gca, 'fontsize', 1.5*fontsize)
    
    print('-dpng', 'Example single participant')
    
    %% also plot bar for Reaction times (reviewer request).
    figure(1);
    
%     goodppants=[1,2,4,6,9:16,18]; %

    clf
    plotbar=[nanmean(onsetRTs(goodppants)); nanmean(offsetRTs(goodppants))];
    plotst = [nanstd(onsetRTs(goodppants)); nanstd(offsetRTs(goodppants))];
    bh=bar(plotbar);
    bh.FaceColor= [0 .5 0];   
    %
    hold on
    bh = bar(2, plotbar(2), 'r');
    shg
    %
    hold on
    
    % adjust within subject error bars (as per Cousineau, 2005).
    tmpX = [onsetRTs(goodppants); offsetRTs(goodppants)];
    [h,p,tstat]= ttest(onsetRTs(goodppants), offsetRTs(goodppants));
    %within subj mean
    mX = nanmean(tmpX,1);
    %overallm
    mmX=nanmean(mX);
    
    %adjust.
    newX= tmpX - repmat(mX,[size(tmpX,1),1]) + repmat(mmX, [size(tmpX,1),size(tmpX,2)]);
    
    plotst = nanstd(newX,0,2)/sqrt(size(newX,2));
    %
    
    e=errorbar(plotbar, plotst, 'color', 'k');
    e.LineWidth=2;
    e.LineStyle='none';
    %
    set(gca, 'fontsize', 30')%, 'xtickmark', {'Onset' 'End'})
    set(gca, 'XTickLabel', {'Onset' 'Offset'})
    ylabel('reaction time [sec]')
    set(gcf, 'color', 'w')
    %     axis square
    shg
    %%
    print('-dpng', 'Catch RT across participants n=13');
end


if job.plotCatchtraces_ObsvsShuff==1 % overall group effects
    cd(basefol)
        
    load('Catchperformance.mat')
    load('ShuffledCatch.mat')
    
    cd(basefol)
    
    load('Catchperformance.mat')
    load('ShuffledCatch.mat')
    
    colis='r';
    
    %%% %%% for these plots, we can use CI or SEM for error bars
    %     useSEMorCI=2;
    
    % can also plot INVIS    
    plotINVIScatchonset=0;    
    %collect onset trace across ppants for each type
    %%
    
    
    RT=[];
    %% plot results of shuffled
    figure(1);
    clf
    for onsetoroffset = 1%:2
        
        if onsetoroffset==1
            %observed data first.
            if plotINVIScatchonset==0;
                % if not plotting the invisibe, plot the filtered (i.e. no buttons pressed).
                stack = catchBP_ProbtraceperPpant_onetothreetargets_mean;
            else
                stack = catchBP_ProbtraceperPpant_onetothreetargets_invis_mean;
            end
            %shuff to plot
            acrossall_mShuffledCatch= squeeze(acrossall_mShuffledCatch_MD_filtered(:,1,:,:)); %has an extra dimension for onsets/offsets.
            x = xtONSET;
            eventis= 'onset' ;
        else
            %not using
            %             stack=mean_catchBP_ProbtraceperppantOFFSET_bynum;
            %             acrossall_mShuffledCatch= acrossall_mShuffledCatchOffset_MD; %new versio of offset shuffled likelihood (separate job - see above).
            %
            %             x = xtOFFSET;
            %                  eventis= 'offset' ;
        end
        %%
        for idatatype=1:2 % shuff. vs obs
            switch idatatype
                case 1 %plot shuffled first
                    %mean across trials, within participants.
                    tmp= squeeze(acrossall_mShuffledCatch(goodppants,:,:));
                    %median across shuffs( 16 ppants).
                    tmp= squeeze(nanmedian(tmp,2));
                    shuffdata=tmp;
                    
                    mShuff = squeeze(nanmedian(tmp,1)); %m across ppants.
                    
                    colsh=[.1 .1 .1];
                    
                case 2
                    shuffdata = stack(goodppants,:);                    
                    mShuff = squeeze(nanmedian(shuffdata,1));
                    colsh='r';
            end
            
            % calculate CI
            Shuffupperbound=[];
            Shufflowerbound=[];
            stErrShuff=[];
            % for each time point, calculate width of distribution:
            for ip = 1:length(shuffdata)
                
                %%
                shd= shuffdata(:,ip)';
                shd(shd==0)=.0001;
                shd(shd==1)=.9999; % otherwise we get inf in computation.
                %note not normally distributed!
                %try logit transform
                lp = log(shd)-log(1-shd);
                
                %now rearrange so that we can place correct scores.
                [lpn,id]=sort(lp);%convert to zscore
                %sort first.
                tz= zscore((lpn));
                
                %find closest to zscores to 99% CI (2.57), 95% CI= 1.96;
                %             AB=dsearchn(tz, [-2.57 2.57]');
                AB=dsearchn(tz', [-1.96 1.96]');
                %this is upper and lower bound.
                CI= [lpn(AB(1)), lpn(AB(2))];
                
                %now convert back from logit.
                %CI upper bound (95%)
                
                Shuffupperbound(1,ip) = exp(CI(2)) ./ (exp(CI(2))+1);
                Shufflowerbound(1,ip) = exp(CI(1)) ./ (exp(CI(1))+1);
                
                
            end
            
            
            %%
            
            plot(x,mShuff, 'color', colsh, 'linew', 3);
            hold on
            %add patch
            %Calculate the error bars
            uE=Shuffupperbound;
            lE=Shufflowerbound;
            
            %Make the patch
            yP=[lE,fliplr(uE)];
            xP=[x,fliplr(x)];
            
            %remove nans otherwise patch won't work
            xP(isnan(yP))=[];
            yP(isnan(yP))=[];
            
            
            H.patch=patch(xP,yP,1,'facecolor',colsh,...
                'edgecolor','none',...
                'facealpha',.15);
            
            
            %Make pretty edges around the patch.
            H.edge(1)=plot(x,lE,'-','color',colsh);
            H.edge(2)=plot(x,uE,'-','color',colsh);
            xlim([-1 4])
            
            %%
            % now repeat the process for observed data.
            hold on
            
        end
        
        ylim([0 1])
        cr=plot([0 0 ], ylim, ['k' ':'], 'linewidth', 4);
        plot(xlim, [0 0 ], ['k', '-'], 'linewidth', 2);
        xlabel(['Time from PMD ' eventis ' [sec]'], 'fontsize', 15)
        
        ylabel('P(Button Press)', 'fontsize', 15)
        
        set(gca, 'fontsize', 1.5*fontsize)
        %     ylim([0 1])
    end
    set(gcf, 'color', 'w')
    %%
    cd(basefol)
    cd ../
    
    cd('Figures')
    cd('Catch analysis')
    %%
    set(gcf,'color', 'w')
    
    print('-dpng', ['Timecourse of each CatchDisap BP, across all with CI'])
end

%%
if job.calcMissedcatchesbyTargperppant==1
    cd(basefol)
    
    %%
    load('MD_AllBP_perppant.mat')
    pD=ppantTrialDatawithDetails;
    catchStruct=[{'ppant'}, {'trial'}, {'catchtype'}, {'failed'}, {'Locations'}, {'RejectTrial'}];
    %     req_frames=60;
    catchData_numbergone=[];
    catchData_location=[];
    
    %manually identified 'false' switches to ignore
    %ppant trial % These are when no fourth button could be pressed, on
    %account of three other buttons already being pressed.
    
    falseswitches=[];
    for ippant=1:nppants
        counter1=1;
        counter2=1;
        counter3=1;
        counter4=1;
        loc1=1;
        loc2=1;loc3=1;loc4=1;
        
        for itrial=1:48
            catchtype=pD(ippant).TrialDetails(itrial).TotalCatchTargetsRemoved;
            %
            
            
            chData=sum(pD(ippant).TrialDetails(itrial).CatchBPs_exact,2);
            
            totalabsence = pD(ippant).TrialDetails(itrial).totalCatchdur;
            
            % for both, adjust for RTs of 1 second.
            totalabsence = totalabsence-60;
            
            
            
            missedcatches=length(find(chData< totalabsence*.50));
            missedindx = find(chData<totalabsence*.5);
            %catch needed, find locations we missed.
            
            catchneeded = [pD(ippant).TrialDetails(itrial).TL_Catch, pD(ippant).TrialDetails(itrial).TR_Catch,pD(ippant).TrialDetails(itrial).BL_Catch,pD(ippant).TrialDetails(itrial).BR_Catch];
            %index of which buttonp presses were catch related
            bpindex=find(catchneeded>0);
            if find(chData<totalabsence*.5)>0 % if we had a miss.
                
                locationsmissed=bpindex(missedindx);
            else
                locationsmissed=0;
            end
            
            
            
            %record for later analysis in matrix.
            switch catchtype
                case 1
                    catchData_numbergone(ippant,counter1, catchtype) = missedcatches;
                    counter1=counter1+1;
                    
                    %Also include an index for whether to reject the trial or not.
                    
                    if missedcatches>0
                        rejectTrial=1; %failed single target catch.
                    else
                        rejectTrial=0;
                    end
                    
                case 2
                    catchData_numbergone(ippant,counter2, catchtype) = missedcatches;
                    counter2=counter2+1;
                    
                    
                    if missedcatches>1
                        rejectTrial=1; %failed dual target catch, but if they reported one, good enough.
                    else
                        rejectTrial=0;
                    end
                case 3
                    catchData_numbergone(ippant,counter3, catchtype) = missedcatches;
                    counter3=counter3+1;
                    
                    if missedcatches>2
                        rejectTrial=1;
                    else
                        rejectTrial=0;
                    end
                    
                case 4
                    catchData_numbergone(ippant,counter4, catchtype) = missedcatches;
                    counter4=counter4+1;
                    
                    if missedcatches>3
                        rejectTrial=1;
                    else
                        rejectTrial=0;
                    end
                    
                case 0
                    missedcatches=0;
                    locationsmissed=0;
                    rejectTrial=0;
            end
            
            
            % input data to structure.
            %override if this trial is a false switch (detected manually);
            %             if any(falseswitches(:,1)==ippant)
            %
            %                 tid=find(falseswitches(:,1)==ippant);
            %                 %check if this trial needs rejecting.
            %                 for ich=1:length(tid)
            %
            %                     if falseswitches(tid(ich),2)==itrial;
            %                         rejectTrial=nan;
            %                     end
            %                 end
            %             end
            catchStruct=[catchStruct; {num2str(ippant)}, {num2str(itrial)}, {num2str(catchtype)}, {num2str(missedcatches)}, {num2str(locationsmissed)}, {num2str(rejectTrial)}];
            
            
            
        end
        %%
    end
    %%
    req_frames= ['.5*totalabsence'];
    
    save('Catchperformance', 'catchStruct', 'catchData_numbergone', 'req_frames', '-append')
end
