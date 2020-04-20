%% File one: making the raw EEG file
%based on IG's Raw_EEG but with corrections for catch timing.

cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/EEG Will SSVEP Data 01')
basefol=pwd;
pdirs=dir([pwd filesep '*_*' 'EEG']);
%load channel data.
getelocs;
%%
for ippant = 1:length(pfiles)
    %%
    
    cd(basefol)
    % find ppant folder
    
    %change into that directory
    cd(pdirs(ippant).name);
    
% allcatchEEG = [];
    allEEG = dir([pwd filesep '*.vhdr']);
for ieeg = 1:length(allEEG)
    
    
    clear trialtrack catchinfo trialstart adjustedcatch EpochedEEG EpochEEG
    loadfile = allEEG(ieeg).name;
    
    EEG= pop_loadbv(pwd, num2str(loadfile));
    % replace channel data
    EEG.chanlocs = elocs(1:64);
 
    trialend = [];
    trialstart = [];
    for ievent = 1:length(EEG.event)
        eventype = EEG.event(ievent).type;
        if strcmp(eventype, 'S 20'); %catchmarker
            trialstart=[trialstart EEG.event(ievent).latency];
        end
        if strcmp(eventype, 'S 88') %trial end
            trialend =  [trialend EEG.event(ievent).latency];
        end
    end
    if length(trialstart)<length(trialend);
    trialstart = trialend - (60*EEG.srate); %60 second trials
    end
   % interpolate bad channels in data set.
 
   %reduce data size
   
   EEG.data = EEG.data(:,abs(trialstart(1)):trialend(end)+1);
      
% identify and interp bad channels using var methods:
%% modified from prepPipeline: Bigdely-Shamo et al., Frontiers 2015
tic
noisyOut=findNoisyChannels(EEG);
toc
% because we have SSVEPs, only remove those with deviations (not low
% channel correlation - as a localized responsed could be due to the entrainment).
%% 
badchans = [noisyOut.noisyChannels.badChannelsFromDeviation, noisyOut.noisyChannels.badChannelsFromDropOuts];

EEG= eeg_interp(EEG, badchans);

%sanity check plot
% clf;
% subplot(121); topoplot(std(EEG.data,0,2), elocs(1:64)); colorbar;
% subplot(122); topoplot(std(EEGn.data,0,2), elocs(1:64)); colorbar;
%%

    eegchans = EEG.data(1:64,:);
    
    Timing = (1:length(eegchans))/EEG.srate; %in seconds.
    %%
    %% Now epoch
    EpochedEEG=zeros(length(trialstart), 64, 60*EEG.srate+1); %60 seconds trials.
    
    %adjust for our shorter windows (since we reduced above, pre chan
    %rejection).
    trialend=trialend - abs(trialstart(1));
    trialstart=trialstart-abs(trialstart(1))+1;
    
    for itrial = 1:length(trialstart)
        
        %beware negatives
        if trialstart(itrial)<0
            %supplement with NaN in duration effected
            trial= EEG.data(:, 1:trialend(itrial)+1);
            nantrain = nan(1,abs(trialstart(itrial))+1);
            newd= zeros(64,60*EEG.srate+1);
            for ichan=1:64
                newd(ichan,1:length(nantrain)) = nantrain;
                newd(ichan,length(nantrain)+1:end) = squeeze(trial(ichan,:));
            end
            
            EpochedEEG(itrial,:,:) = newd;
        else
            trial= EEG.data(:, trialstart(itrial):trialend(itrial)+1);
            EpochedEEG(itrial,:,:) = trial(:, 1:size(EpochedEEG,3));
        end
    end
    
    
        if ieeg==1
            EpochedEEGdata =EpochedEEG;
        else
            
            %concat.
            EpochedEEGdata = cat(1, EpochedEEGdata, EpochedEEG);
            
        end
    
end
    
    
   cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA/EEG')
   try cd(pdirs(ippant).name)
   catch
       mkdir(pdirs(ippant).name)
       cd(pdirs(ippant).name)
   end
   savename =  ['P' num2str(ippant) 'RawEEG'];
    save(savename, 'EpochedEEGdata')
    
    disp(['fin ppant ' num2str(ippant) ])
    clearvars -except basefol pdirs elocs
    print('-dpng', 'interpolated channels')

end