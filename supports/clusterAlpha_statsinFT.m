  %% load alpha
  
      cd(basefol)
      cd('EEG')
      %%
%       load('GFX_LeftvsRight_TFdecomp.mat')
      load('GFX_LeftvsRight_TFdecomp_PMD.mat')
      tstamps = [1:1:size(all_leftON_hilb,3)]/250 - 3;

      Subchans = [12,13,15,16,17:32,47,48,49,50 51:64];

     %same var names
       trim = dsearchn(tstamps', [-2.5 2.5]');
     timeax = tstamps(trim(1):trim(2));
     
     %first remove any Nan pps.
     pnan = isnan(all_leftON_hilb(:,1,1));
     all_leftON_hilb(pnan,:,:) = [];
     pnan = isnan(all_rightON_hilb(:,1,1));
     all_rightON_hilb(pnan,:,:) = [];
     
   GFX_leftPFI = all_leftON_hilb(:,:,trim(1):trim(2));
   GFX_rightPFI = all_rightON_hilb(:,:,trim(1):trim(2));
   
   
     allPFI = cat(4, GFX_leftPFI, GFX_rightPFI);
     allPFI= squeeze(nanmean(allPFI,4));
     
     %perform baseline normalization:
     tmp = zeros(size(allPFI));
     %divide by baseline per subject.
     baselineis = dsearchn(timeax', [-2.5 -2]');
     
     for isub = 1:size(allPFI,1)
         for ichan = 1:64
             tmpchan = squeeze(allPFI(isub, ichan, :));
             tmpbase = nanmean(tmpchan(baselineis(1):baselineis(2)));
             tmp(isub,ichan,:)= (tmpchan./ tmpbase) -1;
         end
     end
     allPFI=tmp;
     
     %% from here, we need to conver to fieldtrip struture (not easy):
     %First to eeglab, then convert to fieldtrip.
     %get electrode locations.
     getelocs;
     
     allPFI_eeg = permute(allPFI, [2,3,1]);
     %%
     % import each sub to fieldtrip
     timelock_PFI=[];
     timelock_baseline=[];
     for isub=1:size(allPFI_eeg,3)
         tmpd = squeeze(allPFI_eeg(:,:,isub));
         
     %initialize EEGlab.
     EEG = pop_importdata('dataformat', 'array', 'nbchan', 64, 'data', 'tmpd', 'srate', 250);
     EEG.times = timeax;
     EEG.chanlocs = elocs(1:64);     
     %% convert to fieldtrip structure.
    data=   eeglab2fieldtrip(EEG, 'preprocessing');
    
  %need to repair time subfield
  for index = 1:EEG.trials
      data.trial{index}  = EEG.data(:,:,index);
      data.time{index}   = timeax;
  end
  % perform an average over the '1' trial.
  cfg=[];
  cfg.layout = 'acticap-64ch-standard2.mat';
  cfg.keeptrials = 'no';
  timelock_PFI{isub} = ft_timelockanalysis(cfg, data);
  
  %also compute baselines within subj, for comparison:
        data2=data;
        tmpD =data2.trial{1};
        bs_Data = zeros(size(tmpD));
        
        for ichan = 1:64
            bs = squeeze(tmpD(ichan,:));
            twin = dsearchn(data2.time{1}', [-2.5 -2]');
            subt = repmat(nanmean(bs(1, twin(1):twin(2))), [1, length(bs)]);
            bs_Data( ichan, :) = subt;
        end
        
        %store as a separate subjCellarray.
        % then feed into ft, for later analysis.
        data2.trial = {bs_Data};
        timelock_baseline{isub} = ft_timelockanalysis(cfg, data2);
  
  
     end
     
     %% now that we have two separate cell structures, per subj, per activity,
     %and per baseline. We can compute stats.
     
     %% need to compute a Grand average structure for fieldtrip functions.
  cfg=[];
  cfg.keepindividual = 'no'; 
  cfg.layout = 'acticap-64ch-standard2.mat';%   
  GrAvg_PFI = ft_timelockgrandaverage(cfg, timelock_PFI{:});
 GrAvg_bl = ft_timelockgrandaverage(cfg, timelock_baseline{:});
  %% sanity check,  plot ERP and topoplot.
  
  %first average across ppants.
  cfg=[];
  cfg.layout = 'acticap-64ch-standard2.mat';
  cfg.keeptrials = 'no';
  GFX_timelock = ft_timelockanalysis(cfg, GrAvg_PFI);
  %%
  cfg=[];
  cfg.layout = 'acticap-64ch-standard2.mat';
  cfg.xlim = [ -.5 1];
  cfg.zlim = [-.1 .1];
  cfg.colormap = 'parula';
  cfg.marker= 'off';
  cfg.style = 'fill';
  cfg.comment = ' ';
  cfg.colorbar= 'yes';
  clf;
  subplot(221); ft_topoplotER(cfg, GrAvg_PFI); shg
  %%
%also plot timecourse:
  cfg=[];
  cfg.xlim = [];
  cfg.channel = Subchans;
  cfg.linewidth = 2; 
  cfg.graphcolor = 'b';
  subplot(212); 
  ft_singleplotER(cfg, GrAvg_PFI); title('');
  %% stat evaluation of ERP:
  
% create a 'neighbourhood-structure', that defines how spatial clusters are
 cfg=[];  
 cfg.method = 'distance';
 cfg.layout= 'acticap-64ch-standard2.mat';
 cfg.neighbourdist = .15; % check!
cfg.channel = 1:64;
  neighbs = ft_prepare_neighbours(cfg);

  
  
%% create a FieldTrip style 'design'-matrix, conditions:

subj = length(timelock_PFI);
design = zeros(2,2*subj);
for i = 1:subj
design(1,i) = i;
end
for i = 1:subj
design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

%%

% compute the difference as a T-statistic, and do an inferential test
cfg                  = [];
cfg.channel          = {'all'};

cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 3;
cfg.neighbours = neighbs;  % same as defined for the between-trials experiment
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 1000;

%
cfg.design = design;
cfg.uvar                = 1; % are unit of observation is the subject (first row)
cfg.ivar                = 2; % %comaring the second.
cfg.spmversion = 'spm12'; % call most recent MEX files to avoid error.
%
staterp = ft_timelockstatistics(cfg, timelock_PFI{:},timelock_baseline{:});


%%

% get relevant (significant) values
neg_cluster_pvals = [staterp.negclusters(:).prob];
signif_clust = find(neg_cluster_pvals < staterp.cfg.alpha);

% find all topo cluster points associated with this cluster:
pos = ismember(staterp.negclusterslabelmat, 1); 
%% how long did the cluster last?
plength = zeros(1,length(pos));
for itime =1:size(pos,2)
    if any(pos(:,itime))
        plength(itime)=1;
    end
end
% how many electrodes max (for mask)? 
maxmask = any(pos, 2);

%% % plot stat output
%prep data for underneath topoplot cluster ID
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_FX = ft_math(cfg,  GrAvg_PFI,GrAvg_bl);
%% 
timestep = 0.25;      %(in seconds)
sampling_rate = 250;
sample_count  = length(staterp.time);
j = [-1:timestep:2];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot

m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in EEG samples

%
%% make sure channel orders match the data, and stat output

[i1,i2] = match_str(GA_FX.label, staterp.label);

%% now plot topo progression
clf
for k = 1:8
   subplot(2,4,k);
   cfg = [];
   cfg.xlim =[j(k) j(k+1)];
   
%    pos_int = zeros(numel(GA_FX.label),1);  
%    pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
%    
   cfg.highlight = 'off';
   cfg.maskparameter = maxmask;
   cfg.comment = 'xlim';
   cfg.commentpos = 'title';
   cfg.layout= 'acticap-64ch-standard2.mat';
   ft_topoplotER(cfg, GA_FX);
end
     
%% plot max effect.
figure(2); clf
twin = find(plength);
% topod = mean(GA_FX.avg(:, twin),2);
topod = mean(staterp.stat(:, twin),2);
adjustm = abs(maxmask-1);
subplot(121);
tp=topoplot(topod, elocs(1:64), 'conv', 'on');
subplot(122);
tp=topoplot(topod, elocs(1:64), 'conv', 'on', 'pmask',maxmask);
set(gcf, 'color', 'w')

     tmp= cbrewer('div', 'RdBu', 100);
     colormap(flipud(tmp));
%      caxis([-.1 .1]); 
caxis([ -3 3])
     c=colorbar;
     set(gca,'fontsize', 18)
     ylabel(c, '\itt- \rmvalue'); 
     shg
     %%
%now we have the SIG electrodes:
clustCHANS = adjustm;
save('Ft-clusterresultsPMD', 'clustCHANS', 'staterp',  'plength', 'timeax', 'maxmask')

%% 
% clf
% % sourcemodelling. 
% sourcemodel = ft_read_headshape('cortex_20484.surf.gii');
% 
% ft_plot_mesh(sourcemodel, 'facecolor', 'brain', 'edgecolor', 'none')
% camlight
% lighting gouraud
% %%
% cfg=[];
% cfg.method = 'mne';
% cfg.headmodel = sourcemodel;
% cfg.sourcemodel.pos = sourcemodel.pos;
% cfg.sourcemodel.tri = sourcemodel.tri;
% cfg.layout= 'acticap-64ch-standard2.mat';
% 
% % lay2= ft_prepare_layout(cfg)
% 
% 
% source = ft_sourceanalysis(cfg, GrAvg_PFI);
% 
% 
% 
% %%
% cfg = [];
% cfg.method         = 'surface';
% cfg.funparameter   = 'avg';
% cfg.maskparameter  = cfg.funparameter;
% cfg.funcolorlim    = [0.0 1.2];
% cfg.funcolormap    = 'jet';
% cfg.opacitylim     = [0.0 1.2];
% cfg.opacitymap     = 'rampup';
% cfg.projmethod     = 'nearest';
% cfg.surffile       = 'surface_white_both.mat'; % Cortical sheet from canonical MNI brain
% cfg.surfdownsample = 10;  % downsample to speed up processing
% ft_sourceplot(cfg, GrAvg_PFI);
% view ([90 0])
