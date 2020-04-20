
%CALLED from s3_EEG_alphalateralization


job.performClusterstat_infieldtrip=0;

% once the above is done, progress:
%% >>>>>>
job.plotHILB_final =1; % paper figure, comparing PFI to PMD>
job.plotAlpha_f1corr =0;
job.plotNormf1 = 0; % added this plot to help interpretation of correlations.

%old plots:
job.plotTOPOdifferences = 0;
job.plotHILBenvelopes =0;

usePFIorPMD=1; % for the above.



%% run field trip cluster correction, to identify significant electrodes.

if job.performClusterstat_infieldtrip==1
    try load('Ft-clusterresults.mat')
    catch
        clusterAlpha_statsinFT;
    end
end
 

 %topoplot first, 
 %%
 if job.plotTOPOdifferences == 1
     
     cd(basefol)
     cd('EEG')
     if usePFIorPMD==1
         load('GFX_LeftvsRight_TFdecomp.mat')
     else
         load('GFX_LeftvsRight_TFdecomp_PMD.mat')
     end
     
     
     
     
     twinsare =[-2.5, 0; 0 2.5];
     for twin= 2%:2 % two time windows, before and after BP
         
         
         topowindow = dsearchn(tstamps',[twinsare(twin,:)]');
         
         
         %plot topo of difference:
         figure(1);  set(gcf, 'color', 'w', 'units', 'normalized', 'position',[.5 .5 .5 .5]);
         
         %% per participant, plot the alpha power
         for itype = 1:4
             switch itype
                 case 1
                     dis = all_leftON_hilb;
                     tis= 'PFI disap left';
                 case 2
                     dis = all_rightON_hilb;
                     tis= 'PFI disap right';
                 case 3
                     dis = all_leftOFF_hilb;
                     tis= 'PFI reap left';
                 case 4
                     dis = all_rightOFF_hilb;
                     tis= 'PFI reap rift';
             end
                 figure(itype);set(gcf, 'color', 'w', 'units', 'normalized', 'position',[.5 .5 .5 .5]);
             
             for ippant = 1:size(all_leftOFF_hilb,1)
             
                 subplot(4,4, ippant)
             
                 plotd = squeeze(mean(dis(ippant,:, topowindow(1):topowindow(2)),3));
                 topoplot(plotd, elocs(1:64)); c= colorbar;
                 title(num2str(ippant));
             
             end
             %
%              figure(twin);
%              subplot(2,2,itype);
%              set(gcf, 'color', 'w', 'units', 'normalized', 'position',[.5 .5 .5 .5]);
%              %GFx.
             plotdn = squeeze(mean(mean(dis(:,:,topowindow(1):topowindow(2)),3),1));
%              topoplot(plotdn, elocs(1:64),'conv', 'on'); c= colorbar;
%              caxis([ 0 40])
%              title({['GFX ' tis]; [num2str(twinsare(twin,:))]});
             
             TOPOfinal(itype,:) = plotdn;
         end
         %% now differences
         figure(10+twin);
         % Grand aveage
         % subplot(3,2,1:2); topoplot(mean(TOPOfinal,1), elocs(1:64),'conv', 'on'); c= colorbar; title('Grand Avg');
         %PFI left vs right
         dp = TOPOfinal(1,:) - TOPOfinal(2,:);
         subplot(2,2,1); topoplot(dp, elocs(1:64),'conv', 'on'); c= colorbar; title('PFI L-Targs - PFI R-Targs');
         caxis([-1 2])
         dp = mean(TOPOfinal(1:2,:),1) - mean(TOPOfinal(3:4,:),1); %
         subplot(2,2,2); topoplot(dp, elocs(1:64),'conv', 'on'); c= colorbar; title('PFI vs reapp');
         caxis([-1 2])
         dp = TOPOfinal(1,:) - TOPOfinal(3,:); %
         subplot(2,2,3); topoplot(dp, elocs(1:64),'conv', 'on'); c= colorbar; title('PFI L - L reapp');
         caxis([-1 2])
         dp = TOPOfinal(2,:) - TOPOfinal(4,:); %
         subplot(2,2,4); topoplot(dp, elocs(1:64),'conv', 'on'); c= colorbar; title('PFI R - R reapp');
         caxis([-1 2])
         %%
         sgtitle(['Effects in window ' num2str(twinsare(twin,:))])
     end
 end
 %%
 
 % LARGE posterior electrode subset.

 Subchans = [12,13,15,16,17:32,47,48,49,50 51:64];

 Ochans = [28:32, 60:64];
 % and define electrode locations for left and right side.
Lchans = [17,18,19,23,24,28,29,51,52,56,57,60,61];
Rchans = [20,21,22,26,27,31,32,54,55,58,59,63,64];
    %%
 if job.plotHILBenvelopes ==1
     
     %% % %% plot trace?
     cd(basefol)
     cd('EEG')
     
     
     if usePFIorPMD==1
         load('GFX_LeftvsRight_TFdecomp.mat')
     else
         load('GFX_LeftvsRight_TFdecomp_PMD.mat')
     end
     %first, trim all, 
      tstamps = [1:1:size(all_leftON_hilb,3)]/250 - 3;

     trim = dsearchn(tstamps', [-2.5 2.5]');
     timeax = tstamps(trim(1):trim(2));
   GFX_leftPFI = all_leftON_hilb(:,:,trim(1):trim(2));
   GFX_rightPFI = all_rightON_hilb(:,:,trim(1):trim(2));
   GFX_leftReap = all_leftOFF_hilb(:,:,trim(1):trim(2));
   GFX_rightReap = all_rightOFF_hilb(:,:,trim(1):trim(2));
   
   
     figure(44); clf
     set(gcf, 'units', 'normalized', 'position', [0 0 1 1 ]);
     clf
     %data should be chansxframesxn (for comparison.
     allPFI = cat(4, GFX_leftPFI, GFX_rightPFI);
     allPFI= squeeze(mean(allPFI,4));
     allReap = cat(4, GFX_leftReap, GFX_rightReap);
     allReap= squeeze(mean(allReap,4));
     
     %now combine for eeg plot:
     pEEG_tmp = cat(4, allPFI, allReap);
     pEEG= squeeze(mean(pEEG_tmp,1));
     %every channel
     
%      plottopo(pEEG, 'chanlocs', elocs(1:64), 'limits', [200, 1300, 10 40], 'ydir', 1)
     
% now plot HILB envelope of this activity.
     %subset.
%      Subchans = [12:32, 47:64];

    

     %plot all DISAP vs REAP. (no lateralization consideration)
     plot_PFIvsReap_Alpha;
    
     % plot lateralized results. Left vs Right Targ disappearance.
%      plot_PFI_LvsR_Alpha;
 end
%%
if  job.plotHILB_final ==1; % paper figure, comparing PFI to PMD>
    
    %laod both PFI and PMD>
     %first, trim all, 
     %%
     cd(basefol)
     cd('EEG')
     for iPFIvsPMD=1:2
     
         if iPFIvsPMD==1
             load('GFX_LeftvsRight_TFdecomp.mat')
             load('Ft-clusterresults.mat', 'clustCHANS');
             %load channels for plots.
             clustCHANS_PFI = clustCHANS;
         else
             load('GFX_LeftvsRight_TFdecomp_PMD.mat')
             load('Ft-clusterresultsPMD.mat', 'clustCHANS');
             clustCHANS_PMD = clustCHANS;
         end
      tstamps = [1:1:size(all_leftON_hilb,3)]/250 - 3;

     %same var names
     trim = dsearchn(tstamps', [-2.5 2.5]');
     timeax = tstamps(trim(1):trim(2));
   GFX_leftPFI = all_leftON_hilb(:,:,trim(1):trim(2));
   GFX_rightPFI = all_rightON_hilb(:,:,trim(1):trim(2));
   GFX_leftReap = all_leftOFF_hilb(:,:,trim(1):trim(2));
   GFX_rightReap = all_rightOFF_hilb(:,:,trim(1):trim(2));
   
   
     %data should be chansxframesxn (for comparison.
     if iPFIvsPMD==1
         
     allPFI = cat(4, GFX_leftPFI, GFX_rightPFI);
     allPFI= squeeze(mean(allPFI,4));
     allReap = cat(4, GFX_leftReap, GFX_rightReap);
     allReap= squeeze(mean(allReap,4));
     else
         allPFI_PMD = cat(4, GFX_leftPFI, GFX_rightPFI);
        allPFI_PMD= squeeze(mean(allPFI_PMD,4));
        allReap_PMD = cat(4, GFX_leftReap, GFX_rightReap);
     allReap_PMD= squeeze(mean(allReap_PMD,4));
     end
     end
     %%
     % now that we have the data, plot both on the same axes:
     clf
     plot_PFIvsPMD_Alpha;
         
         
    
end
%%
if job.plotAlpha_f1corr ==1;

    %% for both, lets compare the same window. Cluster start end in seconds, 
    % is :
    
    %load the RESS f1 data
    
    clf
    for iPFIvsPMD=1%1:2
        cd(basefol)
        cd('EEG')
        cd('GFX_EEG-RESS')        
        %%
        if iPFIvsPMD==1
        load(['GFX_PFIperformance_withSNR_15_min0_RESS'])
         
        cd(basefol)
        cd('EEG')        
            load('Ft-clusterresults.mat')
            load('GFX_LeftvsRight_TFdecomp.mat')
            coln='b';
        else
            
            %load PMD f1 and convert to same size.
            load('GFX_Catch_performance_withSNR_15_min0_RESS.mat')
            %third dimension is the BP locked catch onset data.
               
        cd(basefol)
        cd('EEG')   
            storeacrossPpant_onsetSNR = squeeze(storeacrossPpant_catchEVENTS_SNR(:,3,:,:));
            load('Ft-clusterresultsPMD.mat')            
            load('GFX_LeftvsRight_TFdecomp_PMD.mat')
            coln='r';
        end
        
        twin= [];
        twin(1) = staterp.time(find(plength,1, 'first'));
        twin(2) = staterp.time(find(plength,1, 'last'));
        
%         twin = [-1  1.5]
        %extract the SNR over the cluster window:
        avgg= dsearchn(timeidDYN', twin');
        RESSsnr = mean(squeeze(nanmean(storeacrossPpant_onsetSNR(:, :, avgg(1):avgg(2)),3)),2);
        RESSsnr_t = squeeze(nanmean(storeacrossPpant_onsetSNR,2));
        figure(3);
        hold on
        plot(squeeze(mean(RESSsnr_t)))
        %now do the same for alpha, per subject.
        %%
        % load alpha
        cd(basefol)
        cd('EEG')
        
        
        tstamps = [1:1:size(all_leftON_hilb,3)]/250 - 3;
        
        %same var names for alpha in both cases.
        trim = dsearchn(tstamps', [-2.5 2.5]');
        timeax = tstamps(trim(1):trim(2));
        GFX_leftPFI = all_leftON_hilb(:,:,trim(1):trim(2));
        GFX_rightPFI = all_rightON_hilb(:,:,trim(1):trim(2));
        
        
        allPFI = cat(4, GFX_leftPFI, GFX_rightPFI);
        allPFI= squeeze(mean(allPFI,4));
        
        %convert to log scale to avoid non-normal distribution.
        
        
        
        %perform baseline normalization:
        tmp = zeros(size(allPFI));
        %divide by baseline per subject.
        baselineis = dsearchn(timeax', [-2.5 -2]');
        
        for isub = 1:size(allPFI,1)
            for ichan = 1:64
                tmpchan = squeeze(allPFI(isub, ichan, :));
                tmpbase = nanmean(tmpchan(baselineis(1):baselineis(2)));
                tmp(isub,ichan,:)= (tmpchan./ tmpbase) ;
            end
        end
        allPFI=tmp; % this is the alpha data.
        
        
        
        clustCHANSn = abs(clustCHANS-1);
        clustwindow = find(plength);

        Alphacorr= mean(squeeze(nanmean(allPFI(:,find(clustCHANSn),clustwindow),2)),2);
        
        figure(1); 
%         subplot(1,2, iPFIvsPMD);  
        hold on;
        set(gcf, 'units', 'normalized', 'position', [0 1 .29 .29], 'color', 'w')
        
        %% before scatter, check for outliers:
        rma = isoutlier(Alphacorr);
        rmp = isoutlier(RESSsnr);
        rmall = [find(rma), find(rmp)];
        
        outlierplots = [RESSsnr(rmall), Alphacorr(rmall)];
        if ~isempty(rmall)
            Alphacorr(rmall) = [];
            RESSsnr(rmall) = [];
        end
        
%         %log transform?
%         RESSsnr = log(RESSsnr);
%         Alphacorr= log(Alphacorr);
        %% plot remaining data points.
        H=scatter( (RESSsnr),(Alphacorr), 100, coln); shg
        H.Marker = '*';
        H.LineWidth = 3; 
        hold on 
        %now plot outliers.
         H=scatter(outlierplots(:,1), outlierplots(:,2), 100, coln); shg
        H.Marker = 'o';
        H.LineWidth = 3; 
        %%
        
        [r,p] = corrcoef([Alphacorr, RESSsnr]);
%         [rho, pval] = corr(Alphacorr, RESSsnr, 'type', 'Spearman');
%               title(['rho = ' sprintf('%.2f', rho) ', p=' sprintf('%.2f', pval)])
        P= polyfit(RESSsnr, Alphacorr,1);
        
        Y_fit= polyval(P, RESSsnr);
        
        hold on
        lg= plot(RESSsnr, Y_fit, 'linew', 4, 'color', coln);
        xlabel({['f1 RESS log(SNR)']})
        ylabel({['Evoked Alpha during PFI']});
        set(gca, 'fontsize', 20)

        axis square
if iPFIvsPMD ==1         
    ylim([ 0.75 1.05])
else
    ylim([ 0.55 1.25])
end
        
        legh=legend(lg, ['r = ' sprintf('%.2f', r(1,2)) ', p=' sprintf('%.3f', p(1,2))]);
legh=legend(lg, ['\itr\rm = .58, \itp\rm = .019']);
        legh.FontSize = 20;
        legh.Location = 'SouthWest';
        shg
    end
end
     

if job.plotNormf1 == 1; % added this plot to help interpretation of correlations.
 clf
    for iPFIvsPMD=1:2
        cd(basefol)
        cd('EEG')
        cd('GFX_EEG-RESS')        
        %
        if iPFIvsPMD==1
        load(['GFX_PFIperformance_withSNR_15_min0_RESS'])
         
        cd(basefol)
        cd('EEG')        
            load('Ft-clusterresults.mat')
            load('GFX_LeftvsRight_TFdecomp.mat')
            coln='b';
        else
            
            %load PMD f1 and convert to same size.
            load('GFX_Catch_performance_withSNR_15_min0_RESS.mat')
            %third dimension is the BP locked catch onset data.
               
        cd(basefol)
        cd('EEG')   
            storeacrossPpant_onsetSNR = squeeze(storeacrossPpant_catchEVENTS_SNR(:,3,:,:));
            load('Ft-clusterresultsPMD.mat')            
            load('GFX_LeftvsRight_TFdecomp_PMD.mat')
            coln='r';
        end
        
        twin= [];
        twin(1) = staterp.time(find(plength,1, 'first'));
        twin(2) = staterp.time(find(plength,1, 'last'));
        
        %time series of RESS:
        RESSsnr_t = squeeze(nanmean(storeacrossPpant_onsetSNR,2));
     
        
        %now normalize, as we do for alpha:
        
        %perform baseline normalization:
        tmp = zeros(size(RESSsnr_t));
        %divide by baseline per subject.
        baselineis = dsearchn(timeidDYN', [-2.5 -2]');
        
        for isub = 1:size(RESSsnr_t,1)
            
                tmpchan = squeeze(RESSsnr_t(isub, :));
                tmpbase = nanmean(tmpchan(baselineis(1):baselineis(2)));
                tmp(isub,:)= (tmpchan./ tmpbase) ;
            
        end
        all_f1_norm=tmp; % this is the alpha data.
        %%   
        figure(3);  hold on
        plot(timeidDYN, squeeze(nanmean(all_f1_norm,1)), 'color', coln);
        
        
    end
    %%
end