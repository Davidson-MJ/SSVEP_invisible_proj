 

close all
normON=1;
logscale =0;
plotALPHAdistribat_tZERO=0; % sanity check the data is normally distributed before ttests.
legH=[];


ttestdata=[];
colsare = {'b', 'r'}; % preset colours.
for iPFI = 1:2
    switch iPFI
        case 1
            dtype = allPFI;
            usechans = abs(clustCHANS_PFI-1)';
        case 2
            dtype = allPFI_PMD;
            
    end
    colis= colsare{iPFI}; % set colour
    %%
    if logscale==1
        %log transform.
        dtype= log(dtype);
        ntype=2; % update baseline
        compareTO = 0; % baseline
        ylabis= 'Evoked log(Alpha) amplitude';
        plotsigat = -.2;
        ylimsare=[-.16 0.15];
    else
        compareTO=1; % for ttests
        ntype=1; % perform division by baseline.
        plotsigat = .86; % y value for asterisks.
        ylabis= 'Evoked Alpha amplitude';
        ylimsare=[.83 1.1];
    end
    %%
    if plotALPHAdistribat_tZERO==1
        %check distribution -> log transform?
        timezero = dsearchn(timeax', 0);
        figure(10); hold on;
        % average across channels:
        plotd = squeeze(nanmean(dtype(:,find(usechans),:),2));
        %plot distributiona time zero
        hh=histogram(squeeze(plotd(:,timezero)), size(plotd,1));
        hh.FaceColor = colis;
        
    end
      
         %% first plot the alpha power, baseline window (across channels),
         %before normalization.
         %show topoplot of certain periods, and marked electrodes.
         if iPFI==1
         figure(1); clf
         subplot(1,2,1);
         %which time windows for the topoplots?
         t1 = dsearchn([timeax'],[-2.5 -1.5]');
          
         %take average of all PFI and PMD.
         alld = cat(4, allPFI, allPFI_PMD);
         dtmp = squeeze(nanmean(alld,4));
         
         tp = squeeze(nanmean(nanmean(dtmp(:,:,t1(1):t1(2)),3),1));
         topoplot(tp, elocs(1:64), 'conv', 'on');         
         c=colorbar; 
         ylabel(c, 'Alpha amp. (\rm\muV)', 'fontsize', 15)
         caxis([0 15])
         set(gca, 'fontsize', 15)    
         colormap('viridis')
         end
         
         %% for time-series, normalize within ppant.
         if normON==1 % normalize within ppant.
             tmp = zeros(size(dtype));
             %divide by baseline per subject.
             baselineis = dsearchn(timeax', [-2.5 -2]');
             
             for isub = 1:size(dtype,1)
                 for ichan = 1:64
                     tmpchan = squeeze(dtype(isub, ichan, :));
                     tmpbase = nanmean(tmpchan(baselineis(1):baselineis(2)));
                     if ntype==1
                         tmp(isub,ichan,:)= (tmpchan./ tmpbase) ;
                     elseif ntype==2
                         tmp(isub,ichan,:)= (tmpchan - tmpbase) ; % if log scale.
                     end
                 end
             end
             dtype=tmp;
         end
          
         % plot time-series averaged across channel subset
         plotd = squeeze(nanmean(dtype(:,find(usechans),:),2));       
         plotMean = squeeze(nanmean(plotd,1));
         stE = CousineauSEM(plotd);

         figure(2); hold on
         if iPFI==2
             plotMean = smooth(plotMean, 60);
         end
         
         sh=shadedErrorBar(timeax, plotMean, stE, colis,1);
         sh.mainLine.LineWidth=3;
         legH(iPFI) = sh.mainLine;
         set(gca, 'fontsize', 20)
         hold on
         yvec=[0.8 1.15 1.15 0.8];
         if logscale==1 % shift to zero.
             yvec =yvec-1;
         end
         
        %also add patche for baseline window.
        pch = patch([-2.5 -2.5 -1.5 -1.5],yvec , [.9 .9 .9]);
        pch.FaceAlpha = .3;
        pch.LineStyle = 'none';
        
        ttestdata(iPFI,:,:) = plotd; % save data for ttests at each time point.
     end
 
     % perform sig tests.
  
     for iTEST =1:2 % compare PFI and then PMD, to baseline.
         pvals=[];
         for itime = 1:size(plotd,2)
             [~, pvals(itime)] = ttest(ttestdata(iTEST,:,itime),compareTO); % PFI
             if pvals(itime)<.05
                 % plot sig uncorr
                 text(timeax(itime), plotsigat-(iTEST/100), '*','fontsize', 25, 'color', colsare{iTEST})
             end
         end
         % perform temporal cluster correction for pval chains:
         [pv, pvals_cor ] = fdr(pvals, .05);
         
         if pv~=0 && pv<.05
             textis = ['\itp \rm= ' sprintf('%.3f', pv1)]; % sig output.
         else
             textis = '\itns'; % ns
         end
         %add to figure
         text(timeax(find(pvals<.05, 1, 'last')+50), plotsigat-(iTEST/100),textis ,'fontsize', 15, 'color', colsare{iTEST})
     end
     
     % finish plot labelling     
     ylim(ylimsare)
     xlim([-2.5 2.5])     
     xlabel('Time from button press');
     ylabel(ylabis)     
     hold on;
     plot([ 0 0], ylim, ['k-'])
     plot(xlim, [1 1], ['k-'])
     %legend:
     lg= legend([legH(1) legH(2)], 'PFI', 'PMD');
     lg.Position = [0.55 0.78 .16 .14];
     %
     set(gcf, 'color', 'w')
     