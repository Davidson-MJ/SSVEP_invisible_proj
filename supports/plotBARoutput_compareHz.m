% plotBARoutput_compareHz 
%% after all HZ, plot the bar data.
    % dims are = 15 pfi, 15 catch, 20 pfi, 20 catch, 5pfi 5 catch)
    barax = [{'Target (f1), during PFI'}, {'Target (f1), during Catch'}, {'Background (f2), during PFI'}, {'Background (f2), during Catch'}, {'IM (f2-f1), during PFI'}, {'IM (f2-f1), during Catch'}]; 
    mDelta = squeeze(mean(collectBardata,2));
    stEDelta = std(collectBardata,0,2)./ sqrt(length(allppants));  
    figure(2);
    clf
    ylabel({['Disappearance - Reappearance'];['\Delta RESS log SNR']}); hold on
    bar(1:2, mDelta(1:2), 'b'); hold on; % plot different colors (TG1)
    bar(3:4, mDelta(3:4), 'Facecolor', [.2 .2 .2]); % plot different colors
    bar(5:6, mDelta(13:14), 'm');
 
    shg; hold on; eh=errorbar(1:6, mDelta([1,2,3,4,13,14]), stEDelta([1,2,3,4,13,14]), 'color', 'k', 'LineStyle', ['none'], 'LineWidth', 2);
    %
    ylim([-1 .5])
    %%
    set(gca, 'xtick', 1:6, 'xticklabel', barax, 'XtickLabelRotation', -45, 'ytick', [-1:.5:1], 'fontsize', 35)
    set(gcf, 'color', 'w')
    
    %%
    print('-dpng', ['Change in SNR Barchart summary from ' num2str(collectBarwindow(1)) ' to ' num2str(collectBarwindow(2)) 's'])
    %% alternatively, plot separately per TG or BG specific response (similar to Donner et al., JNS 2016).
    % rearrange all TG, BG and IM responses (data and errorbars)
    TG= mDelta([1:2; 5:6; 9:10]);
    eTG=stEDelta([1:2; 5:6; 9:10]);
    % avoid plotting 60 Hz
    BG= mDelta([3:4; 7:8]);
    eBG=stEDelta([3:4; 7:8]);
%     BG= mDelta([3:4; 7:8; 11:12]);
%     eBG=stEDelta([3:4; 7:8; 11:12]);
    IM= mDelta([13,14]);
    eIM= stEDelta([13,14]);
    % errorbar conditions to place on stacked plot 
    numgroups = 3;
    numbars = 2;   
    groupwidth = min(0.8, numbars/(numbars+1.5));
    
    allfontsize=20;
    yall=[-1 .7];
    
 %   %%%%%%%%%% PLOTTING
    % now plot this version for comparison
    figure(3); clf    
    
    % target specific responses
    
    subplot(1,3,1)
    bh=bar(TG, 'b'); shg; legend('PFI', 'PMD')        
    title('Target specific responses')
%     bh(2).FaceAlpha = .1;
bh(2).FaceColor='r';
    hold on
    
                for i = 1:numbars
                    hold on
                    %          Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
                    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
                    %
                    errb= errorbar(x, TG(:,i), eTG(:,i), 'k', 'linestyle', 'none');
                    errb.LineWidth = 1;
                    %
                end
   
     ylabel({['Disappearance - Reappearance'];['\Delta RESS log SNR']}); hold on
     xtspace=get(gca, 'xtick');   
     set(gca, 'xticklabel', {'Target (f1)', ' (2f1)', ' (3f1)'},  'ytick', [-1:.5:1], 'fontsize', allfontsize)
     ylim([yall])
     
     
%  Background specific responses
%
    subplot(1,3,2)
    bh=bar(BG, 'Facecolor', [.2 .2 .2]);    
%     bh(2).FaceAlpha = .1;
    bh(2).FaceColor='r';
    title('Background specific responses')
    legend('PFI', 'PMD')
    numgroups=2;
     hold on   
                for i = 1:numbars
                    hold on
                    %          Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
                    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
                    %
                    errb= errorbar(x, BG(:,i), eBG(:,i), 'k', 'linestyle', 'none');
                    errb.LineWidth = 1;
                    %
                end
     ylabel({['Disappearance - Reappearance'];['\Delta RESS log SNR']}); hold on
          set(gca, 'xticklabel', {'Background (f2)', ' (2f2)', 'Mix (60 Hz)'},  'ytick', [-1:.5:1], 'fontsize', allfontsize)
ylim(yall)

% IM finally

    subplot(1,3,3)
    % add some nans so we can stack this plot like the others
    IM(1:2,2)=NaN;
    eIM(1:2,2)=NaN;

    bh= bar(IM','m');    
%      bh(2).FaceAlpha = .1;
bh(2).FaceColor='r';
    hold on
    numgroups=2;
    
    for i = 1:numbars
        hold on
        %          Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
        %
        errb= errorbar(x, IM(i,:), eIM(i,:), 'k', 'linestyle', 'none');
        errb.LineWidth = 1;
        %
    end
%     errorbar(1:2, IM, eIM, 'k', 'linestyle','none')
    
    title('Intermodulation response')
    legend('PFI', 'PMD')    
     ylabel({['Disappearance - Reappearance'];['\Delta RESS log SNR']}); hold on
               set(gca, 'xticklabel', {'Intermodulation (f2-f1)'},  'ytick', [-1:.5:1], 'fontsize', allfontsize)

    ylim(yall)
    xlim([.5 1.5])
     shg
    set(gcf, 'color', 'w')
    %%  show ttests;
    tresults = [];
    for ihz = 1:14
        
        [~, p,~,stat]=ttest(collectBardata(ihz,:)); % compare to zero.
        tresults(ihz,1)= stat.tstat;
        tresults(ihz,2)= p;
    end
        %%
    ttype = repmat({'PFI';'Catch'}, 7,1);
    hztype= Interleave(peakfreqsare(1:7), peakfreqsare(1:7));
    pcorr=fdr(tresults(:,2));
    tresults(:,2)=pcorr;
    tresultsTable=table(hztype', ttype, tresults);    
    tresultsTable.Properties.VariableNames(1) = {'Hz'};
    disp(tresultsTable)