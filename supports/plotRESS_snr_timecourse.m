
% plotRESS_snr_timecourse

getelocs
    
    clf
    rmvbase=0; % normalize
    checksigON=1; % check sigs between directions
    checkcluster=1; % clusters only
    
    collectBarwindow = [0, 1]; % seconds to collect SNR difference between disappearance/reappearance.
    
    % we have 6 combinations for Disap - Reap : ( 5/15/20 hz x PFI/Catch)
    
    collectBardata = zeros(length(peakfreqsare)*2, length(allppants));
    
    
    % plot specs:
    plcount=1;
    legendprint=[];
    cd(basefol)
    cd('EEG')
    cd('GFX_EEG-RESS')
    colsare={'b' , 'k','b','k','b','g','m'}; % blue for tg, black for BG.
    hzare={'TG (f1)' , 'BG (f2)','TG (2f1)','BG (2f2)','TG (3f1)','TGBG','IM (f2-f1)'}; % blue for tg, black for BG.
    hzare={'Target (f1)' , 'Backgroud (f2)','TG (2f1)','BG (2f2)','TG (3f1)','TGBG','IM (f2-f1)'}; % blue for tg, black for BG.
    linetypes={'--', '-', [],[],'--'};
    legendis=[];
    lgc=1;
    
    figure(1); hold on;
    clearvars yl;
    
    lsh=[];
    % per frequency, plot both the PFI data AND then CATCH data, and produce
    % barchart output.
    plotBARoutput=0; % plot the output or no, not need all hzis = 1:7 for BAR output.
    barcounter=1;
    %%
    for hzis=1:2%:7% need all 7 for bar data.
        clf
        for idtype=1%1:2
            switch idtype
                case 1
                    useD= 'PFI';
                    
                    falpha = .25;  % alpha value to darken plot.
                    %plot colour:
                    col=colsare{hzis};
                    
                case 2
                    useD= 'Catch_';
                    sigcolor='r';  % color of sig points on figure,
                    falpha = .15;  % alpha value to lighten plot.
                    
                    %** changed catch colour
                    %all now  red
                    col = 'r';
                    
            end
            cd(basefol)
            cd('EEG')
            cd('GFX_EEG-RESS')
            
            usehz=peakfreqsare(hzis);
            
            load(['GFX_' useD 'performance_withSNR_' num2str(usehz) '_min0_RESS'])
            
            sigcolor=col;
            
            tbase = timeidDYN;
            
            hzp= hzare{hzis};
            
            ttestdata=[];
            %         legendprint=[];
            
            diffTEMP = []; % we want the difference between disap/reap for each iteration.
            for itimezero=1:2%1%:2%1:2
                
                switch itimezero
                    case 1
                        % different design for PFI vs Catch
                        if idtype==1
                            useSNR=storeacrossPpant_onsetSNR;
                        else
                            % BP catch onset =
                            useSNR = squeeze(storeacrossPpant_catchEVENTS_SNR(:,3,:,:));
                        end
                        
                        chtype='button press/release';
                        linet='--';
                        
                        % new marker types:
                        if idtype==1
                            markert ='o';
                        else
                            markert='square';
                        end
                    case 2
                        
                        if idtype==1
                            useSNR=storeacrossPpant_offsetSNR;
                        else
                            %BPcatch offset = dim 4
                            useSNR = squeeze(storeacrossPpant_catchEVENTS_SNR(:,4,:,:));
                        end
                        chtype='button press/release';
                        
                        linet='-';
                        markert='none';
                        
                    case 3
                        useSNR=storeacrossPpant_onsetSNR+storeacrossPpant_offsetSNR;
                        chtype = 'subjective report';
                        linet=':';
                        
                end
                
                figure(1);
                
                used=useSNR;
                
                hold on
                
                
                %plot across ppant trace
                ppantMeanSNR= squeeze(mean(used,2));
                %             subplot(2,1,itimezero);
                %             plot(ppantMeanSNR');
                
                
                %adjust standard error as per COusineau(2005)
                %confidence interval for within subj designs.
                % y = x - mXsub + mXGroup,
                
                x = ppantMeanSNR;
                
                mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
                mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
                
                %for each observation, subjtract the subj average, add
                %the group average.
                NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
                
                
                stE = std(NEWdata)/sqrt(length(allppants));
                
                
                if rmvbase==1
                    ppantMeanSNR=ppantMeanSNR-mean(ppantMeanSNR(:));
                end
                %plot
                sh=shadedErrorBar(timeidDYN, mean(ppantMeanSNR,1),stE,[],[1]);
                
                sh.mainLine.LineWidth=3;
                
                sh.mainLine.Color = col;
                sh.mainLine.LineStyle=linet;
                sh.patch.FaceColor = col;
                sh.patch.FaceAlpha= falpha;
                sh.edge(1).Color = col;
                sh.edge(2).Color = col;
                %
                sh.mainLine.Marker = markert;
                sh.mainLine.MarkerSize=25;
                
                lsh(itimezero) = sh.mainLine;
                
                xlabel(['Time from ' chtype])%
                
                
                if itimezero==3
                    
                    ylabel({['\Delta RESS log(SNR)']})
                    
                else
                    
                    if rmvbase~=1
                        ylabel({['RESS log(SNR), ' hzp]})
                        ylabel({['RESS log(SNR)']})
                    else
                        ylabel({['mean subtracted'];['RESS log(SNR)']})
                    end
                    
                end
                
                
                %store output for legend/stats
                ttestdata(itimezero,:,:) = ppantMeanSNR;
                legendprint(lgc)=sh.mainLine;
                %legend based on type
                legendis = [legendis {chtype}];
                
                
                % also collect bar data for plottting
                % over which temporal window?
                twind= dsearchn(timeidDYN', collectBarwindow');
                
                
                
                %average within ppant:
                storedata = squeeze(mean(ppantMeanSNR(:, twind), 2));
                
                diffTEMP(itimezero,:) = storedata; % we want the difference between disap/reap for each iteration.
                
                if itimezero==2; % then perform subtraction and store:
                    % we have 6 combinations (PFIx Catch, 5,15,20 hz)
                    collectBardata(barcounter, :) =diffTEMP(1,:) - diffTEMP(2,:);
                    barcounter=barcounter+1;
                end
                
                
                plcount=plcount+1;
                
                
                lgc=lgc+1;
                
                
            end % by type
            
            axis tight
            if rmvbase==1
                ylim([-.1 .15])
                hold on;
                plot(xlim, [0 0], ['k:'])
            end
            if itimezero==3
                plot(xlim, [0 0], ['k:'])
            end
            
            
            % %%%%%%%% END OF PLOTTING, SIG TESTS BELOW
            if checksigON==1
                temporalclustercheck
            end
            
            
            
        end % PFI vs Catch
        
        hold on
        
        
        %%
        %% adjust ylims.
        
        
        set(gca, 'fontsize', 25)
        
        set(gcf, 'color', 'w')
        
        yl=get(gca, 'ylim');
        
        %Percent of ydim.
        PTy=(yl(2)-yl(1))*.35;
        ylim([yl(1) yl(2)+PTy]) % f1
        
        plot([0 0 ], ylim, ['k-']);
        lg=legend([lsh(1), lsh(2)], [{'Target Invisible (PFI)'}, {'Reappearance'}]);
        
        set(lg, 'Location', 'Northeast'); shg
        %% set same dimensions as the raincloud plots (plotted in plotSmallFigs.m)
        set(gcf, 'units', 'normalized', 'position', [.25 .41 .25 .4], 'color', 'w');
        cd(basefol)
        cd('Figures')
        cd('GFX Dynamic RESS')
        print('-dpng', ['final result at f ' num2str(hzis)]);
        
    end
    
    %%
    if plotBARoutput ==1
        plotBARoutput_compareHz
        
    end
