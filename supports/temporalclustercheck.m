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
                
                %
                for icl=1:size(clusterSTandEND,1)
                    
                    %start and end are now:
                    % change icl to maxClust if we only want the largest
                    % cluster.
                    STC=sigs(clusterSTandEND(icl,1));
                    ENDC=sigs(clusterSTandEND(icl,2)+1);
                    checktimes =STC:ENDC;
                    observedCV = sum(abs(tvals(checktimes)));
                    % now shuffle condition labels to see if this cluster is
                    % sig (compared to chance).
                    
                    
                    
                    %all trials
                    alltr=reshape(ttestdata, [length(allppants)*2, length(Nt)]);
                    alltrialsvec = 1:size(alltr,1);
                    
                    nshuff=2000;
                    
                    sumTestStatsShuff = zeros(1,nshuff);
                    for irand = 1:nshuff
                        
                        if shuffType==1 %null is that no condition differences.
                            
                            
                            shD=zeros(size(ttestdata));
                            %since this is a within subjects design, we permute
                            %the subjet specific averages within each subject
                            %(as per Maris & Oostenveld (2007).
                            
                            % for each subject, randomly permute the averages.
                            %(Dsub1cond1,Datasub1cond2)
                            for ippant = 1:size(ttestdata,2)
                                
                                if mod(randi(100),2)==0 %if random even number
                                    shD(1,ippant,:) = ttestdata(1,ippant,:); % all time points.
                                    shD(2,ippant,:) = ttestdata(2,ippant,:); % all time points.
                                else
                                    shD(1,ippant,:) = ttestdata(2,ippant,:); % all time points.
                                    shD(2,ippant,:) = ttestdata(1,ippant,:); % all time points.
                                end
                                
                                %                             shD(ipartition,ippant,:) = pdata;
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
                        
                        %                     tvalspertimepoint = zeros(1,length(checktimes));
                        %%
                        % figure(3); clf; plot(squeeze(mean(shD(1,:,:),2))); hold on
                        %                     plot(squeeze(mean(shD(2,:,:),2))); ylim([2.4 3])
                        %%
                        testdata = squeeze(shD(1,:,:)) - squeeze(shD(2,:,:));
                        p=[];
                        for itest = 1:length(checktimes)
                            
                            [~, p(itest), ~,stat]= ttest(testdata(:,checktimes(itest)));
                            
                            tvalspertimepoint(1,itest) = stat.tstat;
                        end
                        
                        % the null hypothesis is that these prob distributions
                        % are exchangeable, so retain this permutation cluster-
                        % level stat.
                        sumTestStatsShuff(1,irand) = sum((tvalspertimepoint));
                        
                        
                    end %repeat nshuff times
                    
                    
                    %is the observed greater than CV?
                    % plot histogram:
                    %%
                    figure(2);
                    
                    clf
                    
                    %                    H=histogram(abs(sort(sumTestStatsShuff)));
                    H=histogram((sort(sumTestStatsShuff)));
                    % fit CDF
                    cdf= cumsum(abs(H.Data))/ sum(abs(H.Data));
                    %the X values (actual CV) corresponding to .01
                    [~,cv05uncorr] = (min(abs(cdf-.95)));
                    %                 [~,cv01uncorr] = (min(abs(cdf-.99)));
                    %                 [~,cv001uncorr] = (min(abs(cdf-.999)));
                    hold on
                    pCV=plot([observedCV observedCV], ylim, ['r-']);
                    
                    p05=plot([H.Data(cv05uncorr) H.Data(cv05uncorr)], ylim, ['k:']);
                    %                 plot([H.Data(cv01uncorr) H.Data(cv01uncorr)], ylim, ['k:']);
                    %                 plot([H.Data(cv001uncorr) H.Data(cv001uncorr)], ylim, ['k:']);
                    legend([pCV p05], {['observed'] ['95%'] })
                    
                    %%
                    
                    
                    if observedCV>H.Data(cv05uncorr)
                        title(['sum tvals = ' num2str(observedCV)]);
                        %              title('Spatial Cluster  Significant!')
                        %                     timeidDYN=tgrm-3;
                        
                        figure(1);
                        if icl==1
                            axis tight;
                        yl=get(gca, 'ylim');
                        end
                        if idtype==1 % place sig markers at top of figure for PFI
%                             sigplace = yl(2)+.7*(diff(yl));
                        sigplace = yl(2)+.02;
                        else
                            sigplace = yl(1)-.2*(diff(yl)); % place at bottom for catch
                        end
                        for itime=checktimes
                            
                            hold on
                            
                            
                            plot(timeidDYN(itime), sigplace, ['*' ],'markersize', 15, 'linewidth', 3, 'color', sigcolor)
                            
                            %                         plot(timeidDYN(itime), 3.25, ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'k')
                            %                             plot(tgrm(itime)-3, sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'm')
                        end
                    end
                    
                    %
                end
            end
            
          
            