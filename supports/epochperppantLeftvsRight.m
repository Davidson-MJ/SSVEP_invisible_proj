function epochperppantLeftvsRight(basefol, dirs, allppants)
cd(basefol)
    cd('Behaviour')
    
    %% load data to determine physical catch timing.
    load('MD_AllBP_perppant')
    load('PFI_data')
    %%
    
    window=[-3 3];
    
    srate=250;
    
    allppants=[1,2,4,6,7,9:19]; % new
    
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
        %left side targets disap
        
        
        %
        %start counters per ppant.
        [counterLeft_on,counterLeft_off,counterRight_on, counterRight_off]=deal(ones(1,1));        
        
        %empty variables (per ppant).
        [ppant_EEG_PFI_Lefton,ppant_EEG_PFI_Leftoff]= deal([]);
        [ppant_EEG_PFI_Righton,ppant_EEG_PFI_Rightoff]= deal([]);
        
        for itrial=1:48
            
            trialdata= PFI_only(ifol).Trial(itrial);            
            
            if trialdata.Goodtrial==1
                
                SNRdata = squeeze(EEG_preprocd(:,:,itrial)); %all channels, 20 and 40 hz.
                
                %first establish the direction of changes (disap /reappear
                %targets).
                
                %for each trial, look at the left and right side separately.
                allBP = trialdata.allBPs;
                
                LeftBP= nansum(allBP([1,3],:),1);
                RightBP= nansum(allBP([2,4],:),1);
                
                for iLeftvsRight=1:2
                    
                    % collect disap/reap for left or right hemisphere:                    
                    % convert the buttons on the opposite side, into Nan.
                    
                    %find instances when RIGHT buttons were pressed.
                    
                    if iLeftvsRight==1
                        tmpPress= RightBP;
                        tmpPress(tmpPress>0)=nan;
                        
                        tmpBP = [LeftBP; tmpPress]; % stack nans with remaining 'good' BP periods.
                        plotcol='r';
                    else
                        tmpPress= LeftBP;
                        tmpPress(tmpPress>0)=nan;
                        
                        %now add for accumBP
                        tmpBP = [RightBP;tmpPress];
                        plotcol='b';
                    end
                    %remove those from our trace.
                    accumBP= sum(tmpBP,1); % this leaves only left or right sided BP.
                    
                    
                    % want to remove transients?
                    excludeTransientlength=12; % 60 frames is 1 second.
                    
                    if excludeTransientlength>0%; %frames
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
                        
%                         hold on;
%                         plot(accumBP, plotcol); shg%check method for isolating left from right.
                        
                    end
                    
                    
                    %change nans to zero, to only find disap/reap for true lateralized PFI
                    accumBP(isnan(accumBP))=0;
                    
                    directionBP = diff(accumBP);
                    directionBPtimes = find(directionBP~=0);
                    directionBPtimes=directionBPtimes+1; %account for diff function
                    %now remove zeros, leaving time of all disap/reap
                    directionBP(directionBP==0)=[];
                    
                    
                    %% NOW collect EEG per remaining data BP.
                    for iPFI = 1:length(directionBPtimes)
                        %Collect event info.
                        timeis = directionBPtimes(iPFI);
                        perceptis = accumBP(directionBPtimes(iPFI));
                        perceptwas = accumBP(directionBPtimes(iPFI)-1);
                        
                        
                        outgoingPFI=[];
                        
                        %disap or reap?
                        if perceptwas<perceptis
                            PFIdir=1; %disappearing, an extra button is being pressed
                        else
                            PFIdir=-1; %Reappearing
                        end
                        
                        %gather EEG time index.
                        [~,timePFI]= min(abs(tt-timeis/60));
                        
                        % ignore events at start/end of block.
                        
                        if  timePFI-onsetc>1 && timePFI+onsetc<size(SNRdata,2)
                            %collect data
                            outgoingPFI= SNRdata(:,timePFI-onsetc:timePFI+onsetc);
                            
                            
                            %store according to location, 
                            %
                            if PFIdir==1 %increasing targets absent
                                
                                if iLeftvsRight==1 % left 1st
                                ppant_EEG_PFI_Lefton(counterLeft_on,:,:) = outgoingPFI;
                                counterLeft_on=counterLeft_on+1;
                                else
                                    ppant_EEG_PFI_Righton(counterRight_on,:,:) = outgoingPFI;
                                    counterRight_on=counterRight_on+1;
                                end
                                
                                
                            else %decreasing number targets absent.
                                
                                
                               if iLeftvsRight==1 % left 1st
                                ppant_EEG_PFI_Leftoff(counterLeft_off,:,:) = outgoingPFI;
                                counterLeft_off=counterLeft_off+1;
                                else
                                    ppant_EEG_PFI_Rightoff(counterRight_off,:,:) = outgoingPFI;
                                    counterRight_off=counterRight_off+1;
                               end
                            end
                            
                            
                        end % within trial limits.
                    end % all PFI events
                    
                 
                    
                end % left/right
            end % goodtrials
        end %alltrials
        
      
    disp(['saving participant ' num2str(ifol) ' counts: '...
        num2str(counterLeft_on) ', ' num2str(counterLeft_off) ',',...
         num2str(counterRight_on) ', ' num2str(counterRight_off) ]);
    savename= ['ppant_PFI_Epoched_LeftRight'];
    
    
    save(savename, 'ppant_EEG_PFI_Lefton', 'ppant_EEG_PFI_Leftoff',...
        'ppant_EEG_PFI_Righton', 'ppant_EEG_PFI_Rightoff')
    
    disp(['Finished ppant ' num2str(ifol)])
    end %all ppants