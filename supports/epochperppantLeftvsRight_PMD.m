% function epochperppantLeftvsRight_PMD(basefol, dirs, allppants)
cd(basefol)
    cd('Behaviour')
    
    %% load data to determine physical catch timing.
    load('MD_AllBP_perppant')
load('PFI_data');
    %%
    
    window=[-3 3];
    
    srate=250;
    
    allppants=[1,2,4,6,7,9:19]; % new
    
    epochdur = sum(abs(window))*srate;
    onsetc = ceil(epochdur)/2;
    
    
    
    for ifol =allppants
        cd(basefol)
        cd('EEG')
        cd(dirs(ifol).name)
        %load catch epoched data
        load('ppant_Catch_Epoched.mat');
        
        %% for this ppant group as left or right disappearances.
        LCatches = [ppantTrialDatawithDetails(ifol).TrialDetails(:).TL_Catch;ppantTrialDatawithDetails(ifol).TrialDetails(:).BL_Catch ];
        Lany = find(round(mean(LCatches,1)));
        RCatches = [ppantTrialDatawithDetails(ifol).TrialDetails(:).TR_Catch;ppantTrialDatawithDetails(ifol).TrialDetails(:).BR_Catch ];
        Rany = find(round(mean(RCatches,1)));
        
        %Which catch trials were unique to L or Right side (one or two
        %targets?).
        Lunique = setdiff(Lany, Rany);
        Runique = setdiff(Rany, Lany);
        %%
        %start counters per ppant.
        [counterLeft_on,counterLeft_off,counterRight_on, counterRight_off]=deal(ones(1,1));        
        
        %empty variables (per ppant).
        [ppant_EEG_PMD_Lefton,ppant_EEG_PMD_Leftoff]= deal([]);
        [ppant_EEG_PMD_Righton,ppant_EEG_PMD_Rightoff]= deal([]);
        
         
for iLeftvsRight=1:2
    if iLeftvsRight==1
        usetrials = Lunique; 
    else
        usetrials = Runique; 
    end
    
        for tmptrials = 1:length(usetrials)
            
            itrial = usetrials(tmptrials);
            
            trialdata= PFI_only(ifol).Trial(itrial);
            
            
            if trialdata.Goodtrial==1
                % note that even if is a good trial... (PMD response). some
                % nan can be contained in this EEG, since BP were not quick
                % enough, or simultaneous PFI etc.
                
                %disappearances first.
                outgoingPFI= squeeze(ppant_SNREEG_disapBPwithincatch(itrial, :, :));
             
                if iLeftvsRight==1
                    ppant_EEG_PMD_Lefton(counterLeft_on,:,:) = outgoingPFI;
                    counterLeft_on=counterLeft_on+1;
                else
                    ppant_EEG_PMD_Righton(counterRight_on,:,:) = outgoingPFI;
                    counterRight_on=counterRight_on+1;
                end
                             
                outgoingPFI = squeeze(ppant_SNREEG_reapBPaftercatch(itrial, :, :));
                
                if iLeftvsRight==1 % left 1st
                    ppant_EEG_PMD_Leftoff(counterLeft_off,:,:) = outgoingPFI;
                    counterLeft_off=counterLeft_off+1;
                else
                    ppant_EEG_PMD_Rightoff(counterRight_off,:,:) = outgoingPFI;
                    counterRight_off=counterRight_off+1;
                end
             

            end
                            
                            
        end % L or R unique trials
end % L vs Right.
                    
       
    disp(['saving participant ' num2str(ifol) ' counts: '...
        num2str(counterLeft_on) ', ' num2str(counterLeft_off) ',',...
         num2str(counterRight_on) ', ' num2str(counterRight_off) ]);
    savename= ['ppant_PFI_Epoched_LeftRight_PMD'];
    
    
    save(savename, 'ppant_EEG_PMD_Lefton', 'ppant_EEG_PMD_Leftoff',...
        'ppant_EEG_PMD_Righton', 'ppant_EEG_PMD_Rightoff')
    
    disp(['Finished ppant ' num2str(ifol)])
    end %all ppants