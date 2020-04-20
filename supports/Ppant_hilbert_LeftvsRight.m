function Ppant_hilbert_LeftvsRight(basefol,dirs,allppants)


getelocs;
for ifol=allppants%(1:9);
    cd(basefol)
    cd('EEG')
    cd(dirs(ifol).name)
    
    for iPFIvsPMD =1:2
        if iPFIvsPMD==1
            load('ppant_PFI_Epoched_LeftRight');
        else
            load('ppant_PFI_Epoched_LeftRight_PMD');
        end
        
        
        for idatatype=1:4
            
            if iPFIvsPMD==1
                switch idatatype
                    case 1
                        datatmp=ppant_EEG_PFI_Lefton;
                    case 2
                        datatmp=ppant_EEG_PFI_Leftoff;
                    case 3
                        datatmp=ppant_EEG_PFI_Righton;
                    case 4
                        
                        datatmp=ppant_EEG_PFI_Rightoff;
                end
            else
                
                % note that some trials are Nan arrays, as PMD responses
                % have been sorted, those  that failed (e.g. without BP)
                % are now removed as nan
                switch idatatype
                    case 1
                        datatmp=ppant_EEG_PMD_Lefton;
                    case 2
                        datatmp=ppant_EEG_PMD_Leftoff;
                    case 3
                        datatmp=ppant_EEG_PMD_Righton;
                    case 4
                        
                        datatmp=ppant_EEG_PMD_Rightoff;
                end
            end
            
            
            
            
            TF_outgoing=zeros(64,1501);
           
            tstamps = [1:1:size(datatmp,2)]/250 - 3;
            
            for ichan=1:64
                
                chan_data= squeeze(datatmp(:,ichan,:));
                
                %remove mean per channel
                
                allavg=nanmean(chan_data(:,2));
                dataref=chan_data-repmat(allavg,[1,size(chan_data,2)]);
                
                
                %filter at particular range.
                %              [filtdat,empVals] = filterFGx(data,srate,f,fwhm,showplot)
                
                [filtdat,empVals] = filterFGx(dataref,250,10,4,0);
                
%                 [filtdat,empVals] = filterFGx(dataref,250,15,1,0);
                %%
                % now take the hilbert (envelope) across each trial.
                hilbout=zeros(size(filtdat));
                
                %note there is an edge artefact, so rebaseline using
                %the mean -250:0
                baseremtt= dsearchn(tstamps', [-2 -1]');
                for itrial = 1:size(filtdat,1)
                    
                    
                    %absolute for the envelope.
                    hilbenv=abs(hilbert(filtdat(itrial,:)));
                    
                    %remove basline?
                    baseval = mean(hilbenv(baseremtt(1):baseremtt(2)));
                    
                    %                     hilbenv_brem=hilbenv-baseval;
                    
                    %                     hilbout(itrial,:)=hilbenv_brem;
                    
                    hilbout(itrial,:)=hilbenv;
                    
                    % %
                    %                      plot(chan_data(itrial,:)); hold on;
                    %                      plot(filtdat(itrial,:));
                    %                      plot(hilbenv_brem)
                end
                
                
%                 %remove implausible (artefactual trials).
%                 mtr = max(hilbout(:,100:1401),[],2);
%                 rm = find(mtr>100);
%                 hilbout(rm,:) = nan;
%                 
%                 if max(chan_data(:))>100
%                     
%                     figure(9);
%                     subplot(311); imagesc(chan_data); shg;
%                     subplot(312); plot(chan_data');
%                     subplot(313); plot(hilbout');
% %                     error checkplot
%                 end
                
                % avoid nan trials
                TF_outgoing(ichan,:)=squeeze(nanmean(hilbout,1));
                
                
                
                %
                %                 [s,f]= mtspectrumc(chan_data',params);
                %                 mtspectrumOUTgoing(ichan,:)=log(s);
            end
            
            % %%             occipital alpha?
            %             figure(1); clf
            % %             subplot(121); plot(TF_outgoing');
            %             subplot(122);
            %             topoplot(squeeze(mean(TF_outgoing(:,600:900),2)), elocs(1:64));
            %%
            %             figure(2);
            %             for ichan=1:64
            %                 subplot(8,8,ichan);
            %                 plot(f, mtspectrumOUTgoing(ichan,:));
            %                 xlim([0 40]);
            %                 ylim([0 10]);
            %                 title(elocs(ichan).labels)
            %             end
            %
            
            
            
            %%
            
            switch idatatype
                case 1
                    ppant_HILB_EEG_Lefton=TF_outgoing;
                case 2
                    ppant_HILB_EEG_Leftoff=TF_outgoing;
                case 3
                    ppant_HILB_EEG_Righton=TF_outgoing;
                case 4
                    ppant_HILB_EEG_Rightoff=TF_outgoing;
            end
            
        end % idatatype
        
        
        %check that the av hilberts aren't just shifted in time (using
        %baseline).
        %         d1=squeeze(ppant_HILB_EEG_Lefton(42,:));
        %         d2=squeeze(ppant_HILB_EEG_Leftoff(42,:));
        %         plot(d1); hold on; plot(d2,'r'); shg
        %
        %
        %save in ppant fol:
        if iPFIvsPMD==1
            save('ppant_PFI_Epoched_LeftRight',...
                'ppant_HILB_EEG_Rightoff','ppant_HILB_EEG_Righton',...
                'ppant_HILB_EEG_Leftoff','ppant_HILB_EEG_Lefton', '-append');
        else
            save('ppant_PFI_Epoched_LeftRight_PMD',...
                'ppant_HILB_EEG_Rightoff','ppant_HILB_EEG_Righton',...
                'ppant_HILB_EEG_Leftoff','ppant_HILB_EEG_Lefton', '-append');
        end
    end % PFI or PMD?
end % all ppants