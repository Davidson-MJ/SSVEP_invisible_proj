% PFI_Jackknife_procedure
% confirm latency differences for F1 and F2 (during PFI).


clear all; close all;

%need to load the PFI data:
% copy-paste of BP_SSVEP script (PFI and Catch).

%background dirs
% addpath('/Users/MattDavidson/Desktop/SSVEP-Wills')
% cd('/Users/matthewDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA')
cd('/Users/mDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA')

basefol=pwd;
clearvars -except basefol allppants
dbstop if error
getelocs
cd(basefol)
cd('EEG')
cd('GFX_EEG-RESS')
%
%     cd('newplots-MD')
colsare={'b' , 'k','b','k','b','g','m'}; % blue for tg, black for BG.
hzare={'TG (f1)' , 'BG (f2)','TG (2f1)','BG (2f2)','TG (3f1)','TGBG','IM (f2-f1)'}; % blue for tg, black for BG.
linetypes={'--', '-', [],[],'--'};
legendis=[];
lgc=1;
%     clf
figure(1); hold on;
clearvars yl;

useMiller=0;
useinterp=1;

%% first as sanity check, plot the group results
[Overall_rawsamp, OverallSigdiff]=deal([]);
for hzis=1:2 %for each frequency
    
    switch hzis
        case 1
            load('GFX_PFIperformance_withSNR_15_min0_RESS')
            %rename for next part of script
            storeacrossPpant_onsetSNR15=storeacrossPpant_onsetSNR;
            storeacrossPpant_offsetSNR15=storeacrossPpant_offsetSNR;
            col='b';
        case 2
            
            load('GFX_PFIperformance_withSNR_20_min0_RESS')
            %rename for next part of script
            storeacrossPpant_onsetSNR20=storeacrossPpant_onsetSNR;
            storeacrossPpant_offsetSNR20=storeacrossPpant_offsetSNR;
            col='m';
    end
    
    %%        Jackknife leave one out
    
    useONSETs = squeeze(nanmean(storeacrossPpant_onsetSNR,2)); %participant averages.
    useOFFSETs = squeeze(nanmean(storeacrossPpant_offsetSNR,2));
    
    ttestdata=[];
    if useinterp==1
        % resample to estimate the latencies at ms level:
        
        %find nearest SNR/ closest time point:
        %interpolate to increase resolution and find closest point (ms)
        %
        samps = timeidDYN(1):1/1000:timeidDYN(end);
        %
        data1= pchip(timeidDYN, useONSETs,samps);
        data2= pchip(timeidDYN, useOFFSETs,samps);
        
    else
        data1=useONSETs;
        data2=useOFFSETs;
    end
    
    %which method?
    if useMiller==1 % calculate using 0.5*peak.
        tmp_Monset = squeeze(mean(data1,1));
        tmp_Moffset = squeeze(mean(data2,1));
        %zero this trace:
%         tmp_use = tmp_Monset - tmp_Monset(1);
        %or use difference\
        tmp_use = tmp_Monset-tmp_Moffset;
        tmp_use = tmp_use - tmp_use(1);
        peakat = 0.9* max(tmp_use); % 90% criterion for response locked ERPs (see Miller et al, 1998 psychophys).
        
        
        % now find closest point
        diffV= abs(tmp_use-peakat);
        if useinterp==1
        [~,Xmark] = find(diffV<.005);
         %we want minimum crossover point:
        Xmarkmin=min(Xmark);
        else
            Xmarkmin = dsearchn(diffV(1:7)', min(diffV));
        end
        
       
        if useinterp==1
            STC = samps(Xmarkmin);
        else
            STC= timeidDYN(Xmarkmin);
        end
        OverallSigdiff(hzis) = STC;
        
        Overall_rawsamp(hzis) = Xmarkmin;
    else
        % calculate first sig time point:
        
        %use orig data.
        ttestdata(1,:,:)= data1;
        ttestdata(2,:,:)= data2;
        
        %              %check for sig
        pvals=zeros(1,size(ttestdata,3));
        tvals=zeros(1,size(ttestdata,3));
        for itime = 1:size(ttestdata,3)
            
            [h,pvals(itime),~,stat]=ttest(ttestdata(1,:,itime), ttestdata(2,:,itime));
            
            
            tvals(itime)= stat.tstat;
        end
        sigs=find(pvals<.05);
        %find max cluster.
        % check for cluster of sigs:
        vect1 = diff(sigs);
        v1 = (vect1(:)==1);
        d = diff(v1);
        clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
        [~,maxClust] = max(clusterSTandEND(:,2)-clusterSTandEND(:,1));
        
        % first point in largest cluster
        minSIG=sigs(clusterSTandEND(maxClust,1));
        if useinterp==1
            STC=samps(minSIG);
        else
            STC=timeidDYN(minSIG);
        end
        
        OverallSigdiff(hzis) = STC;
    end
    
    
    
    
    %%%%%% Sanity check plotting
    % new mean SNR:
    tmp_Monset = squeeze(nanmean(data1,1));
    tmp_Moffset = squeeze(nanmean(data2,1));
    
    if useMiller~=1
        figure(1);
        %             subplot(211)
        if useinterp~=1
        plot(timeidDYN, tmp_Monset, col); hold on;
        plot(timeidDYN, tmp_Moffset, col); hold on;
        else
            plot(samps, tmp_Monset, col); hold on;
        plot(samps, tmp_Moffset, col); hold on;
        end
        plot([STC, STC], ylim, ['-'], 'color', 'r');
        title(['min sig at ' num2str(STC)])
    else
        figure(1);
        if useinterp==1
            plot(samps, tmp_use, 'k')
            
        else
            plot(timeidDYN, tmp_use, 'k')
        end
        hold on;
        plot([STC,STC], ylim, ['k:']);
        
    end
    
end
%%


    %preallocate times....
    [LOO_firstSIG,LOO_firstrawsamp]=deal(zeros(size(useONSETs,1),2));
    
    
    nppantsA=1:size(useONSETs,1); % use all participants
    
    for iJK = 1:length(nppantsA)
        
        %remove this participant for averaging.
        useppants= nppantsA;
        useppants(iJK)=[];
        
        
        for hzis=1:2 %for each frequency
            
            %             if iJK==1 % just load the first time.
            switch hzis
                case 1
                    storeacrossPpant_onsetSNR=storeacrossPpant_onsetSNR15;
                    storeacrossPpant_offsetSNR=storeacrossPpant_offsetSNR15;
                    col='b';
                case 2
                    
                    storeacrossPpant_onsetSNR=storeacrossPpant_onsetSNR20;
                    storeacrossPpant_offsetSNR=storeacrossPpant_offsetSNR20;
                    col='k';
            end
            
            %%        Jackknife leave one out
            %
            % need to compute cross-over points of average timecourses, using leave-one-out procedure.
            
            
            useONSETs = squeeze(nanmean(storeacrossPpant_onsetSNR,2)); %participant averages.
            useOFFSETs = squeeze(nanmean(storeacrossPpant_offsetSNR,2));
            
            %Compute first significat point for this subset of ppants:
            
            subuseONSETs= squeeze(useONSETs(useppants,:));
            subuseOFFSETs= squeeze(useOFFSETs(useppants,:));
            %
            if useinterp==1
                % resample to estimate the latencies at ms level:
                
                %find nearest SNR/ closest time point:
                %interpolate to increase resolution and find closest point (ms)
                %
                samps = timeidDYN(1):1/1000:timeidDYN(end);
                %
                data1= pchip(timeidDYN, subuseONSETs,samps);
                data2= pchip(timeidDYN, subuseOFFSETs,samps);
                xlabis = samps;
            else
                data1=subuseONSETs;
                data2=subuseOFFSETs;
                xlabis = timeidDYN;
            end
            
            %which method?
            if useMiller==1 % calculate using 0.5*peak.
                
                tmp_Monset = squeeze(mean(data1,1));
                tmp_Moffset = squeeze(mean(data2,1));
                %zero this trace:
%                 tmp_use = tmp_Monset - tmp_Monset(1);
                tmp_use = tmp_Monset - tmp_Moffset;
                tmp_use=tmp_use - tmp_use(1);
                
                peakat = 0.9* max(tmp_use);
                
                % now find closest point
                diffV= abs(tmp_use-peakat);
                if useinterp==1
                    [~,Xmark] = find(diffV<.005);
                    %we want minimum crossover point:
                    Xmarkmin=min(Xmark);
                else
                    Xmarkmin = dsearchn(diffV(1:7)', min(diffV));
                end
                
                
                if useinterp==1
                    STC = samps(Xmarkmin);
                    xlabis = samps;
                else
                    STC= timeidDYN(Xmarkmin);
                    xlabis=timeidDYN;
                end
                
                LOO_firstSIG(iJK,hzis)= STC;
                LOO_firstrawsamp(iJK, hzis)=Xmarkmin;
            else
                % calculate first sig time point:
                ttestdata=[];
                %use orig data.
                ttestdata(1,:,:)= data1;
                ttestdata(2,:,:)= data2;
                
                %              %check for sig
                pvals=zeros(1,size(ttestdata,3));
                tvals=zeros(1,size(ttestdata,3));
                for itime = 1:size(ttestdata,3)
                    
                    [h,pvals(itime),~,stat]=ttest(ttestdata(1,:,itime), ttestdata(2,:,itime));
                    
                    
                    tvals(itime)= stat.tstat;
                end
                sigs=find(pvals<.05);
                %find max cluster.
                % check for cluster of sigs:
                vect1 = diff(sigs);
                v1 = (vect1(:)==1);
                d = diff(v1);
                clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
                [~,maxClust] = max(clusterSTandEND(:,2)-clusterSTandEND(:,1));
                
                % first point in largest cluster
                minSIG=sigs(clusterSTandEND(maxClust,1));
                if useinterp==1
                    STC=samps(minSIG);
                else
                    STC=timeidDYN(minSIG);
                end
                
                LOO_firstSIG(iJK,hzis)= STC;
            end
            
            %
            
            %% plot also the mean SNR, and first SIG for sanity check.
            figure(2); subplot(4,4, iJK)
            if useMiller==1
                plot(xlabis, tmp_use)
                hold on;
                plot([STC,STC], ylim, 'color', col);
            else
                plot(xlabis, squeeze(mean(data1,1))); hold on
                plot(xlabis, squeeze(mean(data2,1)));
                plot([STC, STC], ylim, ['-'], 'color', 'r');
            end
            %%
            
        end
                title(['diff is ' num2str(LOO_firstSIG(iJK,2)-LOO_firstSIG(iJK,1))]);
        
        
        
    end
    %% SUMMARY STATS, if using Miller Method
    if useMiller==1
    %From Miller et al., Psychophysiology, 1998. The l-o-o variability
    %can be used to estimate standard error, using eqn (2) of their paper:
    
    %se= sqrt((n-1/n) * sum((J - K)^2)); where:
    % J = estimated value in subsample
    % K = average of values in subsample
    % N=subsample size
    
    N=size(LOO_firstSIG,1);
    % we are interested in the difference in our sample beween
    % f1 and f2, when slopes reached 0.5 max amp. This difference:
    
       
%     J= LOO_firstSIG(:,2)-LOO_firstSIG(:,1); %  F2 and F1.
% D= diff(OverallSigdiff);  %calculated  above    

J= LOO_firstrawsamp(:,1)-LOO_firstrawsamp(:,2); %  F2 and F1.
D= Overall_rawsamp(1)-Overall_rawsamp(2)    
    K= mean(J);
    
    StErrHzfromJK =  sqrt(((N-1)/N)  *   sum((J-K).^2));
    
    %note that using this se estimate, we can compute inferential
    %statistics.
    
    %eqn 3 from Miller et al., 1998;
    %jack-knifed t-stat = D/se;
    
    %D=difference in overall sample (true difference in
    %sig), se is computed above.
    
    
    tstatJK= D/StErrHzfromJK;

    pval = tcdf(tstatJK, N-1)
    end
    %% plot ditribution of results
    
    
    figure(3); clf
    subplot(211)
    h1=raincloud_plot(LOO_firstSIG(:,1),'box_on', 1, 'color', [0,0,1], 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
        'box_col_match', 0);
    %increase scatter size:
    h1{2}.SizeData=100;
    h1{2}.LineWidth=2;
    h2=raincloud_plot(LOO_firstSIG(:,2),'box_on', 1, 'color', [0,0,0], 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
        'box_col_match', 1);
    axis tight
    h2{2}.SizeData=100;
    h2{2}.LineWidth=2;
    
    set(gca, 'fontsize', 20, 'ytick', []);
    xlim([-1.5 0.5])
    xlabel('time from button press, [s]')
    
    title({['First significant SNR points,'];['Jackknife (N-1) subsamples']})
    legend([h1{1} h2{1}], 'Target (f1)', 'Surround (f2)')
    %
    subplot(212);
    h1=raincloud_plot((LOO_firstSIG(:,2)-LOO_firstSIG(:,1)),'box_on', 1, 'color', [1,0,0], 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
        'box_col_match', 0);
    %increase scatter size:
    h1{2}.SizeData=100;
    h1{2}.LineWidth=2;
    
    title('Latency difference (f2-f1), each subsample')
    axis tight
    xlim([-1.5 0])
    hold on; plot([0 0], ylim)
    set(gca, 'fontsize', 20);
    xlabel('seconds')
    
 

%%
