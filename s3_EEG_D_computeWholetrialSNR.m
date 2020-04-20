% script to compute whole trial spectrum, and SNR.

try cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA/EEG')
catch
    cd('/Users/MDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA/EEG')
end

basefol=pwd;
%%
pdirs = dir([pwd filesep '*_*' 'EEG']);

allppants=[1,2,4,6,7,9:19]; %

job.crunchacrossppants=0; %saves also within participant folders.


job.plotacrossppants=1; % plots spectrum and SNR at POz
job.plotTOPOplotsperfreq = 1; % plots topoplots of SNR at designated freqs.



if job.crunchacrossppants==1
    params.Fs=250;
    params.tapers=[1 1];
    
    [acrossPPSPEC,acrossPPSNR]=deal(zeros(length(allppants), 64, 8193));
    
    counter=1;
    for ippant = 1:length(allppants)
        cd(basefol)
        cd(pdirs(allppants(ippant)).name)
        
        load(['P' num2str(allppants(ippant)) '_autopreprocd.mat']);
        
        [pSPEC,pSNR]= deal(zeros(64,8193));
        
        for ichan = 1:64
            chdata = squeeze(EEG_preprocd(ichan,:,:));
            
            [s,f]= mtspectrumc(chdata, params);
            
            pSPEC(ichan,:) = squeeze(nanmean(log(s),2));
            
            %% and take SNR
            kernelw = [-1/8 -1/8 -1/8 -1/8 0 0 0 0 0 1  0 0 0 0 0 -1/8 -1/8 -1/8 -1/8];
            %
            pSNR(ichan,:) = conv(pSPEC(ichan,:), 'kernels', 'same');
            
        end
        
        save('WholetrialmSNR_chanXfreq', 'pSPEC','pSPEC', '-append')
        
        acrossPPSPEC(counter,:,:)=pSPEC;
        acrossPPSNR(counter,:,:)=pSNR;
        counter=counter+1;
    end
    %%
    cd(basefol)
    cd('GFX_Pre-RESS')
    %%
    save('GFX_rawSNR_static_wholetrial', 'acrossPPSPEC', 'acrossPPSNR', 'f')
end
%%
if job.plotacrossppants==1
    %%
    cd(basefol)
    cd('GFX_Pre-RESS')
    load('GFX_rawSNR_static_wholetrial');
    %%
    cd ../../
    %%
    cd('Figures')
    cd('SNR spectrum')
    %%
    figure(1);
    clf
    for idtype = 1:2 % plot . SPEC and SNR.
        switch idtype
            case 1
    mSNRchan=squeeze(mean(acrossPPSPEC(:,62,:),1));                
    ylimsare = [2,10];
            case 2
            mSNRchan=squeeze(mean(acrossPPSNR(:,62,:),1)); % chan POz
    ylimsare = [-1,6];
        end
        
        subplot(1,2,idtype)
    plot(f, mSNRchan, 'k')
    hold on
    cols={'b', 'k', 'b', 'k', 'm','m', 'b', 'b', 'k','m','b'};
    counter=1;
    p=[];
    for ifreq=[15,20,30,40,5,35, 45,60,80, 25, 75]
        [~,usef]= min(abs(f-ifreq));
        
        
        
        yis= mSNRchan(usef);
        xis= f(usef);
        p(counter)=plot(xis, yis, 'o', 'markersize', 15, 'linew', 5, 'color', cols{counter}) ;
        %            p(counter)=text(xis, yis, 'o', 'fontsize', 25, 'color', cols{counter}) ;
        
        if ifreq==60
            p(counter)=plot(xis, yis, 'o', 'markersize', 30,'linew', 5, 'color', 'k') ;
        end
        counter=counter+1;
        
        
        
        % significant?
        [~, pvals(ifreq)] = ttest(acrossPPSPEC(:,62,usef), 0); %comparing to zero
        
        %         q=.05/length(hztocheck);
        q=.05/length(f);
        if pvals(ifreq)<q
            tt= plot(f(usef), mSNRchan(usef)+.7 , ['*k'], 'markers', 15,'linew', 2);
        end
        
        
        
    end
    %%
    lg= legend([p(1), p(2), p(5)], {'F1 and harmonics', 'F2 and harmonics', 'Intermodulation'});
    set(lg, 'fontsize', 20)
    
    axis tight
    xlabel('Frequency (Hz)')
    ylabel('log(SNR)')
    set(gcf, 'color', 'w')
    set(gca, 'fontsize', 35)
    %         ylim([0 .4])
    %%
    xlim([ 0 41])
    ylim([ylimsare])
    end % data type.
    %%
    print('-dpng', 'Whole trial SNR, chan POz')
    
    
    if job.plotTOPOplotsperfreq==1;
    
    cd(basefol)
    cd ../
    cd('Figures')
    cd('Topoplots')
    %%
    %topos first.
    getelocs;
    %%
    counter=1;
    colormap('viridis')
    mTOPOs = squeeze(mean(acrossPPSNR,1));
    clf
    titles={'(TG) f1', '(TG) 2f1', '(TG) 3f1',...
        '(BG) f2', '(BG) 2f2', '(BG) 3f2',...
        '(IM) f2-f1', '(IM) 2f2-f1', '(IM) f1 +f2'};
    
    titles={'f1', '2f1', '3f1',...
        'f2', '2f2', '3f2',...
        '(IM) f2-f1', '(IM) 2f2-f1', '(IM) f1 +f2'};
    %     for ifreq=[15,30,45,20,40,60,5,25,35]
    
    %% new shorter version.
    titles={'f1', 'f2', '(IM) f2-f1',...
        '2f1', '2f2', '(IM) f1 +f2'};
    counter=1;
    for ifreq=[15,20, 5, 30,40, 35]
        
        [~,usef]= min(abs(f-ifreq));
        
        
        
        %         topoplot(squeeze(mean(acrossPPspec(:,:,fid),1)), elocs(1:64), 'conv', 'on', 'emarker2', {[64,60], 'o', 'w',15,5});
        pvals = zeros(1,64);
        for ichan=1:64
            [~,pvals(ichan)] = ttest(squeeze(acrossPPSNR(:,ichan,usef))); % compares to zero
        end
        
        q=fdr(pvals);
        pmask = q<1;
        plotD=squeeze(mean(acrossPPSNR(:,:,usef),1));
        
        subplot(2,3,counter)
        topoplot(plotD, elocs(1:64), 'conv', 'on');%, 'pmask', pmask);
        
        
        
        
        c=colorbar;
        %         cm=round(max(plotme));
        %         if cm<1
        %             cm=1;
        %         end
        %         caxis([0 cm])
        if mod(ifreq,20)==0
            caxis([0 3])
        else
            caxis([0 2])
            
        end
        
        %         title([num2str(ifreq) ' Hz'])
        title(titles{counter})
        set(gca, 'fontsize', 25)
        counter=counter+1;
        ylabel(c,'log(SNR)')
    end
    colormap('viridis')
    %
    shg
    %%
    set(gcf, 'color', 'w')
    print('-dpng', 'Topos across all, SNR preRESS, whole trial')
    end
   