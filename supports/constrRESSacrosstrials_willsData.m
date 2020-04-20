function [ressEVEC, ressTS, ressSNR, ressMAPS, hz]= constrRESSacrosstrials_willsData(dataIN, peakfreqsall, usechans, neighbour, window,srate, ifreq)


peakwidt  = .5; % FWHM at peak frequency for stimulus

%%%% these need to be calibrated for our closely spaced freqs.


 

% this bit hard coded, 6 s window length
% epochdur = sum(abs([-3 3]));
% timeid = [0:1/srate:epochdur];
% timeid= timeid-3;
 
% tidx= dsearchn(timeid', [window]');

tidx=[1 size(dataIN,2)] ; %use  all timepoints..
        

%% remove any trials woth nan (skips trials in RESS)
badtrials = find(isnan(dataIN(1,[tidx(1)],:)));
tid=[1:size(dataIN,3)];
tid(badtrials)=[];

%we won't use that trial(s) for covariance matrix construction.
% data = data(:,:,tid);


% extract EEG data for RESS
dataIN = squeeze(dataIN(usechans,:,tid));

useneighs=1;


peakfreq1=peakfreqsall(ifreq);
       
        
        ndet = find(peakfreqsall==peakfreq1);
        try            
neighfreqLo = neighbour(ndet).centrefreqs(1);  % distance of neighboring frequencies away from peak frequency, +/- in Hz
neighwidtLo = neighbour(ndet).sd(1);  % FWHM of the neighboring frequencies

neighfreqHi = neighbour(ndet).centrefreqs(2);  % distance of neighboring frequencies away from peak frequency, +/- in Hz
neighwidtHi = neighbour(ndet).sd(2);  % FWHM of the neighboring frequencies
%recommended by Cohen & Gulbinaite, Neuroimage, 2017.
        catch % in case 20Hz          
        end
        
        % FFT parameters
        nfft = ceil( srate/.1 ); % .1 Hz resolution
        hz    = linspace(0,srate,nfft);

        %% RESS begin

        
        % compute covariance matrix at peak frequency
        

        %filter at interest freq
            fdatAt = filterFGx(dataIN,srate,peakfreq1,peakwidt,0);
        
        
        fdatAt = reshape( fdatAt(:,tidx(1):tidx(2),:), length(usechans),[] );
        fdatAt = bsxfun(@minus,fdatAt,mean(fdatAt,2)); %subtract so mean =0, v.important.
        covAt  = (fdatAt*fdatAt')/diff(tidx); %make sure chans by chans
        
        % are we using neighbourhood for reference, or broadband?
        if useneighs==1
            % compute covariance matrix for lower neighbor
            
            try fdatLo = filterFGx(dataIN(:,:,itrial),srate,neighfreqLo,neighwidtLo);
            catch
                fdatLo = filterFGx(dataIN,srate,neighfreqLo,neighwidtLo,0);
            end
            
            fdatLo = reshape( fdatLo(:,tidx(1):tidx(2),:), length(usechans),[] );
            fdatLo = bsxfun(@minus,fdatLo,mean(fdatLo,2));
            covLo  = (fdatLo*fdatLo')/diff(tidx);
            
            % compute covariance matrix for upper neighbor
            try fdatHi = filterFGx(dataIN(:,:,itrial),srate,neighfreqHi,neighwidtHi,0);
            catch
                fdatHi = filterFGx(dataIN,srate,neighfreqHi,neighwidtHi,0);
            end
            
            fdatHi = reshape( fdatHi(:,tidx(1):tidx(2),:), length(usechans),[] );
            fdatHi = bsxfun(@minus,fdatHi,mean(fdatHi,2));
            covHi  = (fdatHi*fdatHi')/diff(tidx);
            
            %%
            [evecs,evals] = eig(covAt,(covHi+covLo)/2); %finds eig vector of sig vs noise.
            
        else %use broadband signal
            %               fdatBB=filterFGx(dataIN,srate,20,20,1);
            fdatBB = reshape( dataIN(:,tidx(1):tidx(2),:), length(usechans),[] );
            fdatBB = bsxfun(@minus,fdatBB,mean(fdatBB,2));
            covBB  = (fdatBB*fdatBB')/diff(tidx);
            
            [evecs,evals] = eig(covAt,(covBB)); %finds eig vector of sig vs noise.
        end
        
        %%
        [~,comp2plot] = max(diag(evals)); % find maximum component
        evecs = bsxfun(@rdivide,evecs,sqrt(sum(evecs.^2,1)));
        % normalize vectors (not really necessary, but OK) - good if diff length
        % epochs/cross ppants.
        
        % extract components and force sign
        % maps = inv(evecs'); % get maps (this is fine for full-rank matrices)
        maps = covAt * evecs / (evecs' * covAt * evecs); % this works either way
        %
        [~,idx] = max(abs(maps(:,comp2plot))); % find biggest component in map
        maps = maps * sign(maps(idx,comp2plot)); % force to positive sign
        
        %%
%         % reconstruct RESS component time series, for this trial.
%         ress_ts1 = [];
%         
%        % reconstruct RESS component time series
        ress_ts1 = zeros(size(dataIN,2),size(dataIN,3));
        for ti=1:size(dataIN,3)
            ress_ts1(:,ti) = evecs(:,comp2plot)'*squeeze(dataIN(:,:,ti));
        end
      ressTS=ress_ts1;
        %
        ressEVEC=evecs(:,comp2plot);
        
%         
        ressxAll = (abs( fft(ress_ts1(tidx(1):tidx(2),:),nfft,1)/diff(tidx) ).^2);
        ressx= nanmean(ressxAll,2);
%         
%         
%         
%         % compute SNR for this trial:
%         [snrR,snrE] = deal(zeros(size(hz)));
%         
%         snrRall=zeros(size(ressxAll));
%         
%         %%% very important, compare again to the correct frequencies 
%         %we have used
%         %low shoulder
% %         
%         skipbinsLo =  dsearchn(hz', [peakfreq1-neighbour(ndet).snrlow(2)]'); % hard-coded in Hz, how many away before start shoulder.        
%         numbinsLo  = dsearchn(hz', [peakfreq1-neighbour(ndet).snrlow(1)]'); % hard-coded in Hz, how many away before start shoulder.        
%         
%         
%         skipbinsHi =  dsearchn(hz', [neighbour(ndet).snrhigh(1)- peakfreq1]'); % hard-coded in Hz, how many away before start shoulder.        
%         numbinsHi  = dsearchn(hz', [neighbour(ndet).snrhigh(2) - peakfreq1]'); % hard-coded in Hz, how many away before start shoulder.        
% 

snrR=zeros(1,length(ressx));
% 
% for hzi=numbinsLo+1:length(hz)-numbinsHi-1
%     numer = ressx(hzi);
%     denom = mean( ressx([hzi-numbinsLo:hzi-skipbinsLo hzi+skipbinsHi:hzi+numbinsHi]) );
%     snrR(hzi) = numer./denom;
%     
% end


%% alternate method for SNR.
% skipbins =  6; % *.1 - make sure it exceeds the Hbw (0.5Hz for 2s window)
% numbins  = 20+skipbins; %  2 Hz, also hard-coded!
% for hzi=numbins+1:length(hz)-numbins-1
%     numer = ressx(hzi);
%     denom = mean( ressx([hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
%     snrR(hzi) = numer./denom;
%     
% end

% plot(hz, snrR)
% % % hold on
% xlim([ 0 25])
% % % % %%
% plot(hz, 10*log10(ressx.^2))
% hold on
% % plot(hz, snrR2)
% xlim([0 25])
% shg
%%
% end
ressSNR=snrR;
        
        % store map
        % also topo for inspection
        map2plot = maps(:,comp2plot);
        map2save=map2plot./max(map2plot);
        
        if ~isreal(map2save)
            error(' ')
        end
        
        %save original trial by trial data for comparison.
        try ressTS(itrial,ifreq,:) = ress_ts1';
            ressSNR(itrial,ifreq,:) = snrR;
            
            ressMAPS(itrial,ifreq,:) = map2save;
        catch
%             ressTS = ress_ts1';
%             ressSNR = snrR;
            
            ressMAPS = map2save;
            
        end
        
        
    


end