% s3_B_rereferenceEEG
%new basic preprocessing script.
addpath('/Users/MattDavidson/Desktop/SSVEP-Wills')
cd('/Users/MattDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA/EEG')
basefol=pwd;

%%
getelocs % load channel files
dbstop if error

allppants=[1,2,4,6,9:16,18]; %
% goodppants=[1,2,4,6,7,9:19]; % after behavioural based rejection.

dirs = dir([pwd filesep '*_*' 'EEG']);

%%
for ifol = allppants
    cd(basefol)
    cd(dirs(ifol).name)
    
    %load Epoched EEG;
    load(['P' num2str(ifol) 'RawEEG'])
    
    %first trial is a dud, though not always recorded with EEG,
    %retain just 48 in total    
    outgoing=zeros(48,64,15001);    
    
    for itrial=1:48 %since first was a practice.
        
        if size(EpochedEEGdata,1)==49
            realtrial=itrial+1; % if the firs was recorded with EEG, skip.
        elseif  size(EpochedEEGdata,1)==48
            realtrial=itrial;
        end
        tmp=squeeze(EpochedEEGdata(realtrial,1:64,:));
        re_tmp=reref(tmp);        
        %                     %downsample
        re_tmp = downsample(re_tmp',4);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        outgoing(itrial,:,:) = re_tmp';
    end
    
    %perform laplacian
    outgoing = permute(outgoing, [2 3 1]);
    EEG_preprocd=laplacian_perrinX(outgoing, [elocs(1:64).X], [elocs(1:64).Y], [elocs(1:64).Z]);
    
    save(['P' num2str(ifol) '_autopreprocd'], 'EEG_preprocd')
    
    disp(['Fin ppant ' num2str(ifol)])
   
    
end
