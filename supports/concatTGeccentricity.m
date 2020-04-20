cd('/Users/MatthewDavidson/Desktop/SSVEP-feedbackproject/AA_ANALYZED DATA/Behaviour')
load('Raw_Behavioral_AllPP.mat', 'AllTrialsAllPPBottomLeft'); 
%%
targX=[];
ntrials=48;
nppants = size(AllTrialsAllPPBottomLeft,1)/ntrials;
for ippant = 1:nppants;
    
    targX(ippant) = AllTrialsAllPPBottomLeft(ippant*48,12);
end
%%
disp([[1:nppants]',targX'])
%%
figure(1);
clf

%convert to DVA:


% FOR MBI EEG DOWNSTAIRS:
display.dist=50;
display.width = 51;
display.resolution = 1920;
TGdva= zeros(size(targX));
for ippant = 1:nppants
pix = targX(ippant); %pixel width of TGs.

TGdva(ippant)=pix2angle(display,pix);
end
%%
clf
hg=histogram(TGdva, ceil(nppants));
hg.FaceColor = 'b';
set(gca, 'xtick', 0:0.6:6);
xlabel({['Calibrated target eccentricity'];['from center (\circ visual angle)']});
ylabel('Participant count')
set(gcf, 'color', 'w');
set(gca, 'fontsize', 25);
xlim([2.9 6])
hg=histogram(round(TGdva,2), 25);
hg.FaceColor = 'b';
xlabel({['Calibrated target eccentricity (d.v.a.)']});
ylabel('# participants')
set(gcf, 'color', 'w');
%
set(gca, 'fontsize', 25, 'Xminortick', 'on','xtick', 3:.5:7);
shg
%%
cd ../Figures
cd('Behavioural Bar Plots')
%%
print('-dpng', 'Target eccentricity per participant')