%% checks timing for presentation

F.path =        'D:\work\data\SSVEP_volmov\behavior\';
F.sub =         'VP01';
F.blocknum  =   6;
F.RefRate =     85.005; % refresh rate

% open logfile
logfile = open(sprintf('%s%s_timing.mat',F.path,F.sub));

% check timing
for i_bl = 1:F.blocknum
%     logfile.timing.experiment{i_bl}.Missed
%     abs(diff(logfile.timing.experiment{i_bl}.FlipTimestamp))>(1/F.RefRate)
    
    
%     figure; plot(abs(diff(logfile.timing.experiment{i_bl}.FlipTimestamp))); hline((1/F.RefRate),'k')
%     figure; plot(abs(diff(logfile.timing.experiment{i_bl}.VBLTimestamp))); hline((1/F.RefRate),'k')
%     figure; plot(abs(diff(logfile.timing.experiment{i_bl}.StimulusOnsetTime))); hline((1/F.RefRate),'k')
%     figure; plot(abs(diff(logfile.timing.experiment{i_bl}.Missed))); hline((1/F.RefRate),'k')
    figure; plot(abs(diff(logfile.timing.experiment{i_bl}.Beampos))); 
end