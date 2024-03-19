%% script to analyse EOG data
% in relation to button presses
% look at SSVEPs around blinks
% different RESS settings for VP 17 


clearvars

%% parameters
F.PathIn                = 'D:\work\data\SSVEP_volmov\EEG\RAW';
F.PathOut               = 'D:\work\data\SSVEP_volmov\EEG\EOG_Data';
% F.Subjects2Use          = [1:20];
% F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; %based on trial number without blinks
% F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27];
F.Subjects2Use          = [1];
F.EEGChans              = 1:64;
F.EMGChans              = [71 72];
F.VEOGChans             = [65 66];
F.HEOGChans             = [67 68];
F.ChanlocFile           = 'C:\Users\psy05cvd\Dropbox\work\matlab\Auswertungsskripte\Analyzer_G\ChanLocs\BioSemi64_8_1020.epf';
F.RespTrigger           = {'30'};
F.RespEpoch             = [-6 6]; % in s
F.EOGEpoch              = [-3.5 3.5]; % in s
F.SSVEPTrigger          = {'15'};
F.BlinkTrigger          = {'50'};

F.ExpEpochTrigger       = {'11','12'};
F.eptime                = 360; % time of each block in s
F.Blocks                = 6; % number of blocks

F.FrameRate             = 85;
F.TriggerRate           = 6; % SSVEP trigger frequency (e.g. 6 = every 6s)
F.flickframes           = [1 1 1 0 0 0]; % on-off frames for SSVEP



%% loop for subjects
for i_sub = 1:numel(F.Subjects2Use)
    % read in BDF File
    
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_exp.bdf ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
%     EEG=pop_biosig(sprintf('%s\\VP%02.0f_exp.bdf', F.PathIn,F.Subjects2Use(i_sub)));
    if F.Subjects2Use(i_sub) ~= 24
        EEG=pop_readbdf(sprintf('%s\\VP%02.0f_exp.bdf', F.PathIn,F.Subjects2Use(i_sub)),[],F.EMGChans(end)+1,[]);
    else
        EEG = pop_loadset('filename','VP24_exp.set','filepath',F.PathIn);
    end
    
    % troubleshooting for subject 24
%     EEG = pop_loadset('filename','VP24_exp.set','filepath','F:\\work\\data\\SSVEP_volmov\\EEG\\RAW\\');
%     t.idx = strcmp({EEG.event.type},'boundary');
%     EEG = pop_editeventvals(EEG,'delete',find(t.idx));
%     for i_ev = 1:numel(EEG.event)
%         EEG.event(i_ev).type = str2double(EEG.event(i_ev).type);
%     end
    
    % pop_eegplot(EEG,1,1,1);
    %numel(find(cell2mat({EEG.event.type})==15))
    
    % rereference
    EEG = pop_reref( EEG, [],'exclude',[F.EEGChans(end)+1:EEG.nbchan] ); % average
    % pop_eegplot(EEG,1,1,1,1)
    
    EEG = pop_resample(EEG, 256);
    
    %% epoch blocks
    EEG_exp = pop_epoch(EEG, F.ExpEpochTrigger(1), [0 F.eptime], 'epochinfo', 'yes');
    if EEG_exp.trials < 6
        EEG_exp = pop_epoch(EEG, F.ExpEpochTrigger(1), [0 F.eptime-0.5], 'epochinfo', 'yes');
    end
    % pop_eegplot(EEG_exp,1,1,1,1)
    
    %% find blinks
    % bipolarize data
    EEG_exp_EOG=eegF_Bipolarize(EEG_exp);
    EEG_exp_EOG=pop_select(EEG_exp_EOG,'channel',{'VEOG','HEOG'});
    EEG_exp_EOG_f = pop_eegfiltnew(EEG_exp_EOG,0.5, 0, 8*EEG_exp_EOG.srate, 0, [], 0);
    %figure; plot(EEG_exp_EOG_f.times,EEG_exp_EOG_f.data(1,:,1));
    for i_bl = 1:F.Blocks
        [blinks.peaks blinks.peaklocs] = findpeaks(EEG_exp_EOG_f.data(1,:,i_bl),'MinPeakDistance',EEG_exp_EOG_f.srate/2,'MinPeakHeight',100);
        %figure; plot(EEG_exp_EOG_f.times,EEG_exp_EOG_f.data(1,:,i_bl)); 
        %hold on; plot(EEG_exp_EOG_f.times(blinks.peaklocs),EEG_exp_EOG_f.data(1,blinks.peaklocs,i_bl),'or')
        
        % write events of blinks
        for i_ev = 1:numel(blinks.peaklocs)
            EEG_exp = pop_editeventvals(EEG_exp,'append',{1 [] [] [] 1},...
                'changefield',{2 'epoch' i_bl},...
                'changefield',{2 'latency' blinks.peaklocs(i_ev)/EEG_exp_EOG_f.srate*1000},... % in seconds
                'changefield',{2 'type' F.BlinkTrigger{1}});
        end
    end
    EEG_exp=eeg_checkset(EEG_exp);
    % pop_eegplot(EEG_exp,1,1,1)
    
    %% extract blink times
    BLNK.totalsamples = EEG_exp.pnts*EEG_exp.trials;
    BLNK.totaltime = EEG.times(end);
    BLNK.trials = EEG_exp.trials;
    BLNK.samplespertrial = EEG_exp.pnts;
    BLNK.blink_samples = [EEG_exp.event([EEG_exp.event.type]==str2num(F.BlinkTrigger{1})).latency];
    BLNK.blink_times = BLNK.blink_samples./EEG_exp.srate.*1000; % in ms
    BLNK.resp_samples = [EEG_exp.event([EEG_exp.event.type]==str2num(F.RespTrigger{1})).latency];
    BLNK.resp_times = BLNK.resp_samples./EEG_exp.srate.*1000; % in ms
    
    BLNK.blink_trialindex = []; BLNK.blink_times_intrial = [];
    BLNK.resp_trialindex = []; BLNK.resp_times_intrial = [];
    for i_tr = 1:EEG_exp.trials
        t.idx = cell2mat([EEG_exp.epoch(i_tr).eventtype])==str2num(F.BlinkTrigger{1});
        BLNK.blink_trialindex = [BLNK.blink_trialindex repmat(i_tr,1,sum(t.idx))];
        BLNK.blink_times_intrial = [BLNK.blink_times_intrial cell2mat(EEG_exp.epoch(i_tr).eventlatency(t.idx))];
        
        t.idx = cell2mat([EEG_exp.epoch(i_tr).eventtype])==str2num(F.RespTrigger{1});
        BLNK.resp_trialindex = [BLNK.resp_trialindex repmat(i_tr,1,sum(t.idx))];
        BLNK.resp_times_intrial = [BLNK.resp_times_intrial cell2mat(EEG_exp.epoch(i_tr).eventlatency(t.idx))];
    end
    BLNK.blink_samples_intrial = round(BLNK.blink_times_intrial.*EEG_exp.srate./1000);
    BLNK.resp_samples_intrial = round(BLNK.resp_times_intrial.*EEG_exp.srate./1000);
    
    BLNK.blink_samplediffs_intrial = diff(BLNK.blink_samples_intrial);
    BLNK.blink_samplediffs_intrial(BLNK.blink_samplediffs_intrial<0)=[];
    BLNK.blink_timediffs_intrial = diff(BLNK.blink_times_intrial);
    BLNK.blink_timediffs_intrial(BLNK.blink_timediffs_intrial<0)=[];
    
    BLNK.resp_samplediffs_intrial = diff(BLNK.resp_samples_intrial);
    BLNK.resp_samplediffs_intrial(BLNK.resp_samplediffs_intrial<0)=[];
    BLNK.resp_timediffs_intrial = diff(BLNK.resp_times_intrial);
    BLNK.resp_timediffs_intrial(BLNK.resp_timediffs_intrial<0)=[];
    
    % all measures in one dataframe
    BLNK.df.onset_samples = [BLNK.blink_samples BLNK.resp_samples];
    BLNK.df.onset_times = [BLNK.blink_times BLNK.resp_times];
    BLNK.df.onset_samples_intrial = [BLNK.blink_samples_intrial BLNK.resp_samples_intrial];
    BLNK.df.onset_times_intrial = [BLNK.blink_times_intrial BLNK.resp_times_intrial];
    BLNK.df.onset_trialindex = [BLNK.blink_trialindex BLNK.resp_trialindex];
    BLNK.df.onset_type = [repmat({'blink'},1,numel(BLNK.blink_samples)) repmat({'resp'},1,numel(BLNK.resp_samples))];
    
    
    %% visualize data
%     % with gramm
%     g=gramm('x',BLNK.df.onset_samples,'color',BLNK.df.onset_type);
%     g.geom_raster();
%     figure('Position',[100 100 800 550]);
%     g.draw();
%     
%     % manually
%     figure;
%     subplot(4,1,1)
%     t.respidx = strcmp(BLNK.df.onset_type,'resp');
%     h.l1 = line([BLNK.df.onset_times(t.respidx)./1000;BLNK.df.onset_times(t.respidx)./1000],...
%         [ones(1,sum(t.respidx))+0.025; ones(1,sum(t.respidx))+0.5], ...
%         'Color', [0.8 0.1 0.1 0.5],'LineWidth',1.5);
%     hold on
%     t.blnkidx = strcmp(BLNK.df.onset_type,'blink');
%     h.l2 = line([BLNK.df.onset_times(t.blnkidx)./1000;BLNK.df.onset_times(t.blnkidx)./1000],...
%         [ones(1,sum(t.blnkidx))-0.5; ones(1,sum(t.blnkidx))-0.025], ...
%         'Color', [0.2 0.4 1 0.5],'LineWidth',1.5);
%     set(gca,'YTickLabels',[],'YTick',[],'XTickLabels',[],'XTick',[])
%     xlim([0,BLNK.totaltime/1000])
%     ylim([0.45,1.55])
% %     xlabel('time in s')
% %     legend([h.l1(1) h.l2(1)], {'finger taps';'blinks'},'Location','southoutside','Orientation','horizontal')
%     title('distribution of events across experimental run')
%     
%    subplot(4,1,[2 3 4])
%    histogram(BLNK.df.onset_times(t.respidx)./1000,30,'BinLimits',[0 BLNK.totaltime/1000],...
%        'FaceColor', [0.8 0.1 0.1])
%    hold on
%    histogram(BLNK.df.onset_times(t.blnkidx)./1000,30,'BinLimits',[0 BLNK.totaltime/1000],...
%        'FaceColor', [0.2 0.4 1])
%    legend({'finger taps';'blinks'},'Location','southoutside','Orientation','horizontal')
% %    title('distribution of events across experimental run')
%    xlim([0,BLNK.totaltime/1000])
%    xlabel('time in s')
%    
%    % histogram of difference between blinks
%    figure;
%    subplot(2,1,1)
%    histogram(BLNK.resp_timediffs_intrial./1000,50,'FaceColor', [0.8 0.1 0.1])
%    xlabel('time in s')
%    title('time differences between finger taps')
%    subplot(2,1,2)
%    histogram(BLNK.blink_timediffs_intrial./1000,50,'FaceColor', [0.2 0.4 1])
%    xlabel('time in s')
%    title('time differences between blinks')
    
    
       
    %% save
    BLNK.savetime=datestr(now);
    fprintf(1,'|| saving file ||  %s\\VP%02.0f_exp_BLNKRESP.mat ||\n', F.PathOut,F.Subjects2Use(i_sub)')
    save(sprintf('%s\\VP%02.0f_exp_BLNKRESP.mat', F.PathOut,F.Subjects2Use(i_sub)), 'BLNK')   
    
end