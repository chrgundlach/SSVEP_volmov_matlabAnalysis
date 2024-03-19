%% script in order to calculate time frequency analysis of main experiment



clear all

%% parameters
F.PathIn                = 'F:\work\data\SSVEP_volmov\EEG\RAW';
F.PathInEOG             = 'F:\work\data\SSVEP_volmov\EEG\EOG_Data\';
F.PathOut               = 'F:\work\data\SSVEP_volmov\EEG\SET';
% F.Subjects2Use          = [1:20];
%F.Subjects2Use          = [22 23 25 26 27];
% F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; %based on trial number without blinks
F.Subjects2Use          = [25 26 27];
F.EEGChans              = 1:64;
F.EMGChans              = [71 72];
F.VEOGChans             = [];
F.HEOGChans             = [];
F.ChanlocFile           = 'C:\Users\HP-User\matlab\Skripte\Auswertungsskripte\Analyzer_G\ChanLocs\BioSemi64_8_1020.epf';
F.RespTrigger           = {'30'};
F.RespEpoch             = [-5.5 5.5]; % in s
F.SSVEPTrigger          = {'15'};

F.ExpEpochTrigger       = {'11','12'};
F.eptime                = 360; % time of each block in s
F.Blocks                = 6; % number of blocks

F.FrameRate             = 85;
F.TriggerRate           = 6; % SSVEP trigger frequency (e.g. 6 = every 6s)
F.flickframes           = [1 1 1 0 0 0]; % on-off frames for SSVEP

%F.TFAfreqs              = [5:(1/6):40];
F.TFAfreqs              = [4-(1/3):0.25:45];
F.TFAFlag               = [2]; % 1=wavelet; 2=gabor
F.CSD_flag              = 1; % 0 = no; 1 = yes



%% loop for subjects
for i_sub = 1:numel(F.Subjects2Use)
    % read in BDF File
    
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_exp.bdf ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
    EEG=pop_biosig(sprintf('%s\\VP%02.0f_exp.bdf', F.PathIn,F.Subjects2Use(i_sub)));
    % troubleshooting VP24
    EEG = pop_loadset('filename','VP24_exp.set','filepath','F:\\work\\data\\SSVEP_volmov\\EEG\\RAW\\');
    
    % pop_eegplot(EEG,1,1,1);
    %numel(find(cell2mat({EEG.event.type})==15))
    
    % load EOG data
    EOG = open(sprintf('%s\\VP%02.0f_exp_EOG.mat',F.PathInEOG,F.Subjects2Use(i_sub)));
    
    % rereference
    EEG = pop_reref( EEG, [],'exclude',[F.EEGChans(end)+1:EEG.nbchan] ); % average
    % pop_eegplot(EEG,1,1,1,1)
    
    %% TFA power around button press
    % induced and evoked...requires shift of button press to SSVEP cycle beginning
    
    % index block length
    t.i1=cell2mat({EEG.event.type})==str2num(F.ExpEpochTrigger{1});
    t.i2=cell2mat({EEG.event.type})==str2num(F.ExpEpochTrigger{2});
    %.times = [cell2mat({EEG.event(t.i1).latency})' cell2mat({EEG.event(t.i2).latency})'];
    
    EEG_exp = pop_epoch(EEG, F.ExpEpochTrigger(1), [0 F.eptime], 'epochinfo', 'yes');
    if EEG_exp.trials < 6
        EEG_exp = pop_epoch(EEG, F.ExpEpochTrigger(1), [0 F.eptime-0.5], 'epochinfo', 'yes');
    end
    % pop_eegplot(EEG_exp,1,1,1,1)
    
    % do some checking
    %     checking.trials(i_sub) = EEG_exp.trials;
    %     checking.timing_first{i_sub}=...
    %         cell2mat(cellfun(@(x) EEG_exp.epoch(x).eventlatency(find(cell2mat(EEG_exp.epoch(x).eventtype)==15,1,'first')), num2cell(1:EEG_exp.trials)));
    %     checking.timing_last{i_sub}=...
    %         cell2mat(cellfun(@(x) EEG_exp.epoch(x).eventlatency(find(cell2mat(EEG_exp.epoch(x).eventtype)==15,1,'last')), num2cell(1:EEG_exp.trials)));
    %     checking.timing_num{i_sub}=cellfun(@(x) numel(find(cell2mat(EEG_exp.epoch(x).eventtype)==15)), num2cell(1:EEG_exp.trials));
    %     cellfun(@(x) numel(find(cell2mat(EEG_exp.epoch(x).eventtype)==15)), num2cell(1:6))
    %     cellfun(@(x) EEG_exp.epoch(x).eventlatency(find(cell2mat(EEG_exp.epoch(x).eventtype)==15,1,'first')), num2cell(1:6),'UniformOutput',false)
    
    %% loop for blocks
    % create time vectors for visual stimulation
    t.bl=1;
    for i_bl = 1:F.Blocks
        % trigger times for
        t.t_trig=cell2mat({EEG_exp.epoch(i_bl).eventlatency{cell2mat(EEG_exp.epoch(i_bl).eventtype)==str2num(F.SSVEPTrigger{1})}}); % SSVEP trigger
        t.t_btrig=cell2mat({EEG_exp.epoch(i_bl).eventlatency{cell2mat(EEG_exp.epoch(i_bl).eventtype)==str2num(F.RespTrigger{1})}}); % behavior trigger
        
        t.frametimes=[];
        % create time for each frame (fit between triggers)
        for i_tr = 1:numel(t.t_trig)-1 % 1:numel(t.t_trig)-1
            t.frametimes=[ t.frametimes linspace(t.t_trig(i_tr),t.t_trig(i_tr+1),F.FrameRate*F.TriggerRate+1)];
        end
        % comepensate for missing last SSVEP trigger...add time towards the end
        t.frametimes=[ t.frametimes (0:median(diff(t.frametimes)):(median(diff(t.t_trig))))+t.frametimes(end)];
        t.frametimes=sort(unique(t.frametimes));
        % discard frames outside epoch
        exptiming.frametime{i_bl}=t.frametimes(t.frametimes<EEG_exp.times(end));
        
        % create vector with on-off cycles
        exptiming.framemat{i_bl}=[repmat(F.flickframes,1,floor((numel(exptiming.frametime{i_bl}))/numel(F.flickframes)))];
        exptiming.framemat{i_bl}=[exptiming.framemat{i_bl} F.flickframes(1:(numel(exptiming.frametime{i_bl})-numel(exptiming.framemat{i_bl})))];
        
        % create vector with trigger onset
        exptiming.triggermat{i_bl}= zeros(1,numel(exptiming.framemat{i_bl}));
        exptiming.triggers{i_bl}=  t.t_trig(t.t_trig>=exptiming.frametime{i_bl}(1)&t.t_trig<=exptiming.frametime{i_bl}(end));
        exptiming.triggermat{i_bl}(arrayfun(@(x) find(exptiming.frametime{i_bl}<=x,1,'last'), exptiming.triggers{i_bl}))=1;
        
        % create vector with response onsets
        exptiming.respmat{i_bl}=    zeros(1,numel(exptiming.framemat{i_bl}));
        exptiming.responses{i_bl}=  t.t_btrig(t.t_btrig>=exptiming.frametime{i_bl}(1)&t.t_btrig<=exptiming.frametime{i_bl}(end));
        exptiming.respmat{i_bl}(arrayfun(@(x) find(exptiming.frametime{i_bl}<=x,1,'last'), exptiming.responses{i_bl}))=1;
        
        % function to get only first output: subsref(test(input),struct('type','()','subs',{{1}}))
        %     t.shifts = arrayfun(@(x) subsref(findstr(exptiming.framemat{1,1}(x-numel(F.flickframes)+1:end),[1 1]),struct('type','()','subs',{{1}}))-numel(F.flickframes),...
        %         find(exptiming.respmat{i_bl}==1));
        t.shifts = arrayfun(@(x) subsref(findstr(exptiming.framemat{i_bl}(x-round(numel(F.flickframes)/2):end),...
            F.flickframes(1:end/2)),struct('type','()','subs',{{1}}))-round(numel(F.flickframes)/2)-1,...
            find(exptiming.respmat{i_bl}==1));
        exptiming.respmat_sh{i_bl}=zeros(size(exptiming.respmat{i_bl}));
        exptiming.respmat_sh{i_bl}(find(exptiming.respmat{i_bl})+t.shifts)=1;
        
        % t.tt=[exptiming.frametime{i_bl}; exptiming.framemat{i_bl};  exptiming.triggermat{i_bl}; exptiming.respmat{i_bl}; exptiming.respmat_sh{i_bl}];
    end
    
    %% write new events
    for i_bl = 1:6
        t.ev=find(exptiming.respmat_sh{i_bl});
        for i_ev = 1:numel(t.ev)
            EEG_exp = pop_editeventvals(EEG_exp,'append',{1 [] [] [] 1},...
                'changefield',{2 'epoch' i_bl},...
                'changefield',{2 'latency' exptiming.frametime{i_bl}(t.ev(i_ev))},...
                'changefield',{2 'type' 35});
        end
        % numel(find(cell2mat(EEG_exp.epoch(6).eventtype)==30))
        % numel(find(cell2mat(EEG_exp.epoch(6).eventtype)==35))
    end
    EEG_exp=eeg_checkset(EEG_exp);
    EEG_exp = pop_resample(EEG_exp, 256);
    %% extract new epocs + artifact correction
    % pop_eegplot(EEG_exp,1,1,1)
    
        
    % chanlocs
    EEG_exp.chanlocs=pop_chanedit(EEG_exp.chanlocs,'load',{F.ChanlocFile,'filetype','besa (elp)'}); % Load Channel Locations
    %figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
    EEG_ext = pop_select(EEG_exp,'channel',65:72);
    EEG_exp = pop_select(EEG_exp,'channel',1:64);
    
    %% save
    fprintf(1,'saving dataset %s.set here %s\n', sprintf('VP%02.0f_exp_prep',F.Subjects2Use(i_sub)),F.PathOut)
    pop_saveset(EEG_exp,'filename',sprintf('VP%02.0f_exp_prep',F.Subjects2Use(i_sub)),'filepath',F.PathOut);
    
end