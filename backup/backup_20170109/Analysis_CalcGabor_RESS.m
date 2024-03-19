%% script in order to calculate time frequency analysis of main experiment
% based on RESS components



clear all

%% parameters
F.PathIn                = 'F:\work\data\SSVEP_volmov\EEG\RAW';
F.PathOut               = 'F:\work\data\SSVEP_volmov\EEG\Gabor';
F.Subjects2Use          = [1];
% F.Subjects2Use          = [1];
F.EEGChans              = 1:64;
F.EMGChans              = [71 72];
F.VEOGChans             = [];
F.HEOGChans             = [];
F.ChanlocFile           = 'C:\Users\HP-User\matlab\Skripte\Auswertungsskripte\Analyzer_G\ChanLocs\BioSemi64_8_1020.epf';
F.RespTrigger           = {'30'};
F.RespEpoch             = [-4.5 4.5]; % in s
F.RespEpoch             = [-5 5]; % in s
F.SSVEPTrigger          = {'15'};

F.ExpEpochTrigger       = {'11','12'};
F.eptime                = 359.5; % time of each block in s
F.Blocks                = 6; % number of blocks

F.FrameRate             = 85;
F.TriggerRate           = 6; % SSVEP trigger frequency (e.g. 6 = every 6s)
F.flickframes           = [1 1 1 0 0 0]; % on-off frames for SSVEP

Gabor.parameters        = {...
    [14.16667 14.16667] {'RESS'} 'induced' [0 0 0] 2 'vis ind SSVEP' 'RESS';...
    [14.16667 14.16667] {'RESS'} 'evoked' [0.4 0.4 0.4] 2 'vis evo SSVEP' 'RESS';...
    [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'induced' [1 0 0] 1 'vis alpha' 'NoRESS';...
    [15 30] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'induced' [1 0 1] 1 'vis beta' 'NoRESS';...
    [10 14] {'C3';'CP3'} 'induced' [0 1 0] 1 'mot alpha' 'NoRESS';...
    [15 30] {'C3';'CP3'} 'induced' [0 1 1] 1 'mot beta' 'NoRESS';...
    };
Gabor.parameters        = {...
    [14.16667 14.16667] {'RESS'} 'induced' [0 0 0] 2 'vis ind SSVEP' 'RESS';...
    [14.16667 14.16667] {'RESS'} 'evoked' [0.4 0.4 0.4] 2 'vis evo SSVEP' 'RESS';...
    };
Gabor.baseline          = [-2750 -2250];
Gabor.baseline          = [-3250 -2750];

Gabor.FWHM_freq         = 0.5;

ress.SSVEP_freq         = 85/6;

%% loop for subjects
for i_sub = 1:numel(F.Subjects2Use)
    % read in BDF File
    
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_exp.bdf ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
    EEG=pop_biosig(sprintf('%s\\VP%02.0f_exp.bdf', F.PathIn,F.Subjects2Use(i_sub)));
    % pop_eegplot(EEG,1,1,1);
    
    % chanlocs
    EEG.chanlocs=pop_chanedit(EEG.chanlocs,'load',{F.ChanlocFile,'filetype','besa (elp)'}); % Load Channel Locations
    %figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
    
    % rereference
    EEG = pop_reref( EEG, [],'exclude',[F.EEGChans(end)+1:EEG.nbchan] ); % average
    % pop_eegplot(EEG,1,1,1,1)
    
    %% Gabor/Hilbert timecours
    % induced and evoked...requires shift of button press to SSVEP cycle beginning
    
    % index block length
    t.i1=cell2mat({EEG.event.type})==str2num(F.ExpEpochTrigger{1});
    t.i2=cell2mat({EEG.event.type})==str2num(F.ExpEpochTrigger{2});
    %.times = [cell2mat({EEG.event(t.i1).latency})' cell2mat({EEG.event(t.i2).latency})'];
    
    EEG_exp = pop_epoch(EEG, F.ExpEpochTrigger(1), [0 F.eptime], 'epochinfo', 'yes');
    % pop_eegplot(EEG_exp,1,1,1,1)
    
    %% loop for blocks
    % create time vectors for visual stimulation
    t.bl=1;
    for i_bl = 1:F.Blocks
        % trigger times for
        t.t_trig=cell2mat({EEG_exp.epoch(i_bl).eventlatency{cell2mat(EEG_exp.epoch(i_bl).eventtype)==str2num(F.SSVEPTrigger{1})}}); % SSVEP trigger
        t.t_btrig=cell2mat({EEG_exp.epoch(i_bl).eventlatency{cell2mat(EEG_exp.epoch(i_bl).eventtype)==str2num(F.RespTrigger{1})}}); % behavior trigger
        
        t.frametimes=[];
        % create time for each frame (fit between triggers)
        for i_tr = 1:numel(t.t_trig)-1
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
    end
    % resample
    EEG_exp=pop_select(EEG_exp,'channel',1:64);
    EEG_exp = pop_resample(EEG_exp, 256);
    %% extract new epocs
    % pop_eegplot(EEG_exp,1,1,1)
    EEG_expep_sh = pop_epoch( EEG_exp, {'35'}, F.RespEpoch, 'epochinfo', 'yes');
    EEG_expep = pop_epoch( EEG_exp, {'30'}, F.RespEpoch, 'epochinfo', 'yes');
    %EEG_expep = pop_epoch( EEG_exp, {'15'}, F.RespEpoch, 'epochinfo', 'yes');
    
    EEG_expep = eegF_Detrend(EEG_expep);
    EEG_expep_sh = eegF_Detrend(EEG_expep_sh);
    
    %% calculate RESS component
    [RESS.expep_sh] = RESS_Calculate_Component(EEG_expep_sh, ress.SSVEP_freq, 'neighfreqwidth', 0.5);
    
    EEG_expep_sh_ress = RESS.expep_sh{1,1}.EEG_out;
    EEG_expep_sh_noress = RESS.expep_sh{1,1}.EEG_out_noRESS;
    
    % use spatial filters of shifted data for not shifted data
    EEG_expep_ress=EEG_expep;
    EEG_expep_noress=EEG_expep;
    
    EEG_expep_ress = pop_select( EEG_expep_ress,'channel',1);
    EEG_expep_ress=pop_chanedit(EEG_expep_ress, 'changefield',{1 'labels' 'RESS'},...
        'changefield',{1 'theta' ''},'changefield',{1 'radius' ''},'changefield',...
        {1 'X' ''},'changefield',{1 'Y' ''},'changefield',{1 'Z' ''},'changefield',{1 'sph_theta' ''},...
        'changefield',{1 'sph_phi' ''},'changefield',{1 'sph_radius' ''});
    
    % ress data
    ress_ts = zeros(size(EEG_expep_ress.data));
    for i_tr=1:EEG_expep_ress.trials
        ress_ts(1,:,i_tr) = RESS.expep_sh{1}.eigenvectors(:,RESS.expep_sh{1}.CompNum)'*squeeze(EEG_expep.data(:,:,i_tr));
    end
    EEG_expep_ress.data=ress_ts;
    
    % noress data
    noress_ts = zeros(size(EEG_expep_noress.data));
    for i_tr=1:EEG_expep_noress.trials
        noress_ts(:,:,i_tr)  = RESS.expep_sh{1}.ComponentMaps(:,1:end-1) * RESS.expep_sh{1}.eigenvectors(:,1:end-1)' * EEG_expep_noress.data(:,:,i_tr);
    end
    EEG_expep_noress.data=noress_ts;
    
    
    %% loop for parameters
    for i_gab = 1:size(Gabor.parameters,1)
        
        % switch for induced vs evoked
        switch Gabor.parameters{i_gab,3}
            case 'induced'
                switch Gabor.parameters{i_gab,7}
                    case 'RESS'
                        EEG_Gabor=EEG_expep_ress;
                    case 'NoRESS'
                        EEG_Gabor=EEG_expep_noress;
                end
            case 'evoked'
                switch Gabor.parameters{i_gab,7}
                    case 'RESS'
                        EEG_Gabor=EEG_expep_sh_ress;
                    case 'NoRESS'
                        EEG_Gabor=EEG_expep_sh_noress;
                end
        end
        
        % preliminary artifact correction
        t.elec2use_i=find(logical(sum(cell2mat(cellfun(@(x) strcmp({EEG_Gabor.chanlocs.labels},x), Gabor.parameters{i_gab,2}, ...
            'UniformOutput',false)),1)));
        t.gabdata = [];
        
        EEG_Gabor_t=pop_select(EEG_Gabor,'channel',t.elec2use_i);
        
        % index trials below/above criteria
%         t.data = squeeze(EEG_Gabor_t.data(1,:,:));
%         Gabor.trials_used{i_gab,i_el,i_sub}=find(sum((t.data>(mean(t.data(:))+6*std(t.data(:))))|(t.data<(mean(t.data(:))-6*std(t.data(:)))))==0);
%         EEG_Gabor_t=pop_select(EEG_Gabor_t,'trial',Gabor.trials_used{i_gab,i_el,i_sub});
        
        % average across trials for evoked first
        if strcmp(Gabor.parameters{i_gab,3},'evoked')
            EEG_Gabor_t.data(:,:,1)=mean(EEG_Gabor_t.data,3);
            EEG_Gabor_t=pop_select(EEG_Gabor_t,'trial',1);
        end
        
        % calculate gabor transformation
        Gabor.freqs{i_gab}=Gabor.parameters{i_gab,1}(1):Gabor.FWHM_freq:Gabor.parameters{i_gab,1}(2);
        for i_freq = 1:numel(Gabor.freqs{i_gab})
            EEG_Gabor_gab = eegF_Gabor(EEG_Gabor_t,Gabor.freqs{i_gab}(i_freq),Gabor.FWHM_freq);
            t.gabdata(:,i_freq,:)= squeeze(mean(EEG_Gabor_gab.data,3));
        end
        Gabor.data(:,i_gab,i_sub) = squeeze(mean(mean(t.gabdata,1),2));
        
    end
    
    
 
end
Gabor.subjects = F.Subjects2Use;
Gabor.time = EEG_Gabor_gab.times;
Gabor.electrodes=EEG.chanlocs;
Gabor.data_bc=Gabor.data - ...
    repmat(mean(Gabor.data(eeg_time2points(Gabor.baseline(1),Gabor.time):eeg_time2points(Gabor.baseline(2),Gabor.time),:,:),1),...
    [size(Gabor.data,1), 1, 1]);

%% save

Gabor.savetime=datestr(now);
fprintf(1,'|| saving file ||  %s\\exp_gabor_RESS.mat ||\n', F.PathOut)
save(sprintf('%s\\exp_gabor_RESS.mat', F.PathOut),'Gabor')

%% load data
% load(sprintf('%s\\exp_gabor_RESS.mat', F.PathOut))
% load(sprintf('%s\\exp_gabor.mat', F.PathOut))

% Gabor.data_bc=((Gabor.data ./ ...
%     repmat(mean(Gabor.data(eeg_time2points(Gabor.baseline(1),Gabor.time):eeg_time2points(Gabor.baseline(2),Gabor.time),:,:),1),...
%     [size(Gabor.data,1), 1, 1]))-1)*100;

%% plot all lines of interest into one raphics
pl.parameters = Gabor.parameters;
pl.subs2use = F.Subjects2Use;
pl.xlim = [-2250 2500];
pl.xlim = [-4500 4500];

pl.toponum = [ceil(sqrt(size(pl.parameters,1))) round(sqrt(size(pl.parameters,1)))];
pl.plpos = repmat([1:pl.toponum(2)*3],pl.toponum(1),1)+repmat((0:pl.toponum(2)*4:pl.toponum(2)*4*(pl.toponum(1)-1))',1,pl.toponum(2)*3);
pl.topopos = repmat([pl.toponum(2)*3+1:pl.toponum(2)*4],pl.toponum(1),1)+repmat((0:pl.toponum(2)*4:pl.toponum(2)*4*(pl.toponum(1)-1))',1,pl.toponum(2));

% statistical analysis
for i_dat = 1:size(pl.parameters)
    [tt.h tt.p tt.ci tt.stats]=ttest(squeeze(Gabor.data_bc(:,i_dat,:))');
    Gabor.bc_stats.p(:,i_dat)=tt.p;
    Gabor.bc_stats.t(:,i_dat)=tt.stats.tstat;
    Gabor.bc_stats.df(:,i_dat)=tt.stats.df;
end

figure;
subplot(pl.toponum(1),pl.toponum(2)*4,pl.plpos(:))
for i_pl = 1:size(pl.parameters)
    % index time
    [t.t t.ind3]=min(abs(Gabor.time-pl.xlim(1)));
    [t.t t.ind4]=min(abs(Gabor.time-pl.xlim(2)));
    
    pl.data = squeeze(mean(Gabor.data_bc(t.ind3:t.ind4,i_pl,pl.subs2use),3));
    pl.data = squeeze(mean(Gabor.data(t.ind3:t.ind4,i_pl,pl.subs2use),3));
    
    % actual plot
    h.pl(i_pl)=plot(Gabor.time(t.ind3:t.ind4),pl.data,'Color',pl.parameters{i_pl,4},'LineWidth',pl.parameters{i_pl,5});
    hold on
    
    % add significance values
    pl.tdata = nan(size(pl.data));
    try
        pl.tdata(Gabor.bc_stats.p(t.ind3:t.ind4,i_pl)<.05)=pl.data(Gabor.bc_stats.p(t.ind3:t.ind4,i_pl)<.05);
    end
    plot(Gabor.time(t.ind3:t.ind4),pl.tdata,'Color',pl.parameters{i_pl,4},'LineWidth',pl.parameters{i_pl,5}+2)
    
end
xlim(pl.xlim)
set(gca,'ylim',[-1.1 1.1].*max(cellfun(@(x) max(abs(x)), get(get(gca,'Children'),'YData'))))
legend(h.pl,pl.parameters(:,6),'location','NorthWest','FontSize', 8)
ylabel('amplitude in \muV')
xlabel('time in ms')
grid on

for i_pl = 1:size(pl.parameters,1)
    subplot(pl.toponum(1),pl.toponum(2)*4,pl.topopos(i_pl))
    pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({Gabor.electrodes.labels},x), pl.parameters{i_pl,2}, 'UniformOutput',false)),1));
    topoplot([],Gabor.electrodes(1:64),'style','blank','electrodes', 'on', 'emarker2', {find(pl.elec2plot_i),'o','r',2});
    title(pl.parameters{i_pl,6},'Color',pl.parameters{i_pl,4})
end
