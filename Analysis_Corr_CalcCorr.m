%% script in order to calculate time frequency analysis of main experiment
% requires Analysis_Corr_SavePreprocessData.m to be run previously



clear all

%% parameters
F.PathIn                = 'F:\work\data\SSVEP_volmov\EEG\SET';
F.PathInEOG             = 'F:\work\data\SSVEP_volmov\EEG\EOG_Data\';
F.PathOut               = 'F:\work\data\SSVEP_volmov\EEG\CORR_TFA_RESS_noblinks';
% F.Subjects2Use          = [1:20];
%F.Subjects2Use          = [22 23 25 26 27];
F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; %based on trial number without blinks
% F.Subjects2Use          = [1];
F.EEGChans              = 1:64;
F.EMGChans              = [71 72];
F.VEOGChans             = [];
F.HEOGChans             = [];
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
F.TFA_baseline          = [-3000 -2750];


corr.parameters = {...
    [8 12] [-1000 0] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'vis alpha'...
    [8 12] [-1000 0] {'C3';'CP3'} 'motor alpha';...
    [8 12] [-1000 0] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'vis alpha'...
    [10 14] [-1000 0] {'C3';'CP3'} 'motor alpha';...
    [10 14] [-1000 0] {'C3';'CP3'} 'motor alpha'...
    [18 30] [-1000 0] {'C3';'CP3'} 'motor beta';...
    [10 14] [-1000 0] {'C3';'CP3'} 'motor alpha'...
    [10 14] [1500 2500] {'C3';'CP3'} 'motor alpha';...
    [18 30] [-1000 0] {'C3';'CP3'} 'motor beta'...
    [18 30] [1500 2500] {'C3';'CP3'} 'motor beta';...
    };
corr.parameters2 = {...
    [8 12] [-1000 0] [1500 2500] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'vis alpha diff'...
    [8 12] [-1000 0] [1500 2500] {'C3';'CP3'} 'motor alpha diff';...
    [8 12] [-1000 0] [1500 2500] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'vis alpha diff'...
    [10 14] [-1000 0] [1500 2500] {'C3';'CP3'} 'motor alpha diff';...
    [10 14] [-1000 0] [1500 2500] {'C3';'CP3'} 'motor alpha diff'...
    [18 30] [-1000 0] [1500 2500] {'C3';'CP3'} 'motor beta diff';...
    [8 12] [-1000 0] [1500 2500] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'vis alpha diff'...
    [18 30] [-1000 0] [1500 2500] {'C3';'CP3'} 'motor beta diff'...
    };



%% loop for subjects
for i_sub = 1:numel(F.Subjects2Use)
    %% preprocessing
    % read in set file
    
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_exp_prep.set ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
    EEG = pop_loadset('filename',sprintf('VP%02.0f_exp_prep.set',F.Subjects2Use(i_sub)),'filepath',F.PathIn);
    % load EOG data
    EOG = open(sprintf('%s\\VP%02.0f_exp_EOG.mat',F.PathInEOG,F.Subjects2Use(i_sub)));
    
    % set up parameters for calculation
    if F.CSD_flag == 1 & i_sub == 1 % calculate CSD matrix
        CSD.chanmat=ExtractMontage('C:\Users\HP-User\matlab\software\toolboxes\CSD\resource\10-5-System_Mastoids_EGI129.csd',{EEG.chanlocs.labels}');
        [CSD.G,CSD.H] = GetGH(CSD.chanmat);
    end
    
    switch F.TFAFlag
        case 1
            t.ep_time = F.RespEpoch;
        case 2
            t.ep_time = F.RespEpoch+[-0.5 0.5];
    end
    EEG_expep_n = pop_epoch( EEG, {'30'}, t.ep_time, 'epochinfo', 'yes');
    EEG_expep_sh = pop_epoch( EEG, {'35'}, t.ep_time, 'epochinfo', 'yes');
    
    % get timing (adjusted)
    art.resp_time_n=[];
    for i_tr = 1:EEG_expep_n.trials
        t.ix=EEG_expep_n.epoch(i_tr).eventurevent{dsearchn(cell2mat(EEG_expep_n.epoch(i_tr).eventlatency)',0)};
        art.resp_time_n(i_tr)=EEG_expep_n.urevent(t.ix).latency./EEG_expep_n.srate;
    end
    art.resp_time_sh=[];
    for i_tr = 1:EEG_expep_sh.trials
        t.ix=EEG_expep_sh.epoch(i_tr).eventurevent{dsearchn(cell2mat(EEG_expep_sh.epoch(i_tr).eventlatency)',0)};
        art.resp_time_sh(i_tr)=EEG_expep_sh.urevent(t.ix).latency./EEG_expep_sh.srate;
    end
    
    % discard trials with blinks around movement period [-3 +3] s
    t.rejindex_n = false(1,EEG_expep_n.trials);
    t.rejindex_sh = false(1,EEG_expep_sh.trials);
    % index trials in certain time interval
    art.win_blinks = [-3000 3000];
    for i_tr = 1:EEG_expep_n.trials
        if EOG.EOG.art.blinknum_in_resptrials_n(i_tr)
            t.rejindex_n(i_tr)=any((art.win_blinks(1) <= EOG.EOG.art.blink_time_in_resptrials_n{i_tr}) &...
                (art.win_blinks(2) >= EOG.EOG.art.blink_time_in_resptrials_n{i_tr}));
        end
    end
    for i_tr = 1:EEG_expep_sh.trials
        if EOG.EOG.art.blinknum_in_resptrials_sh(i_tr)
            t.rejindex_sh(i_tr)=any((art.win_blinks(1) <= EOG.EOG.art.blink_time_in_resptrials_sh{i_tr}) &...
                (art.win_blinks(2) >= EOG.EOG.art.blink_time_in_resptrials_sh{i_tr}));
        end
    end
    % some later errors with eeglab if trialnumber < 2
    if diff([sum(t.rejindex_n) numel(t.rejindex_n)])<2
        t.rejindex_n([1 2])=false;
    end
    if diff([sum(t.rejindex_sh) numel(t.rejindex_sh)])<2
        t.rejindex_sh([1 2])=false;
    end
    
    % detrend
    EEG_expep_n = eegF_Detrend(EEG_expep_n);
    EEG_expep_sh = eegF_Detrend(EEG_expep_sh);
    
    %% calculate RESS component
    [RESS.expep_sh] = RESS_Calculate_Component(EEG_expep_sh, F.FrameRate/numel(F.flickframes), 'neighfreqdist', 2, 'neighfreqwidth', 2);
            
    % use spatial filters of shifted data for not shifted data
    EEG_expep_n_noress=EEG_expep_n;
    
    % noress data
    noress_ts = zeros(size(EEG_expep_n_noress.data));
    for i_tr=1:EEG_expep_n_noress.trials
        noress_ts(:,:,i_tr)  = RESS.expep_sh{1}.ComponentMaps(:,1:end-1) * RESS.expep_sh{1}.eigenvectors(:,1:end-1)' * EEG_expep_n_noress.data(:,:,i_tr);
    end
    EEG_expep_n_noress.data=noress_ts;
    
    EEG_surrogate =  pop_select( EEG_expep_n_noress,'time',F.RespEpoch);
    
    % reject trials
    EEG_expep_n = pop_rejepoch(EEG_expep_n, t.rejindex_n, 0);
    EEG_expep_n_noress = pop_rejepoch(EEG_expep_n_noress, t.rejindex_n, 0);
    
    %% calculate current source density transform
    fprintf(1,'\n###\ncalculating CSD transform\n###\n')
    if F.CSD_flag == 1
        for i_tr = 1:EEG_expep_n.trials
            % csd of raw data
            EEG_expep_n.data(:,:,i_tr)= CSDTransform(EEG_expep_n.data(:,:,i_tr), CSD.G, CSD.H);
        end
        for i_tr = 1:EEG_expep_n_noress.trials
            % csd of raw data
            EEG_expep_n_noress.data(:,:,i_tr)= CSDTransform(EEG_expep_n_noress.data(:,:,i_tr), CSD.G, CSD.H);
        end
    end
    
    %% calculate induced TFA
    TFA.electrodes = EEG_expep_n.chanlocs;
    TFA.params.srate = EEG_expep_n.srate;
    TFA.art_win_blinks = art.win_blinks;
    
    %TFA.wave_cycles=7;
    TFA.wave_cycles=[6 (F.TFAfreqs(end))/(F.TFAfreqs(1)/6)];
    
    %filter
    TFA.params.filter = {[] []};
    %TFA.params.filter = {0.5 []};
    %EEG_expep_n = pop_eegfiltnew(EEG_expep_n, TFA.params.filter{1}, TFA.params.filter{2}, EEG_expep_n.pnts, 0, [], 0);
    %EEG_expep_sh = pop_eegfiltnew(EEG_expep_sh, TFA.params.filter{1}, TFA.params.filter{2}, EEG_expep_sh.pnts, 0, [], 0);
    
    % index frequencies to calculate
    t.index = false(size(F.TFAfreqs));
    for i_freq = 1:size(corr.parameters,1)
        t.idx = dsearchn(F.TFAfreqs',corr.parameters{i_freq,1}');
        t.index(t.idx(1):t.idx(2))=true;
    end
    
    TFA.freqs_selected = F.TFAfreqs(t.index);
    TFA.data_induced=nan(sum(t.index),EEG_surrogate.pnts/2, numel(F.EEGChans),EEG_expep_n_noress.trials);
        
    
    switch F.TFAFlag
        case 1
            fprintf(1,'calculating %1.0f TFAs: ',numel(F.EEGChans))
            TFA.type = 'wavelet';
            for i_el=1:numel(F.EEGChans) %loop electrodes !needs to be redone!
                fprintf(1,'%1.0f',mod(F.EEGChans(i_el),10))
                TFA.type = 'wavelet';
                
                [t.tfadata,TFA.t,TFA.f,t.D] = mwavelTFA(squeeze(EEG_expep_n_noress.data(i_el,:,:)),F.TFAfreqs,EEG_expep_n_noress.srate,TFA.wave_cycles,1,3.5);
                for i_tr = 1:EEG_expep_n_noress.trials
                    TFA.data_induced(:,:,i_el,i_tr) = imresize(t.tfadata(t.index,:,i_tr),[sum(t.index),EEG_expep_n_noress.pnts/2]);
                end
                TFA.f = F.TFAfreqs(t.index);
            end
        case 2
            fprintf(1,'calculating %1.0f TFAs:\n',sum(t.index))
            TFA.type = 'gabor';
            TFA.f = F.TFAfreqs(t.index);
            TFA.params.gabor_FWHM_freq = 0.5;
            for i_freq = 1:numel(TFA.f)
                fprintf(1,'%1.0f of %1.0f\n',i_freq, sum(t.index))
                % induced without RESS component
                EEG_Gabor_gab = eegF_Gabor(EEG_expep_n_noress, TFA.f(i_freq),TFA.params.gabor_FWHM_freq);
                EEG_Gabor_gab=pop_select( EEG_Gabor_gab,'time',F.RespEpoch);
                for i_tr = 1:EEG_Gabor_gab.trials
                    TFA.data_induced(i_freq,:,:,i_tr)=imresize(EEG_Gabor_gab.data(:,:,i_tr),[EEG_Gabor_gab.nbchan,EEG_Gabor_gab.pnts/2])';
                end
                t.times = EEG_Gabor_gab.times([1 end]);
            end
            
            
    end
    fprintf(1,'...done\n')
    try 
        TFA.t=imresize(EEG_Gabor_gab.times,[1,EEG_Gabor_gab.pnts/2]);
    catch
        TFA.t=imresize(EEG_expep_n.times,[1,EEG_expep_n_noress.pnts/2]);
    end
    
    %% plotting for checking
%     pl.elec2plot = {'Oz'};
%     pl.elec2plot = {'C3';'CP3'};
%     pl.base=F.TFA_baseline; pl.base_i=dsearchn(TFA.t',pl.base');
%     pl.t2pl = [-4000 4000]; pl.t2pl_i = dsearchn(TFA.t',pl.t2pl');
%     pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
%     
%     pl.data = mean(mean(TFA.data_induced(:,pl.t2pl_i(1):pl.t2pl_i(2),pl.elec2plot_i,:),3),4); pl.ti1 = pl.elec2plot; pl.ti2 = 'raw';  pl.ti3 = 'induced'; pl.clims = [0 max(pl.data(:))];
%     pl.ti1 = pl.elec2plot; pl.ti2 = 'raw';  pl.ti3 = 'induced'; pl.clims = [0 1]*max(abs(pl.data(:)));
%     
%     pl.data = mean(mean(100*(...
%         bsxfun(@rdivide, TFA.data_induced(:,pl.t2pl_i(1):pl.t2pl_i(2),pl.elec2plot_i,:),...
%         mean(TFA.data_induced(:,pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:),2))...
%         -1),3),4);
%     pl.data = mean(mean(...
%         TFA.data_induced(:,pl.t2pl_i(1):pl.t2pl_i(2),pl.elec2plot_i,:)-...
%         repmat(mean(TFA.data_induced(:,pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:),2),[1,numel(pl.t2pl_i(1):pl.t2pl_i(2)) ,1,1])...
%         ,3),4);
%     pl.data = mean((100*(...
%         mean(TFA.data_induced(:,pl.t2pl_i(1):pl.t2pl_i(2),pl.elec2plot_i,:),4)./...
%         repmat(mean(mean(TFA.data_induced(:,pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:),2),4),[1,numel(pl.t2pl_i(1):pl.t2pl_i(2)) ,1,1])...
%         -1)),3);
%     
%     
%     pl.ti1 = pl.elec2plot; pl.ti2 = 'baseline corrected';  pl.ti3 = 'induced'; pl.clims = [-1 1]*max(abs(pl.data(:)));
%     
%     figure;
%     imagesc(TFA.t(pl.t2pl_i(1):pl.t2pl_i(2)),TFA.f,pl.data,pl.clims)
%     set(gca,'YDir','normal')
%     title(sprintf('%s %s tfa %s', pl.ti2, pl.ti3, vararg2str(pl.ti1)), 'FontSize',8)
%     colormap(gca,'jet')
%     colorbar
%     xlabel('time in ms')
%     ylabel('frequency in Hz')
%     hline(14.16667,'m')
%     set(gca,'FontSize',8)
%     
    %% extract data for correlation
    for i_corr = 1:size(corr.parameters,1)
        t.fidx = dsearchn(TFA.f',corr.parameters{i_corr,1}');
        t.tidx = dsearchn(TFA.t',corr.parameters{i_corr,2}');
        t.eidx = logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), corr.parameters{i_corr,3}, 'UniformOutput',false)),1));
        t.bidx = dsearchn(TFA.t',F.TFA_baseline');
        % raw
        corr.xdata{i_corr,1} = squeeze(mean(mean(mean(...
            TFA.data_induced(t.fidx(1):t.fidx(2),t.tidx(1):t.tidx(2),t.eidx,:),1),2),3));
        % baseline corrected
        corr.xdata_bc{i_corr,1} = squeeze(mean(mean(mean(...
            bsxfun(@minus, TFA.data_induced(t.fidx(1):t.fidx(2),t.tidx(1):t.tidx(2),t.eidx,:),...
            mean(TFA.data_induced(t.fidx(1):t.fidx(2),t.bidx(1):t.bidx(2),t.eidx,:),2))...
            ,1),2),3));
        
        t.fidx = dsearchn(TFA.f',corr.parameters{i_corr,5}');
        t.tidx = dsearchn(TFA.t',corr.parameters{i_corr,6}');
        t.eidx = logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), corr.parameters{i_corr,7}, 'UniformOutput',false)),1));
        % raw
        corr.ydata{i_corr,1} = squeeze(mean(mean(mean(...
            TFA.data_induced(t.fidx(1):t.fidx(2),t.tidx(1):t.tidx(2),t.eidx,:),1),2),3));
        % baseline corrected
        corr.ydata_bc{i_corr,1} = squeeze(mean(mean(mean(...
            bsxfun(@minus, TFA.data_induced(t.fidx(1):t.fidx(2),t.tidx(1):t.tidx(2),t.eidx,:),...
            mean(TFA.data_induced(t.fidx(1):t.fidx(2),t.bidx(1):t.bidx(2),t.eidx,:),2))...
            ,1),2),3));
        
        % outlier extraction (+- 4 SD)
        t.ind_raw = ...
            corr.xdata{i_corr,1}>=(mean(corr.xdata{i_corr,1})+4*std(corr.xdata{i_corr,1}))|...
            corr.xdata{i_corr,1}<=(mean(corr.xdata{i_corr,1})-4*std(corr.xdata{i_corr,1}))|...
            corr.ydata{i_corr,1}>=(mean(corr.ydata{i_corr,1})+4*std(corr.ydata{i_corr,1}))|...
            corr.ydata{i_corr,1}<=(mean(corr.ydata{i_corr,1})-4*std(corr.ydata{i_corr,1}));
        t.ind_bc = ...
            corr.xdata_bc{i_corr,1}>=(mean(corr.xdata_bc{i_corr,1})+4*std(corr.xdata_bc{i_corr,1}))|...
            corr.xdata_bc{i_corr,1}<=(mean(corr.xdata_bc{i_corr,1})-4*std(corr.xdata_bc{i_corr,1}))|...
            corr.ydata_bc{i_corr,1}>=(mean(corr.ydata_bc{i_corr,1})+4*std(corr.ydata_bc{i_corr,1}))|...
            corr.ydata_bc{i_corr,1}<=(mean(corr.ydata_bc{i_corr,1})-4*std(corr.ydata_bc{i_corr,1}));
             
        
    end
    
    %% extract data for correlation (differences)
    for i_corr = 1:size(corr.parameters2,1)
        t.fidx = dsearchn(TFA.f',corr.parameters2{i_corr,1}');
        t.tidx = dsearchn(TFA.t',corr.parameters2{i_corr,2}');
        t.tidx2 = dsearchn(TFA.t',corr.parameters2{i_corr,3}');
        t.eidx = logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), corr.parameters2{i_corr,4}, 'UniformOutput',false)),1));
        % raw
        corr.xdata_diff{i_corr,1} = squeeze(mean(mean(mean(...
            TFA.data_induced(t.fidx(1):t.fidx(2),t.tidx2(1):t.tidx2(2),t.eidx,:)-...
            TFA.data_induced(t.fidx(1):t.fidx(2),t.tidx(1):t.tidx(2),t.eidx,:),1),2),3));
        
        t.fidx = dsearchn(TFA.f',corr.parameters2{i_corr,6}');
        t.tidx = dsearchn(TFA.t',corr.parameters2{i_corr,7}');
        t.tidx2 = dsearchn(TFA.t',corr.parameters2{i_corr,8}');
        t.eidx = logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), corr.parameters2{i_corr,9}, 'UniformOutput',false)),1));
        % raw
        corr.ydata_diff{i_corr,1} = squeeze(mean(mean(mean(...
            TFA.data_induced(t.fidx(1):t.fidx(2),t.tidx2(1):t.tidx2(2),t.eidx,:)-...
            TFA.data_induced(t.fidx(1):t.fidx(2),t.tidx(1):t.tidx(2),t.eidx,:),1),2),3));
        
        % outlier extraction (+- 4 SD)
        t.ind_raw_diff = ...
            corr.xdata_diff{i_corr,1}>=(mean(corr.xdata_diff{i_corr,1})+4*std(corr.xdata_diff{i_corr,1}))|...
            corr.xdata_diff{i_corr,1}<=(mean(corr.xdata_diff{i_corr,1})-4*std(corr.xdata_diff{i_corr,1}))|...
            corr.ydata_diff{i_corr,1}>=(mean(corr.ydata_diff{i_corr,1})+4*std(corr.ydata_diff{i_corr,1}))|...
            corr.ydata_diff{i_corr,1}<=(mean(corr.ydata_diff{i_corr,1})-4*std(corr.ydata_diff{i_corr,1}));
        
        
    end
    
    %% preliminary plotting
%     for i_corr = 1:size(corr.parameters,1)
%         figure;
%         pl.xdata = corr.xdata{i_corr,1}(~t.ind_raw);
%         pl.ydata = corr.ydata{i_corr,1}(~t.ind_raw);
%         
%         h.sc1=scatter(pl.xdata , pl.ydata ,'k.');
%         hold on
%         h.sc2=scatter(mean(pl.xdata) , mean(pl.ydata) ,'ro');
%         
%         % xlims
%         %xlim(1.1*[-1 1]*max(abs([pl.xdata;pl.ydata])));
%         %ylim(1.1*[-1 1]*max(abs([pl.xdata;pl.ydata])));
%         xlim([min([pl.xdata;pl.ydata]) max([pl.xdata;pl.ydata])]+(diff([min([pl.xdata;pl.ydata]) max([pl.xdata;pl.ydata])])*[-0.05 0.05]));
%         ylim([min([pl.xdata;pl.ydata]) max([pl.xdata;pl.ydata])]+(diff([min([pl.xdata;pl.ydata]) max([pl.xdata;pl.ydata])])*[-0.05 0.05]));
%         
%         % correlation
%         [r.R,r.P]=corrcoef(pl.xdata,pl.ydata);
%         text(mean(get(gca,'xlim')),max(get(gca,'ylim')),...
%             sprintf('R=%1.3f p=%1.4f\n',r.R(2),r.P(2)), ...
%             'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8);
%         % add fitted line
%         r.fit = polyfit(pl.xdata, pl.ydata,1);
%         t.x=linspace(min(get(gca,'xlim')),max(get(gca,'xlim')),100);
%         plot(t.x,r.fit(1)*t.x+r.fit(2),'-','Color',[0.6,0.6,0.6])
%         
%         t.pos=get(gca,'Position');
%         
%         % title and labels
%         title(sprintf('raw correlation | %s X %s | %1.0f outliers', ...
%             vararg2str(corr.parameters{i_corr,4}),vararg2str(corr.parameters{i_corr,8}),sum(t.ind_raw)),'FontSize',8)
%         xlabel(...
%             sprintf('%s | %1.1f to %1.1f Hz | %1.0f to %1.0f ms\n%s', ...
%             corr.parameters{i_corr,4},corr.parameters{i_corr,1},corr.parameters{i_corr,2},vararg2str(corr.parameters{i_corr,3})),...
%             'FontSize',8)
%         ylabel(...
%             sprintf('%s | %1.1f to %1.1f Hz | %1.0f to %1.0f ms\n%s', ...
%             corr.parameters{i_corr,8},corr.parameters{i_corr,5},corr.parameters{i_corr,6},vararg2str(corr.parameters{i_corr,7})),...
%             'FontSize',8)
%         
%         grid on
%         box on
%     end
%     % baseline corrected
%     for i_corr = 1:size(corr.parameters,1)
%         figure;
%         pl.xdata = corr.xdata_bc{i_corr,1}(~t.ind_bc);
%         pl.ydata = corr.ydata_bc{i_corr,1}(~t.ind_bc);
%         
%         h.sc1=scatter(pl.xdata , pl.ydata ,'k.');
%         hold on
%         h.sc2=scatter(mean(pl.xdata) , mean(pl.ydata) ,'ro');
%         
%         % xlims
%         xlim(1.1*[-1 1]*max(abs([pl.xdata;pl.ydata])));
%         ylim(1.1*[-1 1]*max(abs([pl.xdata;pl.ydata])));
%         
%         % correlation
%         [r.R,r.P]=corrcoef(pl.xdata,pl.ydata);
%         text(mean(get(gca,'xlim')),max(get(gca,'ylim')),...
%             sprintf('R=%1.3f p=%1.4f\n',r.R(2),r.P(2)), ...
%             'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8);
%         % add fitted line
%         r.fit = polyfit(pl.xdata, pl.ydata,1);
%         t.x=linspace(min(get(gca,'xlim')),max(get(gca,'xlim')),100);
%         plot(t.x,r.fit(1)*t.x+r.fit(2),'-','Color',[0.6,0.6,0.6])
%         
%         t.pos=get(gca,'Position');
%         
%         % title and labels
%         title(sprintf('baseline corrected correlation | %s X %s | %1.0f outliers', ...
%             vararg2str(corr.parameters{i_corr,4}),vararg2str(corr.parameters{i_corr,8}),sum(t.ind_bc)),'FontSize',8)
%         xlabel(...
%             sprintf('%s | %1.1f to %1.1f Hz | %1.0f to %1.0f ms\n%s', ...
%             corr.parameters{i_corr,4},corr.parameters{i_corr,1},corr.parameters{i_corr,2},vararg2str(corr.parameters{i_corr,3})),...
%             'FontSize',8)
%         ylabel(...
%             sprintf('%s | %1.1f to %1.1f Hz | %1.0f to %1.0f ms\n%s', ...
%             corr.parameters{i_corr,8},corr.parameters{i_corr,5},corr.parameters{i_corr,6},vararg2str(corr.parameters{i_corr,7})),...
%             'FontSize',8)
%         
%         grid on
%         box on
%     end
%     %% difference plotting
%     for i_corr = 1:size(corr.parameters2,1)
%         figure;
%         pl.xdata = corr.xdata_diff{i_corr,1}(~t.ind_raw_diff);
%         pl.ydata = corr.ydata_diff{i_corr,1}(~t.ind_raw_diff);
%         
%         h.sc1=scatter(pl.xdata , pl.ydata ,'k.');
%         hold on
%         h.sc2=scatter(mean(pl.xdata) , mean(pl.ydata) ,'ro');
%         
%         % xlims
%         %xlim(1.1*[-1 1]*max(abs([pl.xdata;pl.ydata])));
%         %ylim(1.1*[-1 1]*max(abs([pl.xdata;pl.ydata])));
%         xlim([min([pl.xdata;pl.ydata]) max([pl.xdata;pl.ydata])]+(diff([min([pl.xdata;pl.ydata]) max([pl.xdata;pl.ydata])])*[-0.05 0.05]));
%         ylim([min([pl.xdata;pl.ydata]) max([pl.xdata;pl.ydata])]+(diff([min([pl.xdata;pl.ydata]) max([pl.xdata;pl.ydata])])*[-0.05 0.05]));
%         
%         % correlation
%         [r.R,r.P]=corrcoef(pl.xdata,pl.ydata);
%         text(mean(get(gca,'xlim')),max(get(gca,'ylim')),...
%             sprintf('R=%1.3f p=%1.4f\n',r.R(2),r.P(2)), ...
%             'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8);
%         % add fitted line
%         r.fit = polyfit(pl.xdata, pl.ydata,1);
%         t.x=linspace(min(get(gca,'xlim')),max(get(gca,'xlim')),100);
%         plot(t.x,r.fit(1)*t.x+r.fit(2),'-','Color',[0.6,0.6,0.6])
%         
%         t.pos=get(gca,'Position');
%         
%         % title and labels
%         title(sprintf('raw correlation | %s X %s | [%1.0f %1.0f]ms - [%1.0f %1.0f]ms | %1.0f outliers', ...
%             vararg2str(corr.parameters2{i_corr,5}),vararg2str(corr.parameters2{i_corr,10}), corr.parameters2{i_corr,3},...
%             corr.parameters2{i_corr,2}, sum(t.ind_raw)),'FontSize',8)
%         xlabel(...
%             sprintf('%s | %1.1f to %1.1f Hz \n%s', ...
%             corr.parameters2{i_corr,5},corr.parameters2{i_corr,1},vararg2str(corr.parameters2{i_corr,4})),...
%             'FontSize',8)
%         ylabel(...
%             sprintf('%s | %1.1f to %1.1f Hz \n%s', ...
%             corr.parameters2{i_corr,10},corr.parameters2{i_corr,6},vararg2str(corr.parameters2{i_corr,9})),...
%             'FontSize',8)
%         
%         grid on
%         box on
%     end
    
    %% save data
    corr.baseline = F.TFA_baseline;
    corr.TFAtype = TFA.type;
    
    save(sprintf('%s\\VP%02.0f_corrmat.mat',F.PathOut,F.Subjects2Use(i_sub)),'corr')
    
    clear EEG* EOG RESS TFA noress_ts
        
    
end