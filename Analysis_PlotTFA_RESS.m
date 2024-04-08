%% plot previously calculated TFA data

    
clearvars

%% parameters
F.PathIn                = 'D:\work\data\SSVEP_volmov\EEG\TFA_Gabor_RESS';
% F.PathIn                = 'D:\work\data\SSVEP_volmov\EEG\TFA_Gabor_RESS_noblinks';
% F.PathIn                = 'F:\work\data\SSVEP_volmov\EEG\TFA_Gabor_RESS_noblinks2'; % latest used
F.PathIn                = 'D:\work\data\SSVEP_volmov\EEG\TFA_Gabor_RESS_noblinks_CSD'; % gabor of 1Hz FWHM
% F.PathIn                = 'D:\work\data\SSVEP_volmov\EEG\TFA_Gabor_RESS_noblinks_CSD2'; % gabor of 0.5Hz FWHM
F.PathIn                = 'D:\work\data\SSVEP_volmov\EEG\TFA_Gabor_RESS_noblinks_CSD3'; % gabor of 1Hz FWHM [no blinks -3 3]
F.PathIn                = 'E:\work\data\SSVEP_volmov\EEG\TFA_Gabor_RESS_noblinks_CSD4'; % gabor of 1Hz FWHM [no blinks -3.5 3.5]
F.PathIn                = 'E:\work\data\SSVEP_volmov\EEG\TFA_Gabor_RESS_noblinks_CSD5'; % gabor of 1Hz FWHM [no blinks -3.5 3.5] RESS altered
F.PathIn                = 'N:\AllgPsy\experimental_data\2016_SSVEP_volmov\EEG\TFA_Gabor_RESS_noblinks_CSD5'; % gabor of 1Hz FWHM [no blinks -3.5 3.5] RESS altered for analysis on Auswerterechner

% F.Subjects2Use          = [1:4];
F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; % based on trial number for each subject
% F.Subjects2Use          = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 25 26 27];

% TFA.baseline            = [-3000 -2750];
% TFA.baseline            = [-2750 -2500];
% TFA.baseline            = [-2250 -2000];
TFA.baseline            = [-3250 -3000];
% TFA.baseline            = [-4000 -3500];
% TFA.baseline            = [-4500 -4000];



pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];
%% loop for subjects
for i_sub=1:numel(F.Subjects2Use)
    %%
    % read in TFA data
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_exp_tfa.mat ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
    
    temp.tfa = open(sprintf('%s\\VP%02.0f_exp_tfa.mat',F.PathIn,F.Subjects2Use(i_sub)));
    
    % convert to single
    temp.tfa.TFA.data_induced = single(temp.tfa.TFA.data_induced);
    temp.tfa.TFA.data_evoked = single(temp.tfa.TFA.data_evoked);
    temp.tfa.TFA.data_RESS_induced = single(temp.tfa.TFA.data_RESS_induced);
    temp.tfa.TFA.data_RESS_evoked = single(temp.tfa.TFA.data_RESS_evoked);
    
    % preallocate memory
    if i_sub == 1
        TFA.data_induced = single(nan([size(temp.tfa.TFA.data_induced) numel(F.Subjects2Use)]));
        TFA.data_evoked = TFA.data_induced;
        TFA.data_induced_bc = TFA.data_induced;
        TFA.data_evoked_bc = TFA.data_induced;
        TFA.data_RESS_induced = single(nan([size(temp.tfa.TFA.data_RESS_induced) numel(F.Subjects2Use)]));
        TFA.data_RESS_evoked = TFA.data_RESS_induced ;
        TFA.data_RESS_induced_bc = TFA.data_RESS_induced ;
        TFA.data_RESS_evoked_bc = TFA.data_RESS_induced ;
        TFA.time = temp.tfa.TFA.t;
        TFA.frequency = temp.tfa.TFA.f;
        TFA.electrodes = temp.tfa.TFA.electrodes;
        TFA.electrodes_RESS = temp.tfa.TFA.electrodes_RESS;
        TFA.RESS_map = nan([length(temp.tfa.TFA.RESS.map) numel(F.Subjects2Use)]);
        TFA.RESS_map_signchan = nan([length(temp.tfa.TFA.RESS.map) numel(F.Subjects2Use)]);
        [TFA.RESS_snr_ind, TFA.RESS_snr_evo] = deal(nan(1,numel(F.Subjects2Use)));
        TFA.srate = temp.tfa.TFA.params.srate/2;
    end
    
    % subject by subject
    TFA.data_induced(:,:,:,i_sub) = temp.tfa.TFA.data_induced; % induced data
%     TFA.data_induced_bc(:,:,:,i_sub) = bsxfun(@minus, temp.tfa.TFA.data_induced, ...
%         mean(temp.tfa.TFA.data_induced(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2));
    TFA.data_induced_bc(:,:,:,i_sub) = 100*((bsxfun(@rdivide, temp.tfa.TFA.data_induced, ...
        mean(temp.tfa.TFA.data_induced(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2)))-1);
   
    TFA.data_evoked(:,:,:,i_sub) = temp.tfa.TFA.data_evoked; % evoked data
%     TFA.data_evoked_bc(:,:,:,i_sub) = bsxfun(@minus, temp.tfa.TFA.data_evoked, ...
%         mean(temp.tfa.TFA.data_evoked(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2));
    TFA.data_evoked_bc(:,:,:,i_sub) = 100*((bsxfun(@rdivide, temp.tfa.TFA.data_evoked, ...
        mean(temp.tfa.TFA.data_evoked(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2)))-1);
    
    TFA.data_RESS_induced(:,:,i_sub) = temp.tfa.TFA.data_RESS_induced; % induced data
%     TFA.data_RESS_induced_bc(:,:,i_sub) = bsxfun(@minus, temp.tfa.TFA.data_RESS_induced, ...
%         mean(temp.tfa.TFA.data_RESS_induced(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2));
    TFA.data_RESS_induced_bc(:,:,i_sub) = 100*((bsxfun(@rdivide, temp.tfa.TFA.data_RESS_induced, ...
        mean(temp.tfa.TFA.data_RESS_induced(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2)))-1);
    
    TFA.data_RESS_evoked(:,:,i_sub) = temp.tfa.TFA.data_RESS_evoked; % evoked data
%     TFA.data_RESS_evoked_bc(:,:,i_sub) = bsxfun(@minus, temp.tfa.TFA.data_RESS_evoked, ...
%         mean(temp.tfa.TFA.data_RESS_evoked(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2));
    TFA.data_RESS_evoked_bc(:,:,i_sub) = 100*((bsxfun(@rdivide, temp.tfa.TFA.data_RESS_evoked, ...
        mean(temp.tfa.TFA.data_RESS_evoked(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2)))-1);
    
    try TFA.trialnum_raw_induced(i_sub)=numel(temp.tfa.TFA.art.resp_time_n);
    end
    try TFA.trialnum_raw_evoked(i_sub)=numel(temp.tfa.TFA.art.resp_time_sh);
    end
    try TFA.trialnum_discarded_induced(i_sub)=numel(temp.tfa.TFA.art.SCADS_Trials2Del_n);
    end
    try TFA.trialnum_discarded_evoked(i_sub)=numel(temp.tfa.TFA.art.SCADS_Trials2Del_sh);
    end
    TFA.totaltrials(i_sub)=temp.tfa.TFA.alltrials; % number of total trials
    %TFA.trials_analyzed_evoked(:,i_sub)=(cellfun(@(x) numel(x),temp.tfa.TFA.trials_evoked)./temp.tfa.TFA.alltrials)*100; % percentage analysed trials
    
    TFA.RESS_map(:,i_sub) = temp.tfa.TFA.RESS.map(:,end);
    TFA.RESS_map_signchan(:,i_sub) = temp.tfa.TFA.RESS.map_signchanged(:,end);
    TFA.RESS_snr_ind(i_sub)=temp.tfa.TFA.RESS.SNR_ind(end);
    TFA.RESS_snr_evo(i_sub)=temp.tfa.TFA.RESS.SNR_evo(end);
    
    clear temp
    
end

%% actual plotting no RESS data
% figure; topoplot([],TFA.electrodes, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);


% plotting parameters
% pl.elec2plot = {'Oz'};
% pl.elec2plot = {'PO3';'POz';'PO4';'O1';'Oz';'O2';'Iz'}; % steady state I
% pl.elec2plot = {'POz';'O1';'Oz';'O2';'I1';'Iz';'I2'}; % steady state II
% pl.elec2plot = {'PO4';'O2';'PO8';'P8';'P10';'I2'}; % vis alpha
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';% vis alpha II
% pl.elec2plot = {'POz'}; % vis alpha II
pl.elec2plot = {'C3';'CP3'};  sav.chan_add = 'MotorSmall';% motor alpha/beta I
% pl.elec2plot = {'C3';'CP3';'C5';'CP5'};  sav.chan_add = 'MotorLarge';% motor alpha/beta II
% pl.elec2plot = {'F1';'F2';'Fz'}; sav.chan_add = 'Frontal1'; % frontal
% pl.elec2plot = {'FP1';'FP2'}; sav.chan_add = 'Frontal2'; % frontal
% pl.elec2plot = {'CP1';'CPz';'CP2';'P1';'Pz';'P2'}; sav.chan_add = 'Parietal'; % parietal
% pl.elec2plot = {'F5';'F3';'F1';'FC3'}; sav.chan_add = 'PreMotor'; % pre-motor
% pl.elec2plot = {'FC3';'FC3';'FC1';'FCz';'Fz';'F5';'F3';'F1';'AF7';'AF3';'AFz';'FP1';'FPz'};  sav.chan_add = 'CentroFrontalLeft'; % pre-motor
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

% pl.xlims=[-2500 2500]; % index time 2 plot
pl.xlims=[-3250 3250]; % index time 2 plot
[t.t t.ind1]=min(abs(TFA.time-pl.xlims(1)));
[t.t t.ind2]=min(abs(TFA.time-pl.xlims(2)));

% pl.flims = TFA.frequency([1 end]); % index frequency 2 plot
% pl.flims = [3 TFA.frequency(end)]; % index frequency 2 plot
% pl.flims = [4 TFA.frequency(end)]; % index frequency 2 plot
% pl.flims = [4 35];
pl.flims = [3 35];
% pl.flims = [8 44];
[t.t t.find1]=min(abs(TFA.frequency-pl.flims(1)));
[t.t t.find2]=min(abs(TFA.frequency-pl.flims(2)));

pl.subs2use = 1:numel(F.Subjects2Use);
% pl.subs2use = [1:14 16 :20];
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 13 14 17 18 19 20]; % subjects discarded due to low trial number version 1
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20]; % subjects discarded due to low trial number version 2
clear figs fig

pl.fre2index = [8 14 15 30];

%%%%%%%%%%%%%
% raw induced no RESS
pl.data_ind = squeeze(nanmean(nanmean(TFA.data_induced(:,:,pl.elec2plot_i,pl.subs2use),3),4));
pl.data_evo = squeeze(nanmean(nanmean(TFA.data_evoked(:,:,pl.elec2plot_i,pl.subs2use),3),4));
t.clims1=[0 1]*max(max(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2)));
t.clims2=[0 1]*max(max(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2)));
% t.clims1=[min(min(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2))) max(max(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2)))];
% t.clims2=[min(min(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2))) max(max(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2)))];
figs{1}=figure;
set(gcf,'Position',[100 100 700 700],'PaperPositionMode','auto')
subplot(2,1,1)
imagesc(TFA.time,TFA.frequency,pl.data_ind,t.clims1)
colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('TFA-amplitude for induced activity\n for channel [%s]', vararg2str(pl.elec2plot)), 'FontSize',8)

xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
hline(14.16667,'m')
hline(pl.fre2index,'c')
set(gca,'FontSize',8)
% set(gca,'ColorScale','log')
cb = colorbar();

% if logscale

subplot(2,1,2)
imagesc(TFA.time,TFA.frequency,pl.data_evo,t.clims2)
colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('TFA-amplitude for evoked activity\n for channel [%s]', vararg2str(pl.elec2plot)), 'FontSize',8)
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
hline(14.16667,'m')
hline(pl.fre2index,'c')
set(gca,'FontSize',8)
% set(gca,'ColorScale','log')
cb = colorbar();
% cb.Ruler.Scale = 'log';
% cb.Ruler.MinorTick = 'on';

% draw topography with electrode positions
h.a1 = axes('position',[0.85 0.45 0.14 0.14],'Visible','off');
topoplot(find(pl.elec2plot_i),TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',5,1});




% normalized induced no RESS
pl.data_ind = squeeze(nanmean(nanmean(TFA.data_induced_bc(:,:,pl.elec2plot_i,pl.subs2use),3),4));
t.clims1=[-1 1]*max(max(abs(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2))));
figs{2}=figure;
set(gcf,'Position',[100 100 700 700],'PaperPositionMode','auto')
h.sp1=subplot(2,1,1);
pl.colmap = [[linspace(0,0.9,128)' linspace(0,0.9,128)' linspace(1,0.9,128)']; [linspace(0.9,1,128)' linspace(0.9,0,128)' linspace(0.9,0,128)']];
% pl.colmap = [[zeros(64,1) linspace(1,0,64)' ones(64,1)];...
%     [linspace(0,1,64)' linspace(0,1,64)' linspace(1,1,64)'];...
%     [linspace(1,1,64)' linspace(1,0,64)' linspace(1,0,64)'];...
%     [ones(64,1) linspace(0,1,64)' zeros(64,1)]];
% colormap(cmap('R2'))
imagesc(TFA.time,TFA.frequency,pl.data_ind,t.clims1)
% colormap(gca,pl.colmap); % cbrewer('div','RdBu',256,'l')
colormap(gca,flipud(cbrewer2('RdBu')))
set(gca,'YDir','normal')
title(sprintf('TFA-amplitude for induced activity; baseline [%1.0f %1.0f]ms\n for channel [%s]',...
    TFA.baseline,vararg2str(pl.elec2plot)),'FontSize',8)
h.cb1=colorbar;
set(h.cb1,'FontSize',8)
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
% freezeColors
% h.cb1=cbfreeze(h.cb1);
hline(14.16667,'m')
hline(pl.fre2index,'c')
set(gca,'FontSize',8)

t.clims1=[0 1]*max(max(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2)));
t.clims2=[0 1]*max(max(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2)));


h.sp2=subplot(2,1,2);
% ttest
t.datain = squeeze(nanmean(TFA.data_induced_bc(:,:,pl.elec2plot_i,pl.subs2use),3));
[tt.h tt.p tt.ci tt.stats]=ttest(permute(t.datain, [3 1 2]));
TFA.pvals_induced_bc=squeeze(tt.p);
tt.stats.tstat = squeeze(tt.stats.tstat);

% do cluster correction?
flag_cluster = 0;
pl.permute_n = 100;
if flag_cluster == 1
    fprintf('running Monte Carlo Simulation with %1.0f repetitions: %7.0f',pl.permute_n,0)
    % find cluster of p-values of original data
    [clust.L, clust.n]=bwlabel(squeeze(tt.h));
    if clust.n>0
        clust.t_sum = []; clust.idx = [];
        for i_cl = 1:clust.n
            clust.t_sum(i_cl)=sum(tt.stats.tstat(clust.L==i_cl));
        end
        % do cluster calculation
        t.dat_all = squeeze(nanmean(TFA.data_induced_bc(:,:,pl.elec2plot_i,pl.subs2use),3));
        t.dat_all(:,:,end+1:end+size(t.dat_all,3))=zeros(size(t.dat_all));
        % create empiric distribution of t-vals
        clust.idx = false(pl.permute_n,size(t.dat_all,3));
        clust.perm_abs_tsum = nan(pl.permute_n,1);
        for i_perm = 1:pl.permute_n
            fprintf('\b\b\b\b\b\b\b%7.0f',i_perm)
            % permute index
            clust.idx(i_perm,randsample(size(t.dat_all,3),size(t.dat_all,3)/2))=true;
            % ttest for permuted data
            [clust.tt.h clust.tt.p clust.tt.ci clust.tt.stats]= ...
                ttest(permute(t.dat_all(:,:,clust.idx(i_perm,:)),[3 1 2]), permute(t.dat_all(:,:,~clust.idx(i_perm,:)),[3 1 2]));
            clust.tt.stats.tstat = squeeze(clust.tt.stats.tstat);
            % cluster identification
            [clust.perm.L, clust.perm.n]=bwlabel(squeeze(clust.tt.h));
            % t-sum value of cluster
            if clust.perm.n>0
                clust.perm.t_sum = nan(clust.perm.n,1);
                for i_cl = 1:clust.perm.n
                    clust.perm.t_sum(i_cl)=sum(clust.tt.stats.tstat(clust.perm.L==i_cl));
                end
                clust.perm_abs_tsum(i_perm) = max(abs(clust.perm.t_sum));
            else
                clust.perm_abs_tsum(i_perm) = 0;
            end
            
        end
        fprintf('...done\n')
    else
        clust.t_sum = 0;
    end
    % check graphically
%     figure; histogram(clust.perm_abs_tsum,50); hold on; vline(abs(clust.t_sum))
    % extract empirical p-value of cluster
    clust.abs_tsum_testdistribution = [clust.perm_abs_tsum; abs(clust.t_sum)'];
    clust.cluster_p = [];
    clust.cluster_p_map = nan(size(t.datain,1), size(t.datain,2));
    for i_clust = 1:numel(clust.t_sum)
        clust.cluster_p(i_clust ) = sum(clust.abs_tsum_testdistribution>=abs(clust.t_sum(i_clust)))./...
            numel(clust.abs_tsum_testdistribution);
        clust.cluster_p_map(clust.L==i_clust) = clust.cluster_p(i_clust );
    end
    clust.cluster_h_map = clust.cluster_p_map<=.05;
end


t.data = abs(log10(TFA.pvals_induced_bc(t.find1:t.find2,t.ind1:t.ind2)));
t.clim = [0 max(t.data(:))];
t.pcriterion = abs(log10(0.05));
if max(t.data(:))<t.pcriterion
    % temp.colormap = repmat([0.5 0.5 0.5],100,1);
    t.colormap = repmat(linspace(1,0.3,1000)',1,3);
else
    t.border = ceil((t.pcriterion / max(t.data(:)))*1000);
    % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
    t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
end
imagesc(TFA.time,TFA.frequency,abs(log10(TFA.pvals_induced_bc)),t.clim)
set(gca,'YDir','normal')
title(sprintf('TFA-p-values for induced activity; baseline [%1.0f %1.0f]ms\n for channel [%s]',...
    TFA.baseline,vararg2str(pl.elec2plot)),'FontSize',8)
% freezeColors
h.cb2=colorbar;
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
t.yticks = get(h.cb2,'YTick');
set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
    'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')),'FontSize',8)
% h.cb2=cbfreeze(h.cb2);
hline(14.16667,'g')
hline(pl.fre2index,'c')
set(gca,'FontSize',8)
colormap(gca,t.colormap)

% plot cluster contour lines
if flag_cluster == 1
    hold on
    contour(TFA.time,TFA.frequency,clust.cluster_h_map,[1,1],'Color','c','LineWidth',2)
end

% draw topography with electrode positions
h.a1 = axes('position',[0.85 0.45 0.14 0.14],'Visible','off');
topoplot(find(pl.elec2plot_i),TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',5,1});



% normalized evoked noRESS
pl.data_evo = squeeze(nanmean(nanmean(TFA.data_evoked_bc(:,:,pl.elec2plot_i,pl.subs2use),3),4));
t.clims2=[-1 1]*max(max(abs(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2))));
figs{3}=figure;
set(gcf,'Position',[100 100 700 700],'PaperPositionMode','auto')
h.sp1=subplot(2,1,1);
imagesc(TFA.time,TFA.frequency,pl.data_evo,t.clims2)
colormap(gca,flipud(cbrewer2('RdBu')))
set(gca,'YDir','normal')
title(sprintf('TFA-amplitude for evoked activity; baseline [%1.0f %1.0f]ms\n for channel [%s]',...
    TFA.baseline,vararg2str(pl.elec2plot)),'FontSize',8)
% freezeColors

h.cb1=colorbar;
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
set(h.cb1,'FontSize',8)
% h.cb1=cbfreeze(h.cb1);
hline(14.16667,'m')
hline(pl.fre2index,'c')
set(gca,'FontSize',8)

h.sp2=subplot(2,1,2);
for i_freq = 1:numel(TFA.frequency)
    [tt.h tt.p tt.ci tt.stats]=ttest(squeeze(nanmean(TFA.data_evoked_bc(i_freq,:,pl.elec2plot_i,pl.subs2use),3))');
    TFA.pvals_induced_bc(i_freq,:)=tt.p;
end
t.data = abs(log10(TFA.pvals_induced_bc(t.find1:t.find2,:)));
t.clim = [0 max(t.data(:))];
t.pcriterion = abs(log10(0.05));
if max(t.data(:))<t.pcriterion
    % temp.colormap = repmat([0.5 0.5 0.5],100,1);
    t.colormap = repmat(linspace(1,0.3,1000)',1,3);
else
    t.border = ceil((t.pcriterion / max(t.data(:)))*1000);
    % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
    t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
end
imagesc(TFA.time,TFA.frequency,abs(log10(TFA.pvals_induced_bc)),t.clim)
set(gca,'YDir','normal')
title(sprintf('TFA-p-values for evoked activity; baseline [%1.0f %1.0f]ms\n for channel [%s]',...
    TFA.baseline,vararg2str(pl.elec2plot)),'FontSize',8)
% freezeColors
h.cb2=colorbar;
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
t.yticks = get(h.cb2,'YTick');
set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
    'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')),'FontSize',8)
% h.cb2=cbfreeze(h.cb2);
hline(14.16667,'g')
hline(pl.fre2index,'c')
set(gca,'FontSize',8)
colormap(gca,t.colormap)

% draw topography with electrode positions
h.a1 = axes('position',[0.85 0.45 0.14 0.14],'Visible','off');
topoplot(find(pl.elec2plot_i),TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',5,1});

sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
% sav.filenames = {'Resp_TFA_Amp_RAW_EvoIndu';'Resp_TFA_Amp_BC_Indu';'Resp_TFA_Amp_BC_Evo'};
sav.filenames = {'Resp_TFA_Amp_RAW_EvoIndu';'Resp_TFA_Amp_BC_Indu'};
for i_fig = 1:2
    print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'-djpeg','-r300')
    saveas(figs{i_fig},fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'fig')
    print(figs{i_fig},fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'-depsc2', '-vector','-r300')
end



%% plot only single fields

% plotting parameters
% pl.elec2plot = {'Oz'};
% pl.elec2plot = {'PO3';'POz';'PO4';'O1';'Oz';'O2';'Iz'}; % steady state I
% pl.elec2plot = {'POz';'O1';'Oz';'O2';'I1';'Iz';'I2'}; % steady state II
% pl.elec2plot = {'PO4';'O2';'PO8';'P8';'P10';'I2'}; % vis alpha
pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; % vis alpha II
% pl.elec2plot = {'POz'}; % vis alpha II
% pl.elec2plot = {'POz'}; % vis alpha IV
% pl.elec2plot = {'C3';'CP3'}; % motor alpha/beta
% pl.elec2plot = {'C3';'CP3';'C5';'CP5'}; % motor alpha/beta II
% pl.elec2plot = {'F1';'F2';'Fz'}; % frontal
% pl.elec2plot = {'Fp1';'Fp2'}; % frontal
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

% pl.xlims=[-2500 2500]; % index time 2 plot
pl.xlims=[-3500 3500]; % index time 2 plot
[t.t t.ind1]=min(abs(TFA.time-pl.xlims(1)));
[t.t t.ind2]=min(abs(TFA.time-pl.xlims(2)));

pl.flims = TFA.frequency([1 end]); % index frequency 2 plot
pl.flims = [3 35];
% pl.flims = [8 44];
[t.t t.find1]=min(abs(TFA.frequency-pl.flims(1)));
[t.t t.find2]=min(abs(TFA.frequency-pl.flims(2)));

pl.subs2use = 1:numel(F.Subjects2Use);
% pl.subs2use = [1:14 16 :20];
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 13 14 17 18 19 20]; % subjects discarded due to low trial number version 1
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20]; % subjects discarded due to low trial number version 2



%%%%%%%%%%%%%
% raw induced no RESS
pl.data_ind = squeeze(nanmean(nanmean(TFA.data_induced(:,:,pl.elec2plot_i,pl.subs2use),3),4));
pl.data_evo = squeeze(nanmean(nanmean(TFA.data_evoked(:,:,pl.elec2plot_i,pl.subs2use),3),4));
t.clims1=[0 1]*max(max(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2)));
t.clims2=[0 1]*max(max(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2)));
figure;
imagesc(TFA.time,TFA.frequency,pl.data_ind,t.clims1)
colormap(gca, fake_parula)
set(gca,'YDir','normal')
title(sprintf('tfa noRESS, for channel %s;\n averaged in frequency domain', vararg2str(pl.elec2plot)), 'FontSize',8)
colorbar
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
hline(14.16667,'m')
set(gca,'FontSize',8)


% raw normalized no RESS
pl.data_ind = squeeze(nanmean(nanmean(TFA.data_induced_bc(:,:,pl.elec2plot_i,pl.subs2use),3),4));
t.clims1=[-1 1]*max(max(abs(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2))));
figure;
colormap(gca,flipud(cbrewer2('RdBu')))
% colormap(gca,cmap('R1'))
% colormap(gca,cmap('D9'))
imagesc(TFA.time,TFA.frequency,pl.data_ind,t.clims1)
set(gca,'YDir','normal')
title(sprintf('normalized tfa noRESS, for channel %s;\n averaged in frequency domain; baseline [%1.0f %1.0f]ms',...
    vararg2str(pl.elec2plot),TFA.baseline),'FontSize',8)
h.cb1=colorbar;
set(h.cb1,'FontSize',8)
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
% freezeColors
% h.cb1=cbfreeze(h.cb1);
hline(14.16667,'m')
set(gca,'FontSize',8)

for i_freq = 1:numel(TFA.frequency)
    [tt.h tt.p tt.ci tt.stats]=ttest(squeeze(nanmean(TFA.data_induced_bc(i_freq,:,pl.elec2plot_i,pl.subs2use),3))');
    TFA.pvals_induced_bc(i_freq,:)=tt.p;
end
pl.alpha = ones(size(TFA.pvals_induced_bc)).*0.85;
pl.alpha(TFA.pvals_induced_bc<.05) = 1;
% alpha(pl.alpha)


%% statistically testing differences in TFA images between signals from different sources
% figure; topoplot([],TFA.electrodes, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);

% plotting parameters
pl.parameters = [];
% vis alpha II large
pl.parameters(1).elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'};
pl.parameters(1).name = 'VisualLarge';
% motor alpha/beta II
pl.parameters(end+1).elec2plot = {'C3';'CP3';'C5';'CP5'};
pl.parameters(end).name = 'MotorLarge';
% parietal
pl.parameters(end+1).elec2plot = {'CP1';'CPz';'CP2';'P1';'Pz';'P2'};
pl.parameters(end).name = 'Parietal';

for i_pl = 1:numel(pl.parameters)
    pl.parameters(i_pl).elec2plot_i = logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.parameters(i_pl).elec2plot, 'UniformOutput',false)),1));
end

% pl.xlims=[-2500 2500]; % index time 2 plot
pl.xlims=[-3250 3250]; % index time 2 plot
[t.t t.ind1]=min(abs(TFA.time-pl.xlims(1)));
[t.t t.ind2]=min(abs(TFA.time-pl.xlims(2)));

% pl.flims = TFA.frequency([1 end]); % index frequency 2 plot
% pl.flims = [3 TFA.frequency(end)]; % index frequency 2 plot
pl.flims = [4 TFA.frequency(end)]; % index frequency 2 plot
% pl.flims = [4 35];
% pl.flims = [3 35];
% pl.flims = [8 44];
[t.t t.find1]=min(abs(TFA.frequency-pl.flims(1)));
[t.t t.find2]=min(abs(TFA.frequency-pl.flims(2)));

pl.baseline = TFA.baseline; %[-3250 -3000];
pl.baseline_i = dsearchn(TFA.time', pl.baseline');

pl.subs2use = 1:numel(F.Subjects2Use);
clear figs fig


%extract data first
pl.testdata = nan(numel(t.find1:t.find2),numel(t.ind1:t.ind2),numel(pl.parameters),numel(F.Subjects2Use),2);
% freq, time, electrode cluster, subjects, [baseline time]
for i_pl = 1:numel(pl.parameters)
    pl.testdata(:,:,i_pl,:,1)= mean(TFA.data_induced(t.find1:t.find2,t.ind1:t.ind2,pl.parameters(i_pl).elec2plot_i,pl.subs2use),3);
    pl.testdata(:,:,i_pl,:,2)= repmat(...
       mean(mean(TFA.data_induced(t.find1:t.find2,pl.baseline_i(1):pl.baseline_i(2),pl.parameters(i_pl).elec2plot_i,pl.subs2use),2),3),...
       [1,size(pl.testdata,2),1,1]);
end

% do liner mixed modelling for each frequency/time combination
prog = 0; i_up = 0;TFA.lme=[];
fprintf('doing %1.0f lme_fits; progress: %5.2f%%\n',size(pl.testdata,1)*size(pl.testdata,2),prog)
for i_freq = 1:size(pl.testdata,1)
    for i_time = 1:size(pl.testdata,2)
        prog = ( 100*(i_up/(size(pl.testdata,1)*size(pl.testdata,2))));
        fprintf(1,'\b\b\b\b\b\b%5.2f%%',prog);
        % data in?
        t.data = squeeze(pl.testdata(i_freq,i_time,:,:,:));
        t.fct1 = repmat({pl.parameters(:).name}',[1 size(t.data,2) size(t.data,3)]);
        t.fct2 = repmat(pl.subs2use,[size(t.data,1) 1 size(t.data,3)]);
        t.fct3 = permute(repmat([1;2], [1 size(t.data,1) size(t.data,2)]),[2 3 1]);
        
        t.tbldata = table(t.data(:),t.fct1(:),t.fct2(:),t.fct3(:),...
            'VariableNames',{'amplitude','source','participant','time'});
        
        % fit lme
        t.lme = fitlme(t.tbldata,'amplitude~source*time+(1|participant)');
        TFA.lme.F(i_freq,i_time).source = t.lme.anova.FStat(2);
        TFA.lme.p(i_freq,i_time).source = t.lme.anova.pValue(2);
        TFA.lme.F(i_freq,i_time).time = t.lme.anova.FStat(3);
        TFA.lme.p(i_freq,i_time).time = t.lme.anova.pValue(3);
        TFA.lme.F(i_freq,i_time).sourceXtime = t.lme.anova.FStat(3);
        TFA.lme.p(i_freq,i_time).sourceXtime = t.lme.anova.pValue(3);
        %         % post-hoc coefficient tests
        %         [t.pVal,t.F,t.DF1,t.DF2] = coefTest(t.lme,[0 0 0 1 0 0])
        i_up = i_up+1;
    end
end


%%%%%%%%%%%%%
% raw induced no RESS
pl.data_ind = squeeze(nanmean(nanmean(TFA.data_induced(:,:,pl.elec2plot_i,pl.subs2use),3),4));
pl.data_evo = squeeze(nanmean(nanmean(TFA.data_evoked(:,:,pl.elec2plot_i,pl.subs2use),3),4));
t.clims1=[0 1]*max(max(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2)));
t.clims2=[0 1]*max(max(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2)));
% t.clims1=[min(min(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2))) max(max(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2)))];
% t.clims2=[min(min(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2))) max(max(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2)))];
figs{1}=figure;
set(gcf,'Position',[100 100 700 700],'PaperPositionMode','auto')
subplot(2,1,1)
imagesc(TFA.time,TFA.frequency,pl.data_ind,t.clims1)
colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('TFA-amplitude for induced activity\n for channel [%s]', vararg2str(pl.elec2plot)), 'FontSize',8)


%% plot RESS components
pl.xlims=[-3250 3250]; % index time 2 plot
% pl.xlims=[-3500 3500]; % index time 2 plot
[t.t t.ind1]=min(abs(TFA.time-pl.xlims(1)));
[t.t t.ind2]=min(abs(TFA.time-pl.xlims(2)));

% pl.flims = TFA.frequency([1 end]); % index frequency 2 plot
% pl.flims = [4 44];
% pl.flims = [4 35];
pl.flims = [3 35];
[t.t t.find1]=min(abs(TFA.frequency-pl.flims(1)));
[t.t t.find2]=min(abs(TFA.frequency-pl.flims(2)));

pl.subs2use = 1:numel(F.Subjects2Use);
% pl.subs2use = [1:14 16 :20];
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 13 14 17 18 19 20]; % subjects discarded due to low trial number version 1
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20]; % subjects discarded due to low trial number version 2
clear figs

% raw induced RESS
pl.data_ind = squeeze(nanmean(TFA.data_RESS_induced(:,:,pl.subs2use),3));
pl.data_evo = squeeze(nanmean(TFA.data_RESS_evoked(:,:,pl.subs2use),3));
t.clims1=[0 1]*max(max(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2)));
t.clims2=[0 1]*max(max(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2)));
figs{1}=figure;
set(gcf,'Position',[100 100 700 700],'PaperPositionMode','auto')

subplot(2,1,1)
imagesc(TFA.time,TFA.frequency,pl.data_ind,t.clims1)
colormap(gca, fake_parula)
set(gca,'YDir','normal')
title(sprintf('tfa RESS component; averaged in frequency domain\nRESS'),'FontSize',8)
colorbar
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
hline(14.16667,'m')
set(gca,'FontSize',8)

subplot(2,1,2)
imagesc(TFA.time,TFA.frequency,pl.data_evo,t.clims2)
colormap(gca, fake_parula)
set(gca,'YDir','normal')
title(sprintf('tfa RESS component; averaged in time domain\nRESS'),'FontSize',8)
colorbar
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
hline(14.16667,'m')
set(gca,'FontSize',8)



% normalized induced RESS
pl.data_ind = squeeze(nanmean(TFA.data_RESS_induced_bc(:,:,pl.subs2use),3));
t.clims1=[-1 1]*max(max(abs(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2))));
figs{2}=figure;
set(gcf,'Position',[100 100 700 700],'PaperPositionMode','auto')

h.sp1=subplot(2,1,1);

% colormap(cmap('R2'))
imagesc(TFA.time,TFA.frequency,pl.data_ind,t.clims1)
colormap(gca,flipud(cbrewer2('RdBu')))
set(gca,'YDir','normal')
title(sprintf('normalized tfa RESS component;\n averaged in frequency domain; baseline [%1.0f %1.0f]ms', TFA.baseline),'FontSize',8)
h.cb1=colorbar;
set(h.cb1,'FontSize',8)
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
% freezeColors
% h.cb1=cbfreeze(h.cb1);
hline(14.16667,'m')
set(gca,'FontSize',8)

h.sp2=subplot(2,1,2);
for i_freq = 1:numel(TFA.frequency)
    [tt.h tt.p tt.ci tt.stats]=ttest(squeeze(TFA.data_RESS_induced_bc(i_freq,:,pl.subs2use))');
    TFA.pvals_RESS_induced_bc(i_freq,:)=tt.p;
end
t.data = abs(log10(TFA.pvals_RESS_induced_bc(t.find1:t.find2,:)));
t.clim = [0 max(t.data(:))];
t.pcriterion = abs(log10(0.05));
if max(t.data(:))<t.pcriterion
    % temp.colormap = repmat([0.5 0.5 0.5],100,1);
    t.colormap = repmat(linspace(1,0.3,1000)',1,3);
else
    t.border = ceil((t.pcriterion / max(t.data(:)))*1000);
    % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
    t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
end
imagesc(TFA.time,TFA.frequency,t.data,t.clim)
set(gca,'YDir','normal')
title(sprintf('p-values RESS component;\n averaged in frequency domain; baseline [%1.0f %1.0f]ms', TFA.baseline),'FontSize',8)
% freezeColors
h.cb2=colorbar;
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
t.yticks = get(h.cb2,'YTick');
set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
    'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')),'FontSize',8)
% h.cb2=cbfreeze(h.cb2);
hline(14.16667,'g')
set(gca,'FontSize',8)
colormap(gca,t.colormap)

% normalized evoked RESS
pl.data_evo = squeeze(nanmean(TFA.data_RESS_evoked_bc(:,:,pl.subs2use),3));
t.clims2=[-1 1]*max(max(abs(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2))));
figs{3}=figure;
set(gcf,'Position',[100 100 700 700],'PaperPositionMode','auto')

h.sp1=subplot(2,1,1);
% colormap(cmap('R2'))
imagesc(TFA.time,TFA.frequency,pl.data_evo,t.clims2)
set(gca,'YDir','normal')
colormap(gca,flipud(cbrewer2('RdBu')))
title(sprintf('normalized tfa RESS component;\n averaged in time domain; baseline [%1.0f %1.0f]ms', TFA.baseline),'FontSize',8)
% freezeColors

h.cb1=colorbar;
set(h.cb1,'FontSize',8)
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)

% h.cb1=cbfreeze(h.cb1);
hline(14.16667,'m')
set(gca,'FontSize',8)

h.sp2=subplot(2,1,2);
for i_freq = 1:numel(TFA.frequency)
    [tt.h tt.p tt.ci tt.stats]=ttest(squeeze(TFA.data_RESS_evoked_bc(i_freq,:,pl.subs2use))');
    TFA.pvals_RESS_evoked_bc(i_freq,:)=tt.p;
end
t.data = abs(log10(TFA.pvals_RESS_evoked_bc(t.find1:t.find2,:)));
t.clim = [0 max(t.data(:))];
t.pcriterion = abs(log10(0.05));
if max(t.data(:))<t.pcriterion
    % temp.colormap = repmat([0.5 0.5 0.5],100,1);
    t.colormap = repmat(linspace(1,0.3,1000)',1,3);
else
    t.border = ceil((t.pcriterion / max(t.data(:)))*1000);
    % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
    t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
end
imagesc(TFA.time,TFA.frequency,t.data,t.clim)
set(gca,'YDir','normal')
title(sprintf('p-values RESS component;\n averaged in time domain; baseline [%1.0f %1.0f]ms', TFA.baseline),'FontSize',8)
% freezeColors
h.cb2=colorbar;
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
t.yticks = get(h.cb2,'YTick');
set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
    'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')),'FontSize',8)
% h.cb2=cbfreeze(h.cb2);
hline(14.16667,'g')
set(gca,'FontSize',8)
colormap(gca,t.colormap)

% sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
% sav.filenames = {'Resp_TFA_Amp_RAW_EvoIndu_RESS';'Resp_TFA_Amp_BC_Indu_RESS';'Resp_TFA_Amp_BC_Evo_RESS'};
% for i_fig = 1:3
%     print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'-djpeg','-r300')
%     saveas(figs{i_fig},fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'fig')
%     print(figs{i_fig},fullfile(sav.pathout,sprintf('%s_%s,',sav.filenames{i_fig},sav.chan_add)),'-depsc2', '-painters','-r300')
% end

%% plot timecourse of frequency
% pl.freq2plot=[14.16667 14.16667];
pl.freq2plot=[8 14];
% pl.freq2plot=[15 25];
% pl.freq2plot=[20 30];
% pl.freq2plot=[8 14];
% pl.freq2plot=[15 30];
% pl.freq2plot=[2 4];
[t.t t.ind3]=min(abs(TFA.frequency-pl.freq2plot(1)));
[t.t t.ind4]=min(abs(TFA.frequency-pl.freq2plot(end)));

% pl.elec2plot = {'Oz';'POz'};
% pl.elec2plot = {'PO3';'POz';'PO4';'O1';'Oz';'O2';'Iz'}; % steady state I
% pl.elec2plot = {'POz';'O1';'Oz';'O2';'I1';'Iz';'I2'}; % steady state II
% pl.elec2plot = {'PO4';'O2';'PO8';'P8';'P10';'I2'}; % vis alpha I
% pl.elec2plot = {'POz'}; % vis alpha II
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; % vis alpha III
% pl.elec2plot = {'POz';'Pz';'P1';'P2'}; % vis alpha II
pl.elec2plot = {'C3';'CP3'}; % motor alpha/beta
% pl.elec2plot = {'C3';'CP3';'C5';'CP5'}; % motor alpha/beta II
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.xlims=[-2750 3000]; % index time 2 plot
[t.t t.ind1]=min(abs(TFA.time-pl.xlims(1)));
[t.t t.ind2]=min(abs(TFA.time-pl.xlims(2)));

pl.subs2use = 1:numel(F.Subjects2Use);

pl.data_ind = squeeze(nanmean(nanmean(TFA.data_induced_bc(:,:,pl.elec2plot_i,pl.subs2use),3),4));
pl.data_evo = squeeze(nanmean(nanmean(TFA.data_evoked_bc(:,:,pl.elec2plot_i,pl.subs2use),3),4));

% ttest for each point
[tt.h pl.pdata_ind tt.ci tt.stats]=ttest(squeeze(nanmean(nanmean(TFA.data_induced_bc(t.ind3:t.ind4,:,pl.elec2plot_i,:),1),3))');
[tt.h pl.pdata_evo tt.ci tt.stats]=ttest(squeeze(nanmean(nanmean(TFA.data_evoked_bc(t.ind3:t.ind4,:,pl.elec2plot_i,:),1),3))');

figure;
subplot(2,1,1)
plot(TFA.time,mean(pl.data_ind(t.ind3:t.ind4,:),1),'k')
hold on
try
   pl.tdata=mean(pl.data_ind(t.ind3:t.ind4,:),1);
   pl.tdata(pl.pdata_ind>=.05)=nan;
   plot(TFA.time,pl.tdata,'b','LineWidth',2)
   try
       pl.tdata(pl.pdata_ind>=.01)=nan;
       plot(TFA.time,pl.tdata,'c','LineWidth',2)
       try
           pl.tdata(pl.pdata_ind>=.001)=nan;
           plot(TFA.time,pl.tdata,'r','LineWidth',2)
           legend({'amplitude';'p<.05';'p<.01';'p<.001'})
       catch
           legend({'amplitude';'p<.05';'p<.01'})
       end
   catch
       legend({'amplitude';'p<.05'})
   end
catch 
    legend('amplitude')
end
ylim([-1.2 1.2]*max(abs(mean(pl.data_ind(t.ind3:t.ind4,t.ind1:t.ind2),1))))
xlabel('time in ms')
ylabel('amplitude in \muV')
xlim(pl.xlims)
title(sprintf('induced timecourse signal %1.2f to %1.2f Hz\n channel %s; baseline [%1.0f %1.0f]ms',...
    pl.freq2plot, vararg2str(pl.elec2plot),TFA.baseline))
grid on
hline(0,'k')

subplot(2,1,2)
plot(TFA.time,mean(pl.data_evo(t.ind3:t.ind4,:),1),'k')
hold on
try
   pl.tdata=mean(pl.data_evo(t.ind3:t.ind4,:),1);
   pl.tdata(pl.pdata_evo>=.05)=nan;
   plot(TFA.time,pl.tdata,'b','LineWidth',2)
   try
       pl.tdata(pl.pdata_evo>=.01)=nan;
       plot(TFA.time,pl.tdata,'c','LineWidth',2)
       try
           pl.tdata(pl.pdata_evo>=.001)=nan;
           plot(TFA.time,pl.tdata,'r','LineWidth',2)
           legend({'amplitude';'p<.05';'p<.01';'p<.001'})
       catch
           legend({'amplitude';'p<.05';'p<.01'})
       end
   catch
       legend({'amplitude';'p<.05'})
   end
catch 
    legend('amplitude')
end
ylim([-1.2 1.2]*max(abs(mean(pl.data_evo(t.ind3:t.ind4,t.ind1:t.ind2),1))))
xlabel('time in ms')
ylabel('amplitude in \muV')
xlim(pl.xlims)
title(sprintf('evoked timecourse signal %1.2f to %1.2f Hz\n channel %s; baseline [%1.0f %1.0f]ms',...
    pl.freq2plot, vararg2str(pl.elec2plot),TFA.baseline))
grid on
hline(0,'k')

%% plot all lines of interest into one graphics

% pl.parameters = {...
%     [14.16667 14.16667] {'RESS'} 'data_RESS_induced_bc' [0.6 0.6 0.6] 2 sprintf('vis ind\nSSVEP') 'RESS';...
%     [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('vis evo\nSSVEP') 'RESS';...
%     [8 12] {'PO4';'O2';'PO8';'P8';'P10';'I2'} 'data_induced_bc' [1 0 0] 1 'vis alpha' 'noRESS';...
%     [18 30] {'PO4';'O2';'PO8';'P8';'P10';'I2'} 'data_induced_bc' [1 0 1] 1 'vis beta' 'noRESS';...
%     [10 14] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0] 1 'mot alpha' 'noRESS';...
%     [18 30] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0.9] 1 'mot beta' 'noRESS';...
%     };
pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'data_RESS_induced_bc' [0.6 0.6 0.6] 2 sprintf('vis ind\nSSVEP') 'RESS';...
    [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('vis evo\nSSVEP') 'RESS';...
    [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [1 0 0] 1 'vis alpha' 'noRESS';...
    [18 30] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [1 0 1] 1 'vis beta' 'noRESS';...
    [10 14] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0] 1 'mot alpha' 'noRESS';...
    [18 30] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0.9] 1 'mot beta' 'noRESS';...
    };
% pl.parameters = {...
%     [14.16667 14.16667] {'RESS'} 'data_RESS_induced_bc' [0.6 0.6 0.6] 2 sprintf('vis ind\nSSVEP') 'RESS';...
%     [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('vis evo\nSSVEP') 'RESS';...
%     [8 12] {'POz'} 'data_induced_bc' [1 0 0] 1 'vis alpha' 'noRESS';...
%     [18 30] {'POz'} 'data_induced_bc' [1 0 1] 1 'vis beta' 'noRESS';...
%     [10 14] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0] 1 'mot alpha' 'noRESS';...
%     [18 30] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0.9] 1 'mot beta' 'noRESS';...
%     };

pl.subs2use = 1:numel(F.Subjects2Use);
% pl.subs2use = [1:14 16 :20];
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 13 14 17 18 19 20]; % subjects discarded due to low trial number version 1
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20]; % subjects discarded due to low trial number version 2
% 
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; % subjects discarded due to low trial number version 2


% pl.xlim = [-3250 3250];
pl.xlim = [-4000 4000];
pl.xlim = [-3500 3500];

pl.toponum = [ceil(sqrt(size(pl.parameters,1))) round(sqrt(size(pl.parameters,1)))];
pl.plpos = repmat([1:pl.toponum(2)*3],pl.toponum(1),1)+repmat((0:pl.toponum(2)*4:pl.toponum(2)*4*(pl.toponum(1)-1))',1,pl.toponum(2)*3);
pl.topopos = repmat([pl.toponum(2)*3+1:pl.toponum(2)*4],pl.toponum(1),1)'+repmat((0:pl.toponum(2)*4:pl.toponum(2)*4*(pl.toponum(1)-1))',1,pl.toponum(2))';

h.pl = [];
figure;
subplot(pl.toponum(1),pl.toponum(2)*4,pl.plpos(:))
for i_pl = 1:size(pl.parameters,1)
    % index frequencies
    [t.t t.ind1]=min(abs(TFA.frequency-pl.parameters{i_pl,1}(1)));
    [t.t t.ind2]=min(abs(TFA.frequency-pl.parameters{i_pl,1}(2)));
    
    % index time
    [t.t t.ind3]=min(abs(TFA.time-pl.xlim(1)));
    [t.t t.ind4]=min(abs(TFA.time-pl.xlim(2)));
    
    % index electrodes
    switch  pl.parameters{i_pl,7}
        case 'noRESS'
            pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.parameters{i_pl,2}, 'UniformOutput',false)),1));
            com=sprintf('pl.data=squeeze(nanmean(nanmean(TFA.%s(t.ind1:t.ind2,t.ind3:t.ind4,pl.elec2plot_i,pl.subs2use),1),3));',pl.parameters{i_pl,3});
        case 'RESS'
            pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes_RESS.labels},x), pl.parameters{i_pl,2}, 'UniformOutput',false)),1));
            com=sprintf('pl.data=squeeze(nanmean(TFA.%s(t.ind1:t.ind2,t.ind3:t.ind4,pl.subs2use),1));',pl.parameters{i_pl,3});
        
    end
    eval(com)
    [tt.h tt.p tt.ci tt.stats]=ttest(pl.data');
       
    % actual plot
    h.pl(i_pl)=plot(TFA.time(t.ind3:t.ind4),mean(pl.data,2),'Color',pl.parameters{i_pl,4},'LineWidth',pl.parameters{i_pl,5});
    hold on
    
    % add significant portions
    pl.tdata = nan(size(pl.data));
    pl.t2data = mean(pl.data,2);
    try
        pl.tdata(tt.p<.05)=pl.t2data(tt.p<.05);
    end
    plot(TFA.time(t.ind3:t.ind4),pl.tdata,'Color',pl.parameters{i_pl,4},'LineWidth',pl.parameters{i_pl,5}+2)
end
xlim(pl.xlim)
set(gca,'ylim',[-1.1 1.1].*max(cellfun(@(x) max(abs(x)), get(get(gca,'Children'),'YData'))))
legend(h.pl, pl.parameters(:,6),'location','NorthWest','FontSize', 8)
grid on
hline(0,'k')
xlabel('time in ms')
ylabel('amplitude')

for i_pl = 1:size(pl.parameters,1)
    subplot(pl.toponum(1),pl.toponum(2)*4,pl.topopos(i_pl))
    pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.parameters{i_pl,2}, 'UniformOutput',false)),1));
    topoplot([],TFA.electrodes(1:64),'style','blank','electrodes', 'on', 'emarker2', {find(pl.elec2plot_i),'o','r',2});
    title(pl.parameters{i_pl,6},'Color',pl.parameters{i_pl,4})
end

%% plot all lines of interest into one graphics (bounded line)
clear figs h


pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
    [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [1 0 0] 2 'vis alpha' 'noRESS';...
    [10 14] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [0 0.5 0] 2 'mot alpha' 'noRESS';...
    [18 30] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [0 0.5 0.9] 2 'mot beta' 'noRESS';...
    [15 23] {'CP1';'CPz';'CP2';'P1';'Pz';'P2'} 'data_induced_bc' [0.8 0 0.8] 2 'parietal beta' 'noRESS';...
    };

% this one for the poster figures
pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
    [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [255 151 0]./255 2 'vis alpha' 'noRESS';...
    [10 14] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [33 92 150]./255 2 'mot alpha' 'noRESS';...
    [18 30] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [170 99 57]./255 2 'mot beta' 'noRESS';...
    [15 23] {'CP1';'CPz';'CP2';'P1';'Pz';'P2'} 'data_induced_bc' [23 150 118]./255 2 'parietal beta' 'noRESS';...
    };


pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
    [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [255 151 0]./255 2 'vis alpha' 'noRESS';...
    [10 14] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [33 92 150]./255 2 'mot alpha' 'noRESS';...
    [15 30] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [170 99 57]./255 2 'mot beta' 'noRESS';...
    [15 23] {'CP1';'CPz';'CP2';'P1';'Pz';'P2'} 'data_induced_bc' [23 150 118]./255 2 'parietal beta' 'noRESS';...
    };

% pl.parameters = {...
%     [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
%     [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [1 0 0] 2 'vis alpha' 'noRESS';...
%     [15 21] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [0.4 0.4 0.4] 2 'vis beta' 'noRESS';...
%     [10 14] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [0 0.5 0] 2 'mot alpha' 'noRESS';...
%     [18 30] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [0 0.5 0.9] 2 'mot beta' 'noRESS';...
%     [16 23] {'CP1';'CPz';'CP2';'P1';'Pz';'P2'} 'data_induced_bc' [0.8 0 0.8] 2 'parietal beta' 'noRESS';...
%     [8 12] {'CP1';'CPz';'CP2';'P1';'Pz';'P2'} 'data_induced_bc' [0 0.8 0.8] 2 'parietal alpha' 'noRESS';...
%     };

% most straight forward?
pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
    [8 14] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [255 151 0]./255 2 'vis alpha' 'noRESS';...
    [8 14] {'C3';'CP3'} 'data_induced_bc' [33 92 150]./255 2 'mot alpha' 'noRESS';...
    [15 20] {'C3';'CP3'} 'data_induced_bc' [170 99 57]./255 2 'mot low beta' 'noRESS';...
    [20 30] {'C3';'CP3'} 'data_induced_bc' [23 150 118]./255 2 'mot high beta' 'noRESS';...
    };




pl.subs2use = 1:numel(F.Subjects2Use);
% pl.subs2use = [1:14 16 :20];
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 13 14 17 18 19 20]; % subjects discarded due to low trial number version 1
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20]; % subjects discarded due to low trial number version 2
% 
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; % subjects discarded due to low trial number version 2

% pl.xlim = [-3250 3250];
pl.xlim = [-4000 4000];
pl.xlim = [-3500 3500];

pl.toponum = [ceil(sqrt(size(pl.parameters,1))) round(sqrt(size(pl.parameters,1)))];
pl.plpos = repmat([1:pl.toponum(2)*3],pl.toponum(1),1)+repmat((0:pl.toponum(2)*4:pl.toponum(2)*4*(pl.toponum(1)-1))',1,pl.toponum(2)*3);
pl.topopos = repmat([pl.toponum(2)*3+1:pl.toponum(2)*4],pl.toponum(1),1)'+repmat((0:pl.toponum(2)*4:pl.toponum(2)*4*(pl.toponum(1)-1))',1,pl.toponum(2))';
pl.permut_n = 10000;
% pl.permut_n = 100;


h.pl = [];
figs{1}=figure;
set(gcf,'Position',[100 100 800 350],'PaperPositionMode','auto')
h.sp(1)=subplot(pl.toponum(1),pl.toponum(2)*4,pl.plpos(:));
clear pl.tdata2
for i_pl = 1:size(pl.parameters,1)
    % index frequencies
    [t.t t.ind1]=min(abs(TFA.frequency-pl.parameters{i_pl,1}(1)));
    [t.t t.ind2]=min(abs(TFA.frequency-pl.parameters{i_pl,1}(2)));
    
    % index time
    [t.t t.ind3]=min(abs(TFA.time-pl.xlim(1)));
    [t.t t.ind4]=min(abs(TFA.time-pl.xlim(2)));
    
    % index electrodes
    switch  pl.parameters{i_pl,7}
        case 'noRESS'
            pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.parameters{i_pl,2}, 'UniformOutput',false)),1));
            com=sprintf('pl.data=squeeze(nanmean(nanmean(TFA.%s(t.ind1:t.ind2,t.ind3:t.ind4,pl.elec2plot_i,pl.subs2use),1),3));',pl.parameters{i_pl,3});
        case 'RESS'
            pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes_RESS.labels},x), pl.parameters{i_pl,2}, 'UniformOutput',false)),1));
            com=sprintf('pl.data=squeeze(nanmean(TFA.%s(t.ind1:t.ind2,t.ind3:t.ind4,pl.subs2use),1));',pl.parameters{i_pl,3});
        
    end
    eval(com)
    [tt.h tt.p tt.ci tt.stats]=ttest(pl.data');
    % do cluster based permutation
    [cluster_runt_diff, timecourse_runt_diff]=eeg_erpStat_clusterP(pl.data,zeros(size(pl.data)),pl.permut_n,2);
    tt.p(timecourse_runt_diff.h_corr==0)=1;
       
    % actual plot with boundedline
    pl.mdata = mean(pl.data,2)';
    pl.ddata = tt.ci - repmat(mean(pl.data,2)',2,1);
    %pl.ddata = tt.ci;
    [h.l(i_pl), h.p(i_pl)] = boundedline(TFA.time(t.ind3:t.ind4),pl.mdata,pl.ddata(1,:));
    set(h.l(i_pl),'Color',pl.parameters{i_pl,4},'LineWidth',pl.parameters{i_pl,5});
    set(h.p(i_pl),'FaceColor',pl.parameters{i_pl,4},'FaceAlpha',0.25)
    hold on
    
%     %%%%%%% vers1 %start%
%     % add significant portions
%     pl.tdata = nan(size(pl.mdata));
%     pl.t2data = mean(pl.data,2);
%     try
%         pl.tdata(tt.p<.05)=pl.t2data(tt.p<.05);
%     end
%     plot(TFA.time(t.ind3:t.ind4),pl.tdata,'Color',[pl.parameters{i_pl,4} 1],'LineWidth',pl.parameters{i_pl,5}+3)
%     %%%%%%% vers1 %end%
    
    pl.tdata2(i_pl,:)=nan(size(pl.mdata));
    pl.tdata2(i_pl,tt.p<.05)=1; 
end
xlim(pl.xlim)
% %%%%%%% vers1 %start%
% set(gca,'ylim',[-1.1 1.1].*max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false))))
% %%%%%%% vers1 %end%

%%%%%%% vers2 %start%
% plot lines for significant effects below
set(gca,'ylim',[-(1+0.15+0.05*size(pl.parameters,1)) 1.1].*max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false))))
pl.maxdata = max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false)));
for i_pl = 1:size(pl.parameters,1)
    pl.y = -1*pl.maxdata*(1+0.05+0.05*i_pl);
    
%     plot(TFA.time(t.ind3:t.ind4),repmat(pl.y,1,numel(TFA.time(t.ind3:t.ind4))),'Color','k','LineWidth',0.1)
    plot(TFA.time(t.ind3:t.ind4),pl.tdata2(i_pl,:).*pl.y,'Color',[pl.parameters{i_pl,4} 1],'LineWidth',2)
    
end
%%%%%%% vers2 %end%

legend(h.l, pl.parameters(:,6),'location','NorthWest','FontSize', 8)
grid on
box on
hline(0,'k')
xlabel('time in ms')
ylabel('amplitude modulation in %')

for i_pl = 1:size(pl.parameters,1)
    subplot(pl.toponum(1),pl.toponum(2)*4,pl.topopos(i_pl))
    pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.parameters{i_pl,2}, 'UniformOutput',false)),1));
%     topoplot([],TFA.electrodes(1:64),'style','blank','electrodes', 'on', 'emarker2', {find(pl.elec2plot_i),'o','r',2});
    topoplot(find(pl.elec2plot_i),TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
        'emarker2',{find(pl.elec2plot_i),'o',pl.parameters{i_pl,4},3,1});
    title(sprintf('%s\n[%1.1f %1.1f]Hz',pl.parameters{i_pl,6},pl.parameters{i_pl,1}),'Color',pl.parameters{i_pl,4})
end

sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
% sav.filenames = {'Resp_AllSignals_Amp_Timecourse'};
% sav.filenames = {'Resp_AllSignals_Amp_Timecourse_v2'};
% sav.filenames = {'Resp_AllSignals_Amp_Timecourse_5'};
sav.filenames = {'Resp_AllSignals_Amp_Timecourse_v2_clustcorr_5'};
% sav.filenames = {'Resp_AllSignals_Amp_Timecourse_v2_clustcorr_5b'};
for i_fig = 1:1
%     print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-dpng','-r300')
%     print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-djpeg','-r300')
%     saveas(figs{i_fig},fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'fig')
%     print(figs{i_fig},fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-depsc2', '-painters','-r300')
end

%% plot all lines of interest into one graphics (separate graphics + single subject data)
clear figs h


pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
    [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [1 0 0] 2 'vis alpha' 'noRESS';...
    [10 14] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [0 0.5 0] 2 'mot alpha' 'noRESS';...
    [18 30] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [0 0.5 0.9] 2 'mot beta' 'noRESS';...
    [15 23] {'CP1';'CPz';'CP2';'P1';'Pz';'P2'} 'data_induced_bc' [0.8 0 0.8] 2 'parietal beta' 'noRESS';...
    };

pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
    [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [255 151 0]./255 2 'vis alpha' 'noRESS';...
    [10 14] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [33 92 150]./255 2 'mot alpha' 'noRESS';...
    [18 30] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [170 99 57]./255 2 'mot beta' 'noRESS';...
    [15 23] {'CP1';'CPz';'CP2';'P1';'Pz';'P2'} 'data_induced_bc' [23 150 118]./255 2 'parietal beta' 'noRESS';...
    [2 4] {'C3';'CP3';'C1';'CP1'} 'data_induced_bc' [0.4 0.4 0.4] 2 'mot delta' 'noRESS';...
    };


% this one for the poster figures
% pl.parameters = {...
%     [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
%     [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [255 151 0]./255 2 'vis alpha' 'noRESS';...
%     [10 14] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [33 92 150]./255 2 'mot alpha' 'noRESS';...
%     [18 30] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [170 99 57]./255 2 'mot beta' 'noRESS';...
%     [15 23] {'CP1';'CPz';'CP2';'P1';'Pz';'P2'} 'data_induced_bc' [23 150 118]./255 2 'parietal beta' 'noRESS';...
%     };



% pl.parameters = {...
%     [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
%     [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [255 151 0]./255 2 'vis alpha' 'noRESS';...
%     [10 14] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [33 92 150]./255 2 'mot alpha' 'noRESS';...
%     [15 30] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [170 99 57]./255 2 'mot beta' 'noRESS';...
%     [15 23] {'CP1';'CPz';'CP2';'P1';'Pz';'P2'} 'data_induced_bc' [23 150 118]./255 2 'parietal beta' 'noRESS';...
%     };


pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
    [8 14] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [255 151 0]./255 2 'vis alpha' 'noRESS';...
    [8 14] {'C3';'CP3'} 'data_induced_bc' [33 92 150]./255 2 'mot alpha' 'noRESS';...
    [15 20] {'C3';'CP3'} 'data_induced_bc' [170 99 57]./255 2 'mot low beta' 'noRESS';...
    [20 30] {'C3';'CP3'} 'data_induced_bc' [23 150 118]./255 2 'mot high beta' 'noRESS';...
    };

pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
    [8 14] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [255 151 0]./255 2 'vis alpha' 'noRESS';...
    [8 14] {'C3';'CP3'} 'data_induced_bc' [33 92 150]./255 2 'mot alpha' 'noRESS';...
    [15 30] {'C3';'CP3'} 'data_induced_bc' [170 99 57]./255 2 'mot beta' 'noRESS';...
    };




pl.subs2use = 1:numel(F.Subjects2Use);
% pl.subs2use = [1:14 16 :20];
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 13 14 17 18 19 20]; % subjects discarded due to low trial number version 1
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20]; % subjects discarded due to low trial number version 2
% 
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; % subjects discarded due to low trial number version 2

% pl.xlim = [-3250 3250];
pl.xlim = [-4000 4000];
pl.xlim = [-3500 3500];

pl.plotnum = [ceil(sqrt(size(pl.parameters,1))) round(sqrt(size(pl.parameters,1)))*3+5];
pl.plpos = reshape(sort(...
    [1:11:prod(pl.plotnum) 2:11:prod(pl.plotnum) 3:11:prod(pl.plotnum) 5:11:prod(pl.plotnum) 6:11:prod(pl.plotnum) 7:11:prod(pl.plotnum)]...
    ),3,[]);
pl.topopos = reshape(sort([8:11:prod(pl.plotnum) 9:11:prod(pl.plotnum) 10:11:prod(pl.plotnum) 11:11:prod(pl.plotnum)]),2,[]);
pl.permut_n = 10000;
% pl.permut_n = 100;


h.pl = [];
figs{1}=figure;
set(gcf,'Position',[100 100 800 550],'PaperPositionMode','auto')
clear pl.tdata2 h.sp h.p h.pm
for i_pl = 1:size(pl.parameters,1)
    h.sp(i_pl,1)=subplot(pl.plotnum(1),pl.plotnum(2),pl.plpos(:,i_pl));
    
    % index frequencies
    [t.t t.ind1]=min(abs(TFA.frequency-pl.parameters{i_pl,1}(1)));
    [t.t t.ind2]=min(abs(TFA.frequency-pl.parameters{i_pl,1}(2)));
    
    % index time
    [t.t t.ind3]=min(abs(TFA.time-pl.xlim(1)));
    [t.t t.ind4]=min(abs(TFA.time-pl.xlim(2)));
    
    % index electrodes
    switch  pl.parameters{i_pl,7}
        case 'noRESS'
            pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.parameters{i_pl,2}, 'UniformOutput',false)),1));
            com=sprintf('pl.data=squeeze(nanmean(nanmean(TFA.%s(t.ind1:t.ind2,t.ind3:t.ind4,pl.elec2plot_i,pl.subs2use),1),3));',pl.parameters{i_pl,3});
        case 'RESS'
            pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes_RESS.labels},x), pl.parameters{i_pl,2}, 'UniformOutput',false)),1));
            com=sprintf('pl.data=squeeze(nanmean(TFA.%s(t.ind1:t.ind2,t.ind3:t.ind4,pl.subs2use),1));',pl.parameters{i_pl,3});
        
    end
    eval(com)
    [tt.h tt.p tt.ci tt.stats]=ttest(pl.data');
    % do cluster based permutation
    [cluster_runt_diff, timecourse_runt_diff]=eeg_erpStat_clusterP(pl.data,zeros(size(pl.data)),pl.permut_n,2);
    tt.p(timecourse_runt_diff.h_corr==0)=1;
       
    % actual plot with boundedline
    pl.mdata = mean(pl.data,2)';
    pl.ddata = tt.ci - repmat(mean(pl.data,2)',2,1);
    %pl.ddata = tt.ci;
    [h.l(i_pl), h.p(i_pl)] = boundedline(TFA.time(t.ind3:t.ind4),pl.mdata,pl.ddata(1,:));
    set(h.l(i_pl),'Color',pl.parameters{i_pl,4},'LineWidth',pl.parameters{i_pl,5});
    set(h.p(i_pl),'FaceColor',pl.parameters{i_pl,4},'FaceAlpha',0.25)
    hold on
    
    % plot single subject data
%     pl.mdata = mean(pl.data,2)';
%     h.p(i_pl,:)=plot(TFA.time(t.ind3:t.ind4),pl.data,'Color',[pl.parameters{i_pl,4} 0.3],'LineWidth',0.2);
%     hold on
%     h.pm(i_pl)=plot(TFA.time(t.ind3:t.ind4),pl.mdata,'Color',[pl.parameters{i_pl,4} 1],'LineWidth',2);
    
    % significant effects?
    pl.tdata2=nan(size(pl.mdata));
    pl.tdata2(tt.p<.05)=1;
    set(gca,'ylim',[-1.1 (1+0.15+0.05*4)]...
        .*max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false))))
    pl.maxdata = max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false)));
    pl.y = 1*pl.maxdata*(1+0.05+0.05*2);
    plot(TFA.time(t.ind3:t.ind4),pl.tdata2.*pl.y,'Color',[pl.parameters{i_pl,4} 1],'LineWidth',6)
    
    
    % graphics
    xlim(pl.xlim)
    title(pl.parameters{i_pl,6},'FontSize', 8)
    grid on
    box on
    hline(0,'k')
    vline(TFA.baseline,'k')
    xlabel('time in ms')
    ylabel('amplitude modulation in %')
    
end
% %%%%%%% vers1 %start%
% set(gca,'ylim',[-1.1 1.1].*max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false))))
% %%%%%%% vers1 %end%

for i_pl = 1:size(pl.parameters,1)
    subplot(pl.plotnum(1),pl.plotnum(2),pl.topopos(:,i_pl))
    pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.parameters{i_pl,2}, 'UniformOutput',false)),1));
%     topoplot([],TFA.electrodes(1:64),'style','blank','electrodes', 'on', 'emarker2', {find(pl.elec2plot_i),'o','r',2});
    topoplot(find(pl.elec2plot_i),TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
        'emarker2',{find(pl.elec2plot_i),'o',pl.parameters{i_pl,4},3,1});
    title(sprintf('%s\n[%1.1f %1.1f]Hz',pl.parameters{i_pl,6},pl.parameters{i_pl,1}),'Color',pl.parameters{i_pl,4})
end

sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
sav.pathout = 'C:\Users\EEG\Documents\MATLAB\christopher\SSVEP_volmov\figures\';
% sav.filenames = {'Resp_AllSignals_Amp_Timecourse'};
% sav.filenames = {'Resp_AllSignals_Amp_Timecourse_v2'};
% sav.filenames = {'Resp_AllSignals_Amp_Timecourse_5'};
sav.filenames = {'Resp_AllSignals_Amp_Timecourse_sep_v2_7'};
% sav.filenames = {'Resp_AllSignals_Amp_Timecourse_v2_clustcorr_5b'};
% for i_fig = 1:1
%     print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-dpng','-r300')
%     print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-djpeg','-r300')
%     saveas(figs{i_fig},fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'fig')
%     print(figs{i_fig},fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-depsc2', '-vector','-r300')
% end

%% do complete data driven TFCE
p.e_h               = [0.66 2]; % tfce parameter
p.Samp              = 128; % data sampling rate (TFA.srate)
p.nperm             = 10000; % number of permutations, if feasible use 100000
p.plevel            = 0.05;
p.time2test         = [-3000 3500]; % time in ms

p.time2test_idx     = dsearchn(TFA.time', p.time2test');

p.tfce_olddate = '22-Mar-2024';
startdate = date;
savnam_tfce_tfa = sprintf('tfce_tfa_%1.2fe_%1.2fh_%s.mat',p.e_h(1),p.e_h(2),startdate);

% preallocate data
data4tfce2_tfa = double(TFA.data_induced_bc(:,p.time2test_idx(1):p.time2test_idx(2),:,:));
data4tfce1_tfa = zeros(size(data4tfce2_tfa));
eloc = TFA.electrodes;

if exist([savnam_tfce_tfa(1:regexp(savnam_tfce_tfa,[date '.mat'])-1) p.tfce_olddate '.mat'],'file') 
    try % try to load previous file
        load([savnam_tfce_tfa(1:regexp(savnam_tfce_tfa,[date '.mat'])-1) p.tfce_olddate '.mat'])       
    catch
        error('error when loading file')
    end
    results_tfce_tfa = Results;
else % do tfce
    % tfce on two dimensions (channels x freq x time), subjects have to be in dim 1
    % channel neighborhood calculated within ept_TFCE function

    results_tfce_tfa= ept_TFCE(permute(data4tfce1_tfa,[4 3 1 2]),permute(data4tfce2_tfa,[4 3 1 2]),...
        eloc,'rsample',p.Samp,...
        'type','d','plots',0,...
        'e_h',p.e_h,'nperm',p.nperm,...
        'flag_tfce',true,'flag_ft',false,... %channel testing
        'savename',savnam_tfce_tfa, ...
        'flag_save',true);
end

ChN = ept_ChN2(eloc);
res.cluster_tfa = ept_calculateClusters(results_tfce_tfa,ChN,0.05);


% what do I want to show
% - frequency spectrum
% - time frequency plot (for peak electrodes)?
% - electrode by time plot showing significant differences
% - topography: plotting significant timepoints by electrode

t.data = permute(data4tfce2_tfa,[4 3 1 2]);
t.data_time = TFA.time(p.time2test_idx(1):p.time2test_idx(2));
clear h
for i_cluster = 1:numel(res.cluster_tfa)
    
    h.fig(i_cluster)=figure;
    set(gcf,'Position',[100 100 1300 300],'PaperPositionMode','auto')
    tiledlayout(1,7,'TileSpacing','compact')

    t.cluster_loc = repmat(res.cluster_tfa(i_cluster).cluster_locations,[size(data4tfce2_tfa,4),1,1,1]);
    t.cluster_chan = squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,[2,3]));
    t.cluster_time = squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,[1,2]));
    t.cluster_freq = squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,[1,3]));
    
    % frequency spectrum for peak electrode
    pl.elec2plot_now_i = res.cluster_tfa(i_cluster).channel_peak;
    t.data2 = t.data; t.data2(~t.cluster_loc) = nan;
    pl.cluster_data = squeeze(mean(t.data2,[1,2,4],"omitnan"));
    pl.clusterElec_data = squeeze(mean(t.data(:,t.cluster_chan,:,t.cluster_time),[1 2 4]));
%     pl.clusterElec_data = squeeze(mean(t.data(:,pl.elec2plot_now_i,:,t.cluster_time),[1 2 4]));

    nexttile([1 2])
    % plot spectrum at frequency
    plot(TFA.frequency,pl.clusterElec_data,'Color',[0.3 0.3 0.3])
    hold on
    plot(TFA.frequency,pl.cluster_data,'Color',[0.9 0.3 0.3])
    xlabel('frequency in Hz')
    ylabel('modulation in %')
    xlim([TFA.frequency(1) TFA.frequency(end)])
    title(sprintf('spectra | comp %1.0f | [%1.0f %1.0f]ms', ...
        i_cluster, t.data_time(find(diff(t.cluster_time)~=0))))
    legend({sprintf('peak electrode %s',TFA.electrodes(pl.elec2plot_now_i).labels);'only sign. data'},'Location','southoutside','Orientation','horizontal')
    grid on
    

    % time-frequency-plot for peak electrode
    pl.elec2plot_now_i = res.cluster_tfa(i_cluster).channel_peak;
    
    pl.plotdata = squeeze(mean(t.data(:,pl.elec2plot_now_i,:,:),[1,2]));
    pl.clim = [-1 1]*max(abs(pl.plotdata),[],'all');

    pl.plotdata_sign = squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,1));
    
    nexttile([1 2])
    imagesc(t.data_time,TFA.frequency,pl.plotdata,pl.clim)
    colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    hold on
    contour(t.data_time,TFA.frequency,pl.plotdata_sign,'EdgeColor','g')
    title(sprintf('TFA at peak elec. %s | cluster in green (for all elecs)', ...
        TFA.electrodes(pl.elec2plot_now_i).labels))
    xlabel('time in ms')
    ylabel('frequency in Hz')
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();

    % time by electrode plot
    pl.plotdata = squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,2));

    nexttile([1 2])
    imagesc(t.data_time,1:64,pl.plotdata,[-1 1])
    colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    set(gca,'YTick',0:5:64,'YTickLabel',{TFA.electrodes(5:5:64).labels})
    title('electrodes in cluster across time')

    % electrode position as summed values

    nexttile([1 1])
    pl.plotdata = sum(squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,2)),2);
    topoplot(pl.plotdata, eloc, ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',[0 max(pl.plotdata)], ...
                'colormap',cbrewer2('PuBuGn'),'whitebk','on');
        title(sprintf('significant electrodes in cluster'))

    colorbar




    
end

%% do complete data driven TFCE | for frequency band alpha
p.e_h               = [0.66 2]; % tfce parameter
p.Samp              = 128; % data sampling rate (TFA.srate)
p.nperm             = 10000; % number of permutations, if feasible use 100000
p.plevel            = 0.05;
p.time2test         = [-3000 3500]; % time in ms
p.freq2test         = [8 14]; % frequency in Hz

p.time2test_idx     = dsearchn(TFA.time', p.time2test');
p.freq2test_idx     = dsearchn(TFA.frequency', p.freq2test');


p.tfce_olddate = '26-Mar-2024';
startdate = date;
savnam_tfce_tfa = sprintf('tfce_tfa_freqrange_alpha%1.2fe_%1.2fh_%s.mat',p.e_h(1),p.e_h(2),startdate);

% preallocate data
data4tfce2_tfa = double(TFA.data_induced_bc(p.freq2test_idx(1):p.freq2test_idx(2),p.time2test_idx(1):p.time2test_idx(2),:,:));
data4tfce1_tfa = zeros(size(data4tfce2_tfa));
eloc = TFA.electrodes;

if exist([savnam_tfce_tfa(1:regexp(savnam_tfce_tfa,[date '.mat'])-1) p.tfce_olddate '.mat'],'file') 
    try % try to load previous file
        load([savnam_tfce_tfa(1:regexp(savnam_tfce_tfa,[date '.mat'])-1) p.tfce_olddate '.mat'])       
    catch
        error('error when loading file')
    end
    results_tfce_tfa = Results;
else % do tfce
    % tfce on two dimensions (channels x freq x time), subjects have to be in dim 1
    % channel neighborhood calculated within ept_TFCE function

    results_tfce_tfa= ept_TFCE(permute(data4tfce1_tfa,[4 3 1 2]),permute(data4tfce2_tfa,[4 3 1 2]),...
        eloc,'rsample',p.Samp,...
        'type','d','plots',0,...
        'e_h',p.e_h,'nperm',p.nperm,...
        'flag_tfce',true,'flag_ft',false,... %channel testing
        'savename',savnam_tfce_tfa, ...
        'flag_save',true);
end

ChN = ept_ChN2(eloc);
res.cluster_tfa = ept_calculateClusters(results_tfce_tfa,ChN,0.05);


% what do I want to show
% - frequency spectrum
% - time frequency plot (for peak electrodes)?
% - electrode by time plot showing significant differences
% - topography: plotting significant timepoints by electrode

t.data = permute(data4tfce2_tfa,[4 3 1 2]);
t.data_time = TFA.time(p.time2test_idx(1):p.time2test_idx(2));
clear h
for i_cluster = 1:numel(res.cluster_tfa)
    
    h.fig(i_cluster)=figure;
    set(gcf,'Position',[100 100 1300 300],'PaperPositionMode','auto')
    tiledlayout(1,7,'TileSpacing','compact')

    t.cluster_loc = repmat(res.cluster_tfa(i_cluster).cluster_locations,[size(data4tfce2_tfa,4),1,1,1]);
    t.cluster_chan = squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,[2,3]));
    t.cluster_time = squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,[1,2]));
    t.cluster_freq = squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,[1,3]));
    
    % frequency spectrum for peak electrode
    pl.elec2plot_now_i = res.cluster_tfa(i_cluster).channel_peak;
    t.data2 = t.data; t.data2(~t.cluster_loc) = nan;
    pl.cluster_data = squeeze(mean(t.data2,[1,2,4],"omitnan"));
    pl.clusterElec_data = squeeze(mean(t.data(:,t.cluster_chan,:,t.cluster_time),[1 2 4]));
%     pl.clusterElec_data = squeeze(mean(t.data(:,pl.elec2plot_now_i,:,t.cluster_time),[1 2 4]));

    nexttile([1 2])
    % plot spectrum at frequency
    plot(TFA.frequency(p.freq2test_idx(1):p.freq2test_idx(2)),pl.clusterElec_data,'Color',[0.3 0.3 0.3])
    hold on
    plot(TFA.frequency(p.freq2test_idx(1):p.freq2test_idx(2)),pl.cluster_data,'Color',[0.9 0.3 0.3])
    xlabel('frequency in Hz')
    ylabel('modulation in %')
    xlim([TFA.frequency(1) TFA.frequency(end)])
    title(sprintf('spectra | freq tested [%1.0f %1.0f]Hz | comp %1.0f | [%1.0f %1.0f]ms', ...
        p.freq2test, i_cluster, t.data_time(find(diff(t.cluster_time)~=0))))
%     legend({sprintf('peak electrode %s',TFA.electrodes(pl.elec2plot_now_i).labels);'only sign. data'},'Location','southoutside','Orientation','horizontal')
    legend({sprintf('all cluster elecs');'only sign. data'},'Location','southoutside','Orientation','horizontal')
    grid on
    

    % time-frequency-plot for peak electrode
%     pl.elec2plot_now_i = res.cluster_tfa(i_cluster).channel_peak;
    pl.elec2plot_now_i = t.cluster_chan;
    
    pl.plotdata = squeeze(mean(t.data(:,pl.elec2plot_now_i,:,:),[1,2]));
    pl.clim = [-1 1]*max(abs(pl.plotdata),[],'all');

    pl.plotdata_sign = squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,1));
    % add zeros around
    pl.plotdata_sign = [zeros(1,size(pl.plotdata_sign,2)); pl.plotdata_sign; zeros(1,size(pl.plotdata_sign,2))];
    pl.freq2plot = TFA.frequency(p.freq2test_idx(1)-1:p.freq2test_idx(2)+1);
    
    nexttile([1 2])
    imagesc(t.data_time,TFA.frequency,pl.plotdata,pl.clim)
    colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    hold on
    contour(t.data_time,pl.freq2plot,pl.plotdata_sign,'EdgeColor','g')
%     title(sprintf('TFA at peak elec. %s | cluster in green (for all elecs)', ...
%         TFA.electrodes(pl.elec2plot_now_i).labels))
    title(sprintf('all cluster channels | cluster in green (for all elecs)'))
    xlabel('time in ms')
    ylabel('frequency in Hz')
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();

    % time by electrode plot
    pl.plotdata = squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,2));

    nexttile([1 2])
    imagesc(t.data_time,1:64,pl.plotdata,[-1 1])
    colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    set(gca,'YTick',0:5:64,'YTickLabel',{TFA.electrodes(5:5:64).labels})
    title('electrodes in cluster across time')

    % electrode position as summed values

    nexttile([1 1])
    pl.plotdata = sum(squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,2)),2);
    pl.colormap = [0.5 0.5 0.5; cbrewer2('PuBuGn')]; % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    
    topoplot(pl.plotdata, eloc, ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',[0 max(pl.plotdata)], ...
                'colormap',pl.colormap,'whitebk','on');
        title(sprintf('significant electrodes in cluster'))

    colorbar

end


%% do complete data driven TFCE | for frequency band beta
p.e_h               = [0.66 2]; % tfce parameter
p.Samp              = 128; % data sampling rate (TFA.srate)
p.nperm             = 10000; % number of permutations, if feasible use 100000
p.plevel            = 0.05;
p.time2test         = [-3000 3500]; % time in ms
p.freq2test         = [15 30]; % frequency in Hz

p.time2test_idx     = dsearchn(TFA.time', p.time2test');
p.freq2test_idx     = dsearchn(TFA.frequency', p.freq2test');


p.tfce_olddate = '26-Mar-2024';
startdate = date;
savnam_tfce_tfa = sprintf('tfce_tfa_freqrange_beta%1.2fe_%1.2fh_%s.mat',p.e_h(1),p.e_h(2),startdate);

% preallocate data
data4tfce2_tfa = double(TFA.data_induced_bc(p.freq2test_idx(1):p.freq2test_idx(2),p.time2test_idx(1):p.time2test_idx(2),:,:));
data4tfce1_tfa = zeros(size(data4tfce2_tfa));
eloc = TFA.electrodes;

if exist([savnam_tfce_tfa(1:regexp(savnam_tfce_tfa,[date '.mat'])-1) p.tfce_olddate '.mat'],'file') 
    try % try to load previous file
        load([savnam_tfce_tfa(1:regexp(savnam_tfce_tfa,[date '.mat'])-1) p.tfce_olddate '.mat'])       
    catch
        error('error when loading file')
    end
    results_tfce_tfa = Results;
else % do tfce
    % tfce on two dimensions (channels x freq x time), subjects have to be in dim 1
    % channel neighborhood calculated within ept_TFCE function

    results_tfce_tfa= ept_TFCE(permute(data4tfce1_tfa,[4 3 1 2]),permute(data4tfce2_tfa,[4 3 1 2]),...
        eloc,'rsample',p.Samp,...
        'type','d','plots',0,...
        'e_h',p.e_h,'nperm',p.nperm,...
        'flag_tfce',true,'flag_ft',false,... %channel testing
        'savename',savnam_tfce_tfa, ...
        'flag_save',true);
end

ChN = ept_ChN2(eloc);
res.cluster_tfa = ept_calculateClusters(results_tfce_tfa,ChN,0.05);


% what do I want to show
% - frequency spectrum
% - time frequency plot (for peak electrodes)?
% - electrode by time plot showing significant differences
% - topography: plotting significant timepoints by electrode

t.data = permute(data4tfce2_tfa,[4 3 1 2]);
t.data_time = TFA.time(p.time2test_idx(1):p.time2test_idx(2));
clear h
for i_cluster = 1:numel(res.cluster_tfa)
    
    h.fig(i_cluster)=figure;
    set(gcf,'Position',[100 100 1300 300],'PaperPositionMode','auto')
    tiledlayout(1,7,'TileSpacing','compact')

    t.cluster_loc = repmat(res.cluster_tfa(i_cluster).cluster_locations,[size(data4tfce2_tfa,4),1,1,1]);
    t.cluster_chan = squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,[2,3]));
    t.cluster_time = squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,[1,2]));
    t.cluster_freq = squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,[1,3]));
    
    % frequency spectrum for peak electrode
    pl.elec2plot_now_i = res.cluster_tfa(i_cluster).channel_peak;
    t.data2 = t.data; t.data2(~t.cluster_loc) = nan;
    pl.cluster_data = squeeze(mean(t.data2,[1,2,4],"omitnan"));
    pl.clusterElec_data = squeeze(mean(t.data(:,t.cluster_chan,:,t.cluster_time),[1 2 4]));
%     pl.clusterElec_data = squeeze(mean(t.data(:,pl.elec2plot_now_i,:,t.cluster_time),[1 2 4]));

    nexttile([1 2])
    % plot spectrum at frequency
    plot(TFA.frequency(p.freq2test_idx(1):p.freq2test_idx(2)),pl.clusterElec_data,'Color',[0.3 0.3 0.3])
    hold on
    plot(TFA.frequency(p.freq2test_idx(1):p.freq2test_idx(2)),pl.cluster_data,'Color',[0.9 0.3 0.3])
    xlabel('frequency in Hz')
    ylabel('modulation in %')
    xlim([TFA.frequency(1) TFA.frequency(end)])
    title(sprintf('spectra | freq tested [%1.0f %1.0f]Hz | comp %1.0f | [%1.0f %1.0f]ms', ...
        p.freq2test, i_cluster, t.data_time(find(diff(t.cluster_time)~=0))))
%     legend({sprintf('peak electrode %s',TFA.electrodes(pl.elec2plot_now_i).labels);'only sign. data'},'Location','southoutside','Orientation','horizontal')
    legend({sprintf('all cluster elecs');'only sign. data'},'Location','southoutside','Orientation','horizontal')
    grid on
    

    % time-frequency-plot for peak electrode
%     pl.elec2plot_now_i = res.cluster_tfa(i_cluster).channel_peak;
    pl.elec2plot_now_i = t.cluster_chan;
    
    pl.plotdata = squeeze(mean(t.data(:,pl.elec2plot_now_i,:,:),[1,2]));
    pl.clim = [-1 1]*max(abs(pl.plotdata),[],'all');

    pl.plotdata_sign = squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,1));
    % add zeros around
    pl.plotdata_sign = [zeros(1,size(pl.plotdata_sign,2)); pl.plotdata_sign; zeros(1,size(pl.plotdata_sign,2))];
    pl.freq2plot = TFA.frequency(p.freq2test_idx(1)-1:p.freq2test_idx(2)+1);
    
    nexttile([1 2])
    imagesc(t.data_time,TFA.frequency,pl.plotdata,pl.clim)
    colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    hold on
    contour(t.data_time,pl.freq2plot,pl.plotdata_sign,'EdgeColor','g')
%     title(sprintf('TFA at peak elec. %s | cluster in green (for all elecs)', ...
%         TFA.electrodes(pl.elec2plot_now_i).labels))
    title(sprintf('all cluster channels | cluster in green (for all elecs)'))
    xlabel('time in ms')
    ylabel('frequency in Hz')
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();

    % time by electrode plot
    pl.plotdata = squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,2));

    nexttile([1 2])
    imagesc(t.data_time,1:64,pl.plotdata,[-1 1])
    colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    set(gca,'YTick',0:5:64,'YTickLabel',{TFA.electrodes(5:5:64).labels})
    title('electrodes in cluster across time')

    % electrode position as summed values

    nexttile([1 1])
    pl.plotdata = sum(squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,2)),2);
    pl.colormap = [0.5 0.5 0.5; cbrewer2('PuBuGn')]; % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    
    topoplot(pl.plotdata, eloc, ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',[0 max(pl.plotdata)], ...
                'colormap',pl.colormap,'whitebk','on');
        title(sprintf('significant electrodes in cluster'))

    colorbar

end


%% do frequency band collapsed data driven TFCE
p.e_h               = [0.66 2]; % tfce parameter
p.Samp              = 128; % data sampling rate (TFA.srate)
p.nperm             = 10000; % number of permutations, if feasible use 100000
p.plevel            = 0.05;
p.time2test         = [-3000 3250]; % time in ms
% p.freq2test         = [8 12]; % frequency in Hz alpha
p.freq2test         = [8 14]; % frequency in Hz alpha
% p.freq2test         = [10 14]; % frequency in Hz alpha
% p.freq2test         = [15 20]; % frequency in Hz low beta
% p.freq2test         = [20 30]; % frequency in Hz high beta


p.time2test_idx     = dsearchn(TFA.time', p.time2test');
p.freq2test_idx     = dsearchn(TFA.frequency', p.freq2test');

p.overwrite         = false; % redo calculation


% p.tfce_olddate = '26-Mar-2024';
startdate = date;
savnam_tfce_tfa = sprintf('tfce_tfa_freqcollapsed_%1.0f_to_%1.0f_Hz_%1.2fe_%1.2fh_%s.mat',p.freq2test,p.e_h(1),p.e_h(2),startdate);

% preallocate data
data4tfce2_tfa= double(squeeze(mean(TFA.data_induced_bc(p.freq2test_idx(1):p.freq2test_idx(2),p.time2test_idx(1):p.time2test_idx(2),:,:),1)));
data4tfce1_tfa = zeros(size(data4tfce2_tfa));
eloc = TFA.electrodes;

if ~isempty(dir([savnam_tfce_tfa(1:regexp(savnam_tfce_tfa,[date '.mat'])-1) '*.mat'])) && ~p.overwrite 
    try % try to load previous file
        % load most recent file
        t.dir = dir([savnam_tfce_tfa(1:regexp(savnam_tfce_tfa,[date '.mat'])-1) '*.mat']);
        [~, t.fileidx] = max([t.dir.datenum]);
        load(t.dir(t.fileidx).name)       
    catch
        error('error when loading file')
    end
    results_tfce_tfa = Results;
else % do tfce
    % tfce on two dimensions (channels x time), subjects have to be in dim 1
    % channel neighborhood calculated within ept_TFCE function

    results_tfce_tfa= ept_TFCE(permute(data4tfce1_tfa,[3 2 1 ]),permute(data4tfce2_tfa,[3 2 1]),...
        eloc,'rsample',p.Samp,...
        'type','d','plots',0,...
        'e_h',p.e_h,'nperm',p.nperm,...
        'flag_tfce',true,'flag_ft',false,... %channel testing
        'savename',savnam_tfce_tfa, ...
        'flag_save',true);
end

ChN = ept_ChN2(eloc);
res.cluster_tfa = ept_calculateClusters(results_tfce_tfa,ChN,0.05);


% what do I want to show
% - plot (for peak electrodes)?
% - electrode by time plot showing significant differences
% - topography: plotting significant timepoints by electrode

t.data = permute(data4tfce2_tfa,[3 2 1 ]);
t.data_tfa = squeeze(TFA.data_induced_bc(:,p.time2test_idx(1):p.time2test_idx(2),:,:));
t.data_time = TFA.time(p.time2test_idx(1):p.time2test_idx(2));
clear h
for i_cluster = 1:numel(res.cluster_tfa)
    
    h.fig(i_cluster)=figure;
    set(gcf,'Position',[100 100 1300 300],'PaperPositionMode','auto')
    tiledlayout(1,7,'TileSpacing','compact')

    t.cluster_loc = repmat(res.cluster_tfa(i_cluster).cluster_locations,[size(data4tfce2_tfa,3),1,1]);
    t.cluster_chan = squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,[2]));
    t.cluster_time = squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,[1]));
    

    % show TFA for cluster in
    %     pl.elec2plot_now_i = res.cluster_tfa(i_cluster).channel_peak;
    pl.elec2plot_now_i = t.cluster_chan;
    
    pl.plotdata = squeeze(mean(t.data_tfa(:,:,pl.elec2plot_now_i,:),[3,4]));
    pl.clim = [-1 1]*max(abs(pl.plotdata),[],'all');

    pl.plotdata_sign = squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,1));
    % add zeros around
    pl.plotdata_sign = [zeros(1,size(pl.plotdata_sign,2)); ...
        repmat(pl.plotdata_sign,diff(p.freq2test_idx)+1,1); zeros(1,size(pl.plotdata_sign,2))];
    pl.freq2plot = TFA.frequency(p.freq2test_idx(1)-1:p.freq2test_idx(2)+1);
    
    nexttile([1 2])
    imagesc(t.data_time,TFA.frequency,pl.plotdata,pl.clim)
    colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    hold on
    contour(t.data_time,pl.freq2plot,pl.plotdata_sign,'EdgeColor','g')
%     title(sprintf('TFA at peak elec. %s | cluster in green (for all elecs)', ...
%         TFA.electrodes(pl.elec2plot_now_i).labels))
    title(sprintf('all cluster channels | cluster in green (for all elecs)'))
    xlabel('time in ms')
    ylabel('frequency in Hz')
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();
    
        

    % time-frequency-plot for peak electrode
%     pl.elec2plot_now_i = res.cluster_tfa(i_cluster).channel_peak;
    pl.elec2plot_now_i = t.cluster_chan;
    
    pl.plotdata = squeeze(mean(t.data(:,pl.elec2plot_now_i,:,:),[1,2]));
    pl.clim = [-1 1]*max(abs(pl.plotdata),[],'all');

    pl.plotdata_sign = squeeze(any(res.cluster_tfa(i_cluster).cluster_locations,1));
    pl.plotdata_sign_data = pl.plotdata;  pl.plotdata_sign_data(~pl.plotdata_sign)=nan;
    
    nexttile([1 2])
    plot(t.data_time,pl.plotdata,'Color',[0.3 0.3 0.3])
    hold on
    plot(t.data_time,pl.plotdata_sign_data,'Color',[0.9 0.3 0.3],'LineWidth',2)
    xlabel('time in ms')
    ylabel('amplitude in \muV/m²')
%     title(sprintf('time course at peak elec. %s | freq [%1.0f %1.0f]Hz', ...
%         TFA.electrodes(pl.elec2plot_now_i).labels, p.freq2test ))
    title(sprintf('time course at cluster electrodes | freq [%1.0f %1.0f]Hz collapsed', ...
        p.freq2test ))
    grid on

    % time by electrode plot
    pl.plotdata = res.cluster_tfa(i_cluster).cluster_locations;

    nexttile([1 2])
    imagesc(t.data_time,1:64,pl.plotdata,[-1 1])
    colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YTick',0:5:64,'YTickLabel',{TFA.electrodes(5:5:64).labels})
    set(gca,'YDir','normal')
    title('electrodes in cluster across time')

    % electrode position as summed values

    nexttile([1 1])
    pl.plotdata = sum(res.cluster_tfa(i_cluster).cluster_locations,2);
    pl.colormap = [0.5 0.5 0.5; cbrewer2('PuBuGn')]; % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    topoplot(pl.plotdata, eloc, ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',[0 max(pl.plotdata)], ...
                'colormap',pl.colormap,'whitebk','on');
        title(sprintf('significant electrodes in cluster'))

    colorbar




    
end




%% do electrode collapsed data driven cluser correction
p.e_h               = [0.66 2]; % tfce parameter
p.Samp              = 128; % data sampling rate (TFA.srate)
p.nperm             = 10000; % number of permutations, if feasible use 100000
p.plevel            = 0.05;
p.time2test         = [-3000 3250]; % time in ms
p.time2plot         = [-3250 3250]; % time in ms
% p.elec2test         = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; pl.elecname = 'visuooccipital';%visual 
% p.elec2test         = {'P9';'P10';'PO7';'P7';'P8';'P5';'P6';'PO3';'PO4';'O1';'O2';'I1';'I2'}; pl.elecname = 'visuooccipittemporal';%visual 
% p.elec2test         = {'C3';'CP3';'C5';'CP5'}; pl.elecname = 'leftmotor'; % motor
p.elec2test         = {'C3';'CP3'}; pl.elecname = 'smallleftmotor'; % motor

pl.freq2test        = [3 35];
pl.freq2test_idx    = dsearchn(TFA.frequency',pl.freq2test');

p.time2test_idx     = dsearchn(TFA.time', p.time2test');
pl.elec2test_i      = logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), p.elec2test, 'UniformOutput',false)),1));

p.overwrite         = false; % redo calculation
% p.overwrite         = true; % redo calculation


% p.tfce_olddate = '26-Mar-2024';
startdate = datetime("today");
savnam_tfce_tfa = sprintf('clustercorrection_tfa_eleccollapsed_%s_%s.mat',pl.elecname,startdate);

% preallocate data
data4tfce2_tfa= double(squeeze(mean( ...
    TFA.data_induced_bc(pl.freq2test_idx(1):pl.freq2test_idx(2),p.time2test_idx(1):p.time2test_idx(2),pl.elec2test_i,:),3)));
data4tfce1_tfa = zeros(size(data4tfce2_tfa));
eloc = TFA.electrodes(1);

clear res

if ~isempty(dir([savnam_tfce_tfa(1:regexp(savnam_tfce_tfa,[date '.mat'])-1) '*.mat'])) && ~p.overwrite 
    try % try to load previous file
        % load most recent file
        t.dir = dir([savnam_tfce_tfa(1:regexp(savnam_tfce_tfa,[date '.mat'])-1) '*.mat']);
        [~, t.fileidx] = max([t.dir.datenum]);
        load(t.dir(t.fileidx).name)       
    catch
        error('error when loading file')
    end
else % do cluster correction
    % first comparison of real data
    res.data4cc2_tfa = data4tfce2_tfa;
    res.data4cc1_tfa = data4tfce1_tfa;
    
    % t.test
    [tt.h,tt.p,tt.ci,tt.stats] = ttest(res.data4cc2_tfa,res.data4cc1_tfa,'Dim',3);

    % extract data
    res.data.tt.h = tt.h; res.data.tt.p = tt.p; res.data.tt.tstat = tt.stats.tstat;

    % do cluster calculaten
    res.data.tt.clusterlabel = bwlabel(res.data.tt.h);
    t.clusternum = 1:max(res.data.tt.clusterlabel,[],'all');

    % extract cluster information
    for i_cluster = 1:max(t.clusternum)
        % index cluster
        res.data.clust(i_cluster).clusternum = i_cluster;
        res.data.clust(i_cluster).clustersize = sum(res.data.tt.clusterlabel==i_cluster,"all");
        [res.data.clust(i_cluster).freqrange(1), res.data.clust(i_cluster).freqrange(2)] = ...
            bounds(TFA.frequency(any(res.data.tt.clusterlabel==i_cluster,2)));
        [res.data.clust(i_cluster).timerange(1), res.data.clust(i_cluster).timerange(2)] = ...
            bounds(TFA.time(find(any(res.data.tt.clusterlabel==i_cluster,1))+p.time2test_idx(1)-1));
        res.data.clust(i_cluster).tsum = sum(res.data.tt.tstat(res.data.tt.clusterlabel==i_cluster),"all");
        if res.data.clust(i_cluster).tsum<0
            res.data.clust(i_cluster).tpeak = min(res.data.tt.tstat(res.data.tt.clusterlabel==i_cluster),[],"all");
        else
            res.data.clust(i_cluster).tpeak = max(res.data.tt.tstat(res.data.tt.clusterlabel==i_cluster),[],"all");
        end
        res.data.clust(i_cluster).ppeak = min(res.data.tt.p(res.data.tt.clusterlabel==i_cluster),[],"all");
    end

    % now do the same thing for the permutation
    % first append data (to sample from null distribution)
    res.data4cc_permute = cat(3,res.data4cc2_tfa, res.data4cc1_tfa);
    % all combinations of reaarrangement
    t.permidx = [false(size(res.data4cc_permute,3)/2,1); true(size(res.data4cc_permute,3)/2,1)];
    t.permidxs = repmat(false(size(res.data4cc_permute,3),1),1,p.nperm);
    for i_rand = 1:p.nperm
        t.permidxs(:,i_rand) = t.permidx(randperm(numel(t.permidx)));
        % create unique solutions
        while i_rand > 1 & ismember(t.permidxs(:,i_rand)',t.permidxs(:,1:i_rand-1)','rows')
            t.permidxs(:,i_rand) = t.permidx(randperm(numel(t.permidx)))
        end
    end
    % now do the loop across these permutations
    res.perm.clust = [];
    fprintf('running %1.0f permutations | progress in %%: 0',p.nperm)
    t.update_num = linspace(0,p.nperm,11);
    t.update_perc = linspace(0,100,11);
    for i_rand = 1:p.nperm
        % display progress
        if any(i_rand == t.update_num)
            fprintf('%3.0f ',t.update_perc(i_rand == t.update_num))
        end
        % do ttest of data
        [tt.h,tt.p,tt.ci,tt.stats] = ttest( ...
            res.data4cc_permute(:,:,t.permidxs(:,i_rand)), ...
            res.data4cc_permute(:,:,~t.permidxs(:,i_rand)), ...
            'Dim',3);
        
        % do cluster calculaten
        t.tt.clusterlabel = bwlabel(tt.h);
        t.clusternum = 1:max(t.tt.clusterlabel,[],'all');

        % extract cluster information
        t.clust = [];
        for i_cluster = 1:max(t.clusternum)
            % index cluster
            t.clust(i_cluster).clusternum = i_cluster;
            t.clust(i_cluster).clustersize = sum(t.tt.clusterlabel==i_cluster,"all");
            t.clust(i_cluster).tsum = sum(tt.stats.tstat(t.tt.clusterlabel==i_cluster),"all");
            if t.clust(i_cluster).tsum<0
                t.clust(i_cluster).tpeak = min(tt.stats.tstat(t.tt.clusterlabel==i_cluster),[],"all");
            else
                t.clust(i_cluster).tpeak = max(tt.stats.tstat(t.tt.clusterlabel==i_cluster),[],"all");
            end
            t.clust(i_cluster).ppeak = min(tt.p(t.tt.clusterlabel==i_cluster),[],"all");
        end

        % extract important info
        res.perm.clust(i_rand).maxclustersize = max([t.clust.clustersize]);
        res.perm.clust(i_rand).maxtsum = max([t.clust.tsum]);
        res.perm.clust(i_rand).maxtpeak = max([t.clust.tpeak]);
        res.perm.clust(i_rand).mintsum = min([t.clust.tsum]);
        res.perm.clust(i_rand).mintpeak = min([t.clust.tpeak]);
        res.perm.clust(i_rand).minppeak = min([t.clust.ppeak]);

        % graphically check values
        %figure; histogram([res.perm.clust.maxtsum],50)
    end
    fprintf('...done!\n')

    % now compare each of the clusters and create empirical p-values
    for i_cluster = 1:size(res.data.clust,2)
        % empirical p-value for clustersize
        res.data.clust(i_cluster).pMC_clustersize = ...
            sum([res.perm.clust.maxclustersize]>res.data.clust(i_cluster).clustersize)/ ...
            (p.nperm + 1);
        res.data.clust(i_cluster).hMC_clustersize = ...
             res.data.clust(i_cluster).pMC_clustersize < p.plevel;
        % empirical p-value for tsum
        if res.data.clust(i_cluster).tsum > 0
            res.data.clust(i_cluster).pMC_tsum = ...
                sum([res.perm.clust.maxtsum]>res.data.clust(i_cluster).tsum)/ ...
                (p.nperm + 1);
        else
            res.data.clust(i_cluster).pMC_tsum = ...
                sum([res.perm.clust.mintsum]<res.data.clust(i_cluster).tsum)/ ...
                (p.nperm + 1);
        end
        res.data.clust(i_cluster).hMC_tsum = ...
             res.data.clust(i_cluster).pMC_tsum < p.plevel;
        % empirical p-value for tpeak
        if res.data.clust(i_cluster).tsum > 0
            res.data.clust(i_cluster).pMC_tpeak = ...
                sum([res.perm.clust.maxtpeak]>res.data.clust(i_cluster).tpeak)/ ...
                (p.nperm + 1);
        else
            res.data.clust(i_cluster).pMC_tpeak = ...
                sum([res.perm.clust.mintpeak]<res.data.clust(i_cluster).tpeak)/ ...
                (p.nperm + 1);
        end
        res.data.clust(i_cluster).hMC_tpeak = ...
             res.data.clust(i_cluster).pMC_tpeak < p.plevel;
    end

    % graphical check
    % figure; histogram([res.perm.clust.mintpeak],100); hold on; histogram([res.perm.clust.maxtpeak],100)
    % figure; histogram([res.perm.clust.mintsum],100); hold on; histogram([res.perm.clust.maxtsum],100)
    
    % save current file
    res.p = p;
    res.date = startdate;
    res = rmfield(res,'data4cc_permute');
    save(savnam_tfce_tfa,'res')
    
end

% plot TFA with clursters
figure;
pl.data_ind = mean(res.data4cc2_tfa,3);
pl.clim = [-1 1]*max(abs(pl.data_ind),[],'all');
imagesc(TFA.time(p.time2test_idx(1):p.time2test_idx(2)),TFA.frequency(pl.freq2test_idx(1):pl.freq2test_idx(2)),pl.data_ind,pl.clim)
colormap(gca, flipud(cbrewer2('RdBu')))
set(gca,'YDir','normal')
hold on
% index significant cluster 
t.signum = [res.data.clust([res.data.clust.hMC_tsum]).clusternum]; % based on t sum
% t.signum = [res.data.clust([res.data.clust.hMC_tpeak]).clusternum]; % based on t peak
% t.signum = [res.data.clust.clusternum]; % uncorrected
pl.plotdata_cluster = ismember(res.data.tt.clusterlabel,t.signum);
contour(TFA.time(p.time2test_idx(1):p.time2test_idx(2)),TFA.frequency(pl.freq2test_idx(1):pl.freq2test_idx(2)),pl.plotdata_cluster,'EdgeColor','g')

title(sprintf('bc tfa noRESS, for ROI %s | cluster corrected;\n %s', ...
    pl.elecname, vararg2str(p.elec2test)), 'FontSize',8)
colorbar
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(p.time2test)
% ylim(pl.flims)
hline(14.16667,'m')
set(gca,'FontSize',8)










    


%% topoplot
% pl.time2plot=[500 2000];
% pl.time2plot=[-1000 0];
pl.time2plot=[-500 300];
% pl.time2plot=[-2000 500];
% pl.time2plot=[-2000 750];
[t.t t.ind1]=min(abs(TFA.time-pl.time2plot(1)));
[t.t t.ind2]=min(abs(TFA.time-pl.time2plot(2)));

% pl.freq2plot=[14.16667 14.16667];
% pl.freq2plot=[10 14];
% pl.freq2plot=[15 25];
% pl.freq2plot=[20 30];
pl.freq2plot=[8 12];
[t.t t.ind3]=min(abs(TFA.frequency-pl.freq2plot(1)));
[t.t t.ind4]=min(abs(TFA.frequency-pl.freq2plot(2)));

pl.data_raw=squeeze(nanmean(nanmean(TFA.data_induced(t.ind3:t.ind4,t.ind1:t.ind2,:,pl.subs2use),1),2));
pl.data_raw(:,:,2)=squeeze(nanmean(nanmean(TFA.data_evoked(t.ind3:t.ind4,t.ind1:t.ind2,:,pl.subs2use),1),2));
pl.data_bc=squeeze(nanmean(nanmean(TFA.data_induced_bc(t.ind3:t.ind4,t.ind1:t.ind2,:,pl.subs2use),1),2));
pl.data_bc(:,:,2)=squeeze(nanmean(nanmean(TFA.data_evoked_bc(t.ind3:t.ind4,t.ind1:t.ind2,:,pl.subs2use),1),2));

pl.subs2use = 1:numel(F.Subjects2Use);

pl.titleadd = {'induced';'evoked'};

figure;
% loop for evoked/induced
for i_dat = 1:2
    h.s1=subplot(2,3,1+3*(i_dat-1));
    pl.data=squeeze(nanmean(pl.data_raw(:,:,i_dat),2));
    topoplot( pl.data, TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',[0 max(pl.data)],'conv','on','colormap',fake_parula);
    title(sprintf('raw %s amp %1.2f to %1.2f Hz;\n%1.0f to %1.0f ms',pl.titleadd{i_dat},pl.freq2plot, pl.time2plot))
%     freezeColors
    h.c1 = colorbar;
%     h.c1=cbfreeze(h.c1);
    
    h.s2=subplot(2,3,2+3*(i_dat-1));
    pl.data=squeeze(nanmean(pl.data_bc(:,:,i_dat),2));
    topoplot( pl.data, TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits','absmax','conv','on','colormap',colormap(gca,flipud(cbrewer2('RdBu'))));
    title(sprintf('normalized %s amp %1.2f to %1.2f Hz;\n%1.0f to %1.0f ms',pl.titleadd{i_dat},pl.freq2plot, pl.time2plot))
%     freezeColors
    h.c2 = colorbar;
%     h.c2=cbfreeze(h.c2);
    
    [tt.h tt.p tt.ci tt.stats]=ttest(pl.data_bc(:,:,i_dat)');
    t.data = abs(log10(tt.p));
    t.clim = [0 max(t.data(:))];
    t.pcriterion = abs(log10(0.05));
    if max(t.data(:))<t.pcriterion
        % temp.colormap = repmat([0.5 0.5 0.5],100,1);
        t.colormap = repmat(linspace(1,0.3,1000)',1,3);
    else
        t.border = ceil((t.pcriterion / max(t.data(:)))*1000);
        % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
        t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
    end
    h.s2=subplot(2,3,3+3*(i_dat-1));
    topoplot( t.data, TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',t.clim,'conv','on','colormap',t.colormap);
    title(sprintf('p-values %s %1.2f to %1.2f Hz;\n%1.0f to %1.0f ms',pl.titleadd{i_dat},pl.freq2plot, pl.time2plot))
%     freezeColors
    h.cb2 = colorbar;
    t.yticks = get(h.cb2,'YTick');
    set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
        'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')))
%     h.cb2=cbfreeze(h.cb2);
    
end

%% topoplot of only single topos
clear figs
% specific topos for effects
pl.time2plot=[-974.6 166]; pl.freq2plot=[8 12]; %vis. alpha
pl.time2plot=[-1482 447.3]; pl.freq2plot=[10 14]; %motor. alpha
pl.time2plot=[1400 2729]; pl.freq2plot=[10 14]; %motor. alpha
pl.time2plot=[-1553 252]; pl.freq2plot=[18 30]; %motor. beta
pl.time2plot=[384.8 3250]; pl.freq2plot=[18 30]; %motor. beta
pl.time2plot=[-1654 -287.1]; pl.freq2plot=[15 23]; % parietal beta
pl.time2plot=[377 1838]; pl.freq2plot=[15 23]; % parietal beta


% specific topos for effects [revisited
pl.time2plot=[-1037.11 244.14]; pl.freq2plot=[8 14]; %vis. alpha
pl.time2plot=[-1427.73 439.45]; pl.freq2plot=[8 14]; %motor. alpha
pl.time2plot=[1291 2728.52]; pl.freq2plot=[8 14]; %motor. alpha
pl.time2plot=[-1333.98 228.52]; pl.freq2plot=[15 30]; %motor. beta
pl.time2plot=[400.39 3250]; pl.freq2plot=[15 30]; %motor. beta


% pl.time2plot=[500 2500];
% pl.time2plot=[500 2000];
% pl.time2plot=[500 1000];
% pl.time2plot=[-1000 0];
% pl.time2plot=[-3000 -2750];
% pl.time2plot=[-500 300];
% pl.time2plot=[-2000 500];
% pl.time2plot=[-1500 0];
% pl.time2plot=[1500 2500];
% pl.time2plot=[-3500 3500];
[t.t t.ind1]=min(abs(TFA.time-pl.time2plot(1)));
[t.t t.ind2]=min(abs(TFA.time-pl.time2plot(2)));

% pl.freq2plot=[14.16667 14.16667];
% pl.freq2plot=[10 14];
% pl.freq2plot=[16 23];
% pl.freq2plot=[15 23];
% pl.freq2plot=[18 30];
% pl.freq2plot=[8 12];
% pl.freq2plot=[7 9];
[t.t t.ind3]=min(abs(TFA.frequency-pl.freq2plot(1)));
[t.t t.ind4]=min(abs(TFA.frequency-pl.freq2plot(2)));

pl.elec2plot = {'C3';'CP3'};
pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'};
pl.elec2plot = {};
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.subs2use = 1:numel(F.Subjects2Use);            

pl.data_raw=squeeze(nanmean(nanmean(TFA.data_induced(t.ind3:t.ind4,t.ind1:t.ind2,:,pl.subs2use),1),2));
pl.data_raw(:,:,2)=squeeze(nanmean(nanmean(TFA.data_evoked(t.ind3:t.ind4,t.ind1:t.ind2,:,pl.subs2use),1),2));
pl.data_bc=squeeze(nanmean(nanmean(TFA.data_induced_bc(t.ind3:t.ind4,t.ind1:t.ind2,:,pl.subs2use),1),2));
pl.data_bc(:,:,2)=squeeze(nanmean(nanmean(TFA.data_evoked_bc(t.ind3:t.ind4,t.ind1:t.ind2,:,pl.subs2use),1),2));

% colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))



figs{1}=figure;
set(gcf,'Position',[100 100 300 200],'PaperPositionMode','auto')
pl.data=squeeze(nanmean(pl.data_raw(:,:,1),2));
topoplot( pl.data, TFA.electrodes(1:64), ...
    'shading', 'flat', 'numcontour', 0, 'maplimits',[0 max(pl.data)], 'colormap', fake_parula,...
    'conv','on','whitebk','on','emarker2', {find(pl.elec2plot_i),'o','m',8});
% topoplot( pl.data, TFA.electrodes(1:64), ...
%     'shading', 'flat', 'numcontour', 0, 'maplimits',[0 max(pl.data)], 'colormap', flipud(cbrewer2('greys')),...
%     'conv','on','whitebk','on','emarker2', {find(pl.elec2plot_i),'o','m',8});
% topoplot( pl.data, TFA.electrodes(1:64), ...
%     'shading', 'flat', 'numcontour', 0, 'maplimits','minmax', 'colormap', cbrewer2('greys'),...
%     'conv','on','whitebk','on','emarker2', {find(pl.elec2plot_i),'o','m',8});
title(sprintf('raw; %1.2f to %1.2f Hz; %1.0f to %1.0f ms',pl.freq2plot, pl.time2plot))
%     freezeColors
h.cb = colorbar;
t.pos = get(h.cb,'Position');
set(h.cb,'Position',[t.pos(1)+0.08 t.pos(2)+(1/6)*t.pos(4) t.pos(3) t.pos(4)*2/3 ])
% set(h.cb,'Position',[t.pos(1)+0.1 t.pos(2) t.pos(3) t.pos(4)])
%     h.c1=cbfreeze(h.c1);


figs{2}=figure;
pl.data=squeeze(nanmean(pl.data_bc(:,:,1),2));
set(gcf,'Position',[100 100 300 200],'PaperPositionMode','auto')
topoplot( pl.data, TFA.electrodes(1:64),...
    'shading', 'flat', 'numcontour', 0,'maplimits','absmax', 'colormap', flipud(cbrewer2('RdBu')),...
    'whitebk','on','conv','on','emarker2', {find(pl.elec2plot_i),'o','m',8});
% topoplot( pl.data, TFA.electrodes(1:64),...
%     'shading', 'flat', 'numcontour', 0,'maplimits',[-10 10], 'colormap', flipud(cbrewer2('RdBu')),...
%     'whitebk','on','conv','on','emarker2', {find(pl.elec2plot_i),'o','m',8});
title(sprintf('bc; %1.2f to %1.2f Hz; %1.0f to %1.0f ms',pl.freq2plot, pl.time2plot))
%     freezeColors
h.cb = colorbar;
t.pos = get(h.cb,'Position');
set(h.cb,'Position',[t.pos(1)+0.08 t.pos(2)+(1/6)*t.pos(4) t.pos(3) t.pos(4)*2/3 ])


sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
sav.pathout = 'C:\Users\EEG\Documents\MATLAB\christopher\SSVEP_volmov\';
sav.filenames = {'TOPO_raw_effects';'TOPO_bc_effects'};
for i_fig = 1:2
    print(figs{i_fig}, fullfile(sav.pathout,...
        sprintf('%s_%1.0f_%1.0fHz_%1.0f_%1.0fms',sav.filenames{i_fig},pl.freq2plot,pl.time2plot)),'-dpng','-r300')
%     print(figs{i_fig}, fullfile(sav.pathout,...
%         sprintf('%s_%1.0f_%1.0fHz_%1.0f_%1.0fms',sav.filenames{i_fig},pl.freq2plot,pl.time2plot)),'-djpeg','-r300')
    saveas(figs{i_fig},fullfile(sav.pathout,...
        sprintf('%s_%1.0f_%1.0fHz_%1.0f_%1.0fms',sav.filenames{i_fig},pl.freq2plot,pl.time2plot)),'fig')
    print(figs{i_fig},fullfile(sav.pathout,...
        sprintf('%s_%1.0f_%1.0fHz_%1.0f_%1.0fms',sav.filenames{i_fig},pl.freq2plot,pl.time2plot)),'-depsc2', '-painters','-r300')
end

%% topoplot differences...rudimentary
clear figs
% time windows of to be subtracted data
pl.time2plot = {[500 2000];[500 2000]};
% pl.time2plot = {[-1000 0];[-1000 0]};
pl.time2plot_i=cell2mat(cellfun(@(x) dsearchn(TFA.time',x'),pl.time2plot,'UniformOutput',false)');

pl.freq2plot={[16 23];[18 30]};
% pl.freq2plot={[8 12];[10 14]};
pl.freq2plot_i=cell2mat(cellfun(@(x) dsearchn(TFA.frequency',x'),pl.freq2plot,'UniformOutput',false)');

pl.elec2plot = {'C3';'CP3'};
pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'};
pl.elec2plot = {};
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.subs2use = 1:numel(F.Subjects2Use);            

pl.data_raw=squeeze(nanmean(nanmean(TFA.data_induced(pl.freq2plot_i(1,1):pl.freq2plot_i(2,1),pl.time2plot_i(1,1):pl.time2plot_i(2,1),:,pl.subs2use),1),2))-...
    squeeze(nanmean(nanmean(TFA.data_induced(pl.freq2plot_i(1,2):pl.freq2plot_i(2,2),pl.time2plot_i(1,2):pl.time2plot_i(2,2),:,pl.subs2use),1),2));
pl.data_raw(:,:,2)=squeeze(nanmean(nanmean(TFA.data_evoked(pl.freq2plot_i(1,1):pl.freq2plot_i(2,1),pl.time2plot_i(1,1):pl.time2plot_i(2,1),:,pl.subs2use),1),2))-...
    squeeze(nanmean(nanmean(TFA.data_evoked(pl.freq2plot_i(1,2):pl.freq2plot_i(2,2),pl.time2plot_i(1,2):pl.time2plot_i(2,2),:,pl.subs2use),1),2));
pl.data_bc=squeeze(nanmean(nanmean(TFA.data_induced_bc(pl.freq2plot_i(1,1):pl.freq2plot_i(2,1),pl.time2plot_i(1,1):pl.time2plot_i(2,1),:,pl.subs2use),1),2))-...
    squeeze(nanmean(nanmean(TFA.data_induced_bc(pl.freq2plot_i(1,2):pl.freq2plot_i(2,2),pl.time2plot_i(1,2):pl.time2plot_i(2,2),:,pl.subs2use),1),2));
pl.data_bc(:,:,2)=squeeze(nanmean(nanmean(TFA.data_evoked_bc(pl.freq2plot_i(1,1):pl.freq2plot_i(2,1),pl.time2plot_i(1,1):pl.time2plot_i(2,1),:,pl.subs2use),1),2))-...
    squeeze(nanmean(nanmean(TFA.data_evoked_bc(pl.freq2plot_i(1,2):pl.freq2plot_i(2,2),pl.time2plot_i(1,2):pl.time2plot_i(2,2),:,pl.subs2use),1),2));
% colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))



figs{1}=figure;
set(gcf,'Position',[100 100 300 200],'PaperPositionMode','auto')
pl.data=squeeze(nanmean(pl.data_raw(:,:,1),2));
h.tp=topoplot( pl.data, TFA.electrodes(1:64),...
    'shading', 'flat', 'numcontour', 0,'maplimits','absmax', 'colormap', flipud(cbrewer2('RdBu')),...
    'whitebk','on','conv','on','emarker2', {find(pl.elec2plot_i),'o','m',8});
t.pos = get(gca,'Position');
set(gca,'Position',[t.pos(1) t.pos(2)-0.1 t.pos(3) t.pos(4) ])
title(sprintf('raw diff; [%1.1f:%1.1f]-[%1.1f:%1.1f] Hz\n [%1.0f:%1.0f]-[%1.0f:%1.0f] ms',...
    pl.freq2plot{1},pl.freq2plot{2}, pl.time2plot{1},pl.time2plot{1}))
%     freezeColors
h.cb = colorbar;
t.pos = get(h.cb,'Position');
set(h.cb,'Position',[t.pos(1)+0.08 t.pos(2)+(1/6)*t.pos(4) t.pos(3) t.pos(4)*2/3 ])
% set(h.cb,'Position',[t.pos(1)+0.1 t.pos(2) t.pos(3) t.pos(4)])
%     h.c1=cbfreeze(h.c1);


figs{2}=figure;
pl.data=squeeze(nanmean(pl.data_bc(:,:,1),2));
set(gcf,'Position',[100 100 300 200],'PaperPositionMode','auto')
h.tp=topoplot( pl.data, TFA.electrodes(1:64),...
    'shading', 'flat', 'numcontour', 0,'maplimits','absmax', 'colormap', flipud(cbrewer2('RdBu')),...
    'whitebk','on','conv','on','emarker2', {find(pl.elec2plot_i),'o','m',8});
h.tp=topoplot( pl.data, TFA.electrodes(1:64),...
    'shading', 'flat', 'numcontour', 0,'maplimits',[-5 5], 'colormap', flipud(cbrewer2('RdBu')),...
    'whitebk','on','conv','on','emarker2', {find(pl.elec2plot_i),'o','m',8});
t.pos = get(gca,'Position');
set(gca,'Position',[t.pos(1) t.pos(2)-0.1 t.pos(3) t.pos(4) ])
title(sprintf('bc diff; [%1.1f:%1.1f]-[%1.1f:%1.1f] Hz\n [%1.0f:%1.0f]-[%1.0f:%1.0f] ms',...
    pl.freq2plot{1},pl.freq2plot{2}, pl.time2plot{1},pl.time2plot{1}))
%     freezeColors
h.cb = colorbar;
t.pos = get(h.cb,'Position');
set(h.cb,'Position',[t.pos(1)+0.08 t.pos(2)+(1/6)*t.pos(4) t.pos(3) t.pos(4)*2/3 ])


sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
sav.filenames = {'TOPO_raw';'TOPO_bc'};
% for i_fig = 1:2
% %     print(figs{i_fig}, fullfile(sav.pathout,...
% %         sprintf('%s_%1.0f_%1.0fHz_%1.0f_%1.0fms',sav.filenames{i_fig},pl.freq2plot,pl.time2plot)),'-dpng','-r300')
%     print(figs{i_fig}, fullfile(sav.pathout,...
%         sprintf('%s_%1.0f_%1.0fHz_%1.0f_%1.0fms',sav.filenames{i_fig},pl.freq2plot,pl.time2plot)),'-djpeg','-r300')
%     saveas(figs{i_fig},fullfile(sav.pathout,...
%         sprintf('%s_%1.0f_%1.0fHz_%1.0f_%1.0fms',sav.filenames{i_fig},pl.freq2plot,pl.time2plot)),'fig')
%     print(figs{i_fig},fullfile(sav.pathout,...
%         sprintf('%s_%1.0f_%1.0fHz_%1.0f_%1.0fms',sav.filenames{i_fig},pl.freq2plot,pl.time2plot)),'-depsc2', '-painters','-r300')
% end

%% plot RESS.topography

% pl.data = mean(TFA.RESS_map(:,:),2);
% pl.data = mean(TFA.RESS_map(:,:).*-1,2);
pl.data = mean(TFA.RESS_map_signchan(:,:),2);
figs{1}=figure;
set(gcf,'Position',[100 100 300 200],'PaperPositionMode','auto')
topoplot( pl.data, TFA.electrodes(1:64),'shading', 'flat', 'numcontour', 0,'maplimits','absmax',...
    'whitebk','on','conv','on','colormap',fake_parula);
h.cb = colorbar;
t.pos = get(h.cb,'Position');
set(h.cb,'Position',[t.pos(1)+0.08 t.pos(2)+(1/6)*t.pos(4) t.pos(3) t.pos(4)*2/3 ])

sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
sav.filenames = {'TOPO_RESS'};
for i_fig = 1:1
%     print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-dpng','-r300')
    print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-djpeg','-r300')
    saveas(figs{i_fig},fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'fig')
    print(figs{i_fig},fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-depsc2', '-painters','-r300')
end

%% plot electrode head
% pl.elec2plot = {'C3';'CP3'};
pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'};


figure;
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
% topoplot(find(pl.elec2plot_i),TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on');
topoplot([],TFA.electrodes(1:64),'whitebk','on','style','blank','electrodes', 'on', 'emarker2', {find(pl.elec2plot_i),'o','r',8});

%% plot topographies across time


% define parameters
t.p_time = [100 100 -3000 3000]; % width step min max
t.p_time = [100 100 -3000 200]; % width step min max
t.p_time = [100 100 -3000 500]; % width step min max
% t.p_time = [100 100 -3000 3200]; % width step min max
t.p_time = [100 100 -2500 3000]; % width step min max
t.posScale = 1.1;

% pl.freq2plot=[14.16667 14.16667];
% pl.freq2plot=[10 14];
% pl.freq2plot=[16 23];
% pl.freq2plot=[15 19];
% pl.freq2plot=[8 12];
% pl.freq2plot=[15 30];
% pl.freq2plot=[18 30];
pl.freq2plot=[2 4];
% pl.freq2plot=[3 8];
[t.t t.ind3]=min(abs(TFA.frequency-pl.freq2plot(1)));
[t.t t.ind4]=min(abs(TFA.frequency-pl.freq2plot(2)));

pl.data=squeeze(nanmean(nanmean(TFA.data_induced(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
pl.data(:,:,2)=squeeze(nanmean(nanmean(TFA.data_evoked(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
pl.data(:,:,3)=squeeze(nanmean(nanmean(TFA.data_induced_bc(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
pl.data(:,:,4)=squeeze(nanmean(nanmean(TFA.data_evoked_bc(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
t.conlabel={'raw induced';'raw evoked';'normalized induced';'normalized evoked'};
t.con_lims = [1 1 0 0]; % 1 = 0 to max; % 0 = minmax

t.time=[];
t.timedot=[];
for i_st = 1:floor((t.p_time(4)-t.p_time(1)-t.p_time(3))/t.p_time(2))+1
    t.time(i_st,:)=t.p_time(3)+t.p_time(2)*(i_st-1)+[0 t.p_time(1)];
    [t.t t.timedot(i_st,1)]=min(abs(TFA.time-t.time(end,1)));
    [t.t t.timedot(i_st,2)]=min(abs(TFA.time-t.time(end,2)));
end
% t.t = get(0,'MonitorPositions');
% t.row = round(sqrt((size(t.time,1)+2)/(1/(t.t(1,4)/t.t(1,3)))));
% t.col = ceil(sqrt((size(t.time,1)+2)/(t.t(1,4)/t.t(1,3))));
t.t = [1920 1080];
t.col = ceil(sqrt((size(t.time,1)+2)/(1/(t.t(1)/t.t(2)))));
t.row = ceil(sqrt((size(t.time,1)+2)/(t.t(1)/t.t(2))));
% create plotdata

plotdata=[];
for i_pl = 1:size(t.time,1)
    plotdata(:,:,i_pl)=squeeze(mean(pl.data(t.timedot(i_pl,1):t.timedot(i_pl,2),:,:),1))';
end


for i_fig = 1:size(plotdata,1)
    h.fig(i_fig)=figure;
    if t.con_lims(i_fig) == 1
        t.lims = [0 1]*max(max(abs(plotdata(i_fig,:,:))));
    else
        t.lims = [-1 1]*max(max(abs(plotdata(i_fig,:,:))));
    end
    for i_spl = 1:size(plotdata,3)
        h.sp(i_spl)=subplot(t.row,t.col,i_spl);
        if i_fig < 3
            topoplot( plotdata(i_fig,:,i_spl), TFA.electrodes(1:64), ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims,'electrodes','off','colormap',fake_parula);
        else
            topoplot( plotdata(i_fig,:,i_spl), TFA.electrodes(1:64), ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims,'electrodes','off','colormap',flipud(cbrewer2('RdBu')));
        end
        title(sprintf('[%1.0f %1.0f]',t.time(i_spl,1),t.time(i_spl,2)),'FontSize',8)
        t.pos = get(h.sp(i_spl),'Position');
        set(h.sp(i_spl),'Position',[t.pos(1:2)-(t.pos(3:4).*((t.posScale-1)/2)) t.pos(3:4).*t.posScale])
    end
    h.sp(i_spl+1)=subplot(t.row,t.col,i_spl+2);
    topoplot( [], TFA.electrodes(1:64),  ...
        'style','blank');
    title(sprintf('%s\n[%1.2f %1.2f]Hz',t.conlabel{i_fig},pl.freq2plot),'FontSize',8)
    t.pos = get(h.sp(i_spl),'Position');
        set(h.sp(i_spl),'Position',[t.pos(1:2)-(t.pos(3:4).*((t.posScale-1)/2)) t.pos(3:4).*t.posScale])
    
    t.pos2 = get(h.sp(i_spl+1),'Position');
    t.pos3 = get(h.sp(i_spl+1),'OuterPosition');
    h.a1 = axes('position',[t.pos3(1) t.pos2(2) t.pos3(3) t.pos2(4)],'Visible','off');
    if i_fig < 3
        colormap(gca, fake_parula)
    else
        colormap(gca,flipud(cbrewer2('RdBu')))
    end
    caxis(t.lims);
    h.c3 = colorbar();
    t.pos4 = get(h.c3,'Position');
    set(h.c3,'Position',[t.pos4(1)+0.065 t.pos4(2)+(t.pos4(4)*(1/6)) t.pos4(3)/2 t.pos4(4)*(2/3)])
end
axcopy(h.fig(i_fig))

%% plot topographies across time [movie]


% define parameters
t.p_time = [100 -3000 3000]; % width step min max

% pl.freq2plot=[14.16667 14.16667];
pl.freq2plot=[10 14];
% pl.freq2plot=[16 25];
% pl.freq2plot=[8 12];
% pl.freq2plot=[15 30];
% pl.freq2plot=[18 30];
[t.t t.ind3]=min(abs(TFA.frequency-pl.freq2plot(1)));
[t.t t.ind4]=min(abs(TFA.frequency-pl.freq2plot(2)));

pl.data=squeeze(nanmean(nanmean(TFA.data_induced(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
pl.data(:,:,2)=squeeze(nanmean(nanmean(TFA.data_evoked(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
pl.data(:,:,3)=squeeze(nanmean(nanmean(TFA.data_induced_bc(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
pl.data(:,:,4)=squeeze(nanmean(nanmean(TFA.data_evoked_bc(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
t.conlabel={'raw induced';'raw evoked';'normalized induced';'normalized evoked'};
t.con_lims = [1 1 0 0]; % 1 = 0 to max; % 0 = minmax


t.timedot=[];
t.t_lims = dsearchn(TFA.time', t.p_time(2:3)');
t.t_steps = t.t_lims(1):t.t_lims(2);
for i_st = 1:numel(t.t_steps)
    % index time windows to be averaged
    t.timedot(i_st,:)=dsearchn(TFA.time',TFA.time(t.t_steps(i_st))+[-1; 1]*t.p_time(1)/2);
end

% create data
plotdata=nan(numel(t.t_lims(1):t.t_lims(2)),size(pl.data,2),size(pl.data,3));
for i_pl = 1:size(t.timedot,1)
    plotdata(i_pl,:,:)=squeeze(mean(pl.data(t.timedot(i_pl,1):t.timedot(i_pl,2),:,:),1));
end

% create movie
for i_dat = 3
    figure;
    if t.con_lims(i_dat) == 1
        pl.clims = [0 1]*max(max(abs(plotdata(:,:,i_dat))));
    else
        pl.clims = [-1 1]*max(max(abs(plotdata(:,:,i_dat))));
    end
    if i_dat < 3
        [mov.Movie,mov.Colormap] = eegmovie(plotdata(:,:,i_dat)',TFA.srate,TFA.electrodes,...
            'topoplotopt',{'maplimits',pl.clims,'numcontour', 0, 'conv','on', 'shading', 'flat','colormap',fake_parula},...
            'time','on','startsec',-3);
    else
        [mov.Movie,mov.Colormap] = eegmovie(plotdata(:,:,i_dat)',TFA.srate,TFA.electrodes,...
            'topoplotopt',{'maplimits',pl.clims,'numcontour', 0, 'conv','on', 'shading', 'flat','colormap',flipud(cbrewer2('RdBu'))},...
            'time','on','startsec',-3);
    end
end
figure; seemovie(mov.Movie,-5,mov.Colormap);


%% plot topographies across time [save single images]
pl.dat2plot = 3; %1=induced raw; 2=evoked raw; 3=induced bc; 4=evoked bc

% sav.pathout = 'C:\Users\psy05cvd\Pictures\gifski\SSVEP_volmov\upper_alpha_topo';
% sav.pathout = 'C:\Users\psy05cvd\Pictures\gifski\SSVEP_volmov\lower_alpha_topo';
sav.pathout = 'C:\Users\psy05cvd\Pictures\gifski\SSVEP_volmov\lower_beta_topo';
% sav.pathout = 'C:\Users\psy05cvd\Pictures\gifski\SSVEP_volmov\higher_beta_topo';


% define parameters
t.p_time = [100 100 -3000 3000]; % width step min max
t.p_time = [100 100 -3000 200]; % width step min max
t.p_time = [100 100 -3000 500]; % width step min max
% t.p_time = [100 100 -3000 3200]; % width step min max
t.p_time = [100 50 -3000 3000]; % width step min max
t.posScale = 1.1;

% pl.freq2plot=[14.16667 14.16667];
% pl.freq2plot=[10 14];
% pl.freq2plot=[16 23];
pl.freq2plot=[15 23];
% pl.freq2plot=[15 19];
% pl.freq2plot=[8 12];
% pl.freq2plot=[15 30];
% pl.freq2plot=[18 30];
% pl.freq2plot=[2 4];
% pl.freq2plot=[3 8];
[t.t t.ind3]=min(abs(TFA.frequency-pl.freq2plot(1)));
[t.t t.ind4]=min(abs(TFA.frequency-pl.freq2plot(2)));

pl.data=squeeze(nanmean(nanmean(TFA.data_induced(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
pl.data(:,:,2)=squeeze(nanmean(nanmean(TFA.data_evoked(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
pl.data(:,:,3)=squeeze(nanmean(nanmean(TFA.data_induced_bc(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
pl.data(:,:,4)=squeeze(nanmean(nanmean(TFA.data_evoked_bc(t.ind3:t.ind4,:,:,pl.subs2use),1),4));
t.conlabel={'raw_induced';'raw_evoked';'normalized_induced';'normalized_evoked'};
t.con_lims = [1 1 0 0]; % 1 = 0 to max; % 0 = minmax

t.time=[];
t.timedot=[];
for i_st = 1:floor((t.p_time(4)-t.p_time(1)-t.p_time(3))/t.p_time(2))+1
    t.time(i_st,:)=t.p_time(3)+t.p_time(2)*(i_st-1)+[0 t.p_time(1)];
    [t.t t.timedot(i_st,1)]=min(abs(TFA.time-t.time(end,1)));
    [t.t t.timedot(i_st,2)]=min(abs(TFA.time-t.time(end,2)));
end
% t.t = get(0,'MonitorPositions');
% t.row = round(sqrt((size(t.time,1)+2)/(1/(t.t(1,4)/t.t(1,3)))));
% t.col = ceil(sqrt((size(t.time,1)+2)/(t.t(1,4)/t.t(1,3))));
t.t = [1920 1080];
t.col = ceil(sqrt((size(t.time,1)+2)/(1/(t.t(1)/t.t(2)))));
t.row = ceil(sqrt((size(t.time,1)+2)/(t.t(1)/t.t(2))));
% create plotdata

plotdata=[];
for i_pl = 1:size(t.time,1)
    plotdata(:,:,i_pl)=squeeze(mean(pl.data(t.timedot(i_pl,1):t.timedot(i_pl,2),:,:),1))';
end


fig1=figure;
set(gcf,'Position',[100 100 200 230],'PaperPositionMode','auto')

if t.con_lims(pl.dat2plot) == 1
    t.lims = [0 1]*max(max(abs(plotdata(pl.dat2plot,:,:))));
else
    t.lims = [-1 1]*max(max(abs(plotdata(pl.dat2plot,:,:))));
    t.lims = [-1 1]*10;
end
for i_spl = 1:size(plotdata,3)
    clf(fig1,'reset')
    h.sp1=subplot(4,1,1:3);
    if pl.dat2plot < 3
        topoplot( plotdata(pl.dat2plot,:,i_spl), TFA.electrodes(1:64), ...
            'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims,'electrodes','off',...
            'colormap',fake_parula,'whitebk','on');
    else
        topoplot( plotdata(pl.dat2plot,:,i_spl), TFA.electrodes(1:64), ...
            'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims,'electrodes','off',...
            'colormap',flipud(cbrewer2('RdBu')),'whitebk','on');
    end
    h.cb = colorbar;
    
    %
    h.sp2=subplot(4,1,4);
    plot([min(t.time(:)) max(t.time(:))], [0 0],'Color','k','LineWidth',1)
    hold on
    line([0;0],[-1; 1],'Color',[0 0 0],'LineWidth', 1)
    line([-3000 -2000 -1000 1000 2000 3000;-3000 -2000 -1000 1000 2000 3000],...
        [-0.5 -0.5 -0.5 -0.5 -0.5 -0.5; 0.5  0.5  0.5  0.5  0.5  0.5],'Color',[0 0 0],'LineWidth', 1)
    line([t.time(i_spl,1) t.time(i_spl,2) t.time(i_spl,1) t.time(i_spl,1); t.time(i_spl,2) t.time(i_spl,2) t.time(i_spl,2) t.time(i_spl,1)],...
        [-1 -1 1 -1; -1 1 1 1],'Color',[0.8 0 0],'LineWidth', 2)
    xlim([min(t.time(:)) max(t.time(:))])
    ylim([-1.5 1.5])
    %     set(gca,'YColor','none','XAxisLocation','Origin','Box','Off')
    set(gca,'YColor','none','XColor',[1 1 1],'Box','Off','XTick',[])
    xlabel(sprintf('[%1.0f %1.0f]ms',t.time(i_spl,:)),'Color','k')
    title(sprintf('[%1.0f %1.0f]ms | [%1.0f %1.0f]Hz',[min(t.time(:)) max(t.time(:))], pl.freq2plot))
    
    t.pos = get(h.sp1,'Position');
%     set(h.sp1,'Position',[t.pos(1:2)-(t.pos(3:4).*((t.posScale-1)/2)) t.pos(3:4).*t.posScale])
%     set(h.sp1,'Position',[t.pos(1:2)-(t.pos(3:4).*((t.posScale-1)*2)) t.pos(3:4).*t.posScale])
    set(h.sp1,'Position',[t.pos(1)-(t.pos(3).*((t.posScale-1)*2.5)) t.pos(2)-(t.pos(4).*((t.posScale-1)/2)) t.pos(3:4).*t.posScale])
    t.pos2 = get(h.cb,'Position');
    set(h.cb,'Position',[t.pos2(1) t.pos2(2)+(1/6)*t.pos2(4) t.pos2(3) t.pos2(4)*2/3 ])
    
    % save file
    print(fig1, ...
        fullfile(sav.pathout,sprintf('%s_%1.0f_%1.0fHz_%04.0f',t.conlabel{pl.dat2plot},pl.freq2plot,i_spl)),...
        '-dpng','-r300')
end


%% calculate correlations between amplitudes for frequencies and timewindows of interest
pl.parameters = {...
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

pl.subs2use = 1:numel(F.Subjects2Use);
% pl.subs2use = [1:14 16 :20];
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 13 14 17 18 19 20]; % subjects discarded due to low trial number version 1
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20]; % subjects discarded due to low trial number version 2
% 
% pl.subs2use = [1 2 3 4 5 7 8 10 11 12 14 17 18 19 20 24 25 26 27]; % subjects discarded due to low trial number version 2


h.pl = [];
for i_pl = 1:size(pl.parameters,1)
    figure;
    % index frequencies
    t.ind1 = dsearchn(TFA.frequency',pl.parameters{i_pl,1}');
    t.ind2 = dsearchn(TFA.frequency',pl.parameters{i_pl,5}');
   
    % index time
    t.ind3 = dsearchn(TFA.time',pl.parameters{i_pl,2}');
    t.ind4 = dsearchn(TFA.time',pl.parameters{i_pl,6}');
    
    % index electrodes
    pl.elec2plot_i1=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.parameters{i_pl,3}, 'UniformOutput',false)),1));
    pl.elec2plot_i2=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.parameters{i_pl,7}, 'UniformOutput',false)),1));
    
    % get xdata
    pl.xdata = squeeze(mean(mean(mean(TFA.data_induced_bc(t.ind1(1):t.ind1(2),t.ind3(1):t.ind3(2),pl.elec2plot_i1,:),1),2),3));
    
    % get ydata
    pl.ydata = squeeze(mean(mean(mean(TFA.data_induced_bc(t.ind2(1):t.ind2(2),t.ind4(1):t.ind4(2),pl.elec2plot_i2,:),1),2),3));
    
    % plotting
    
    h.sc1=scatter(pl.xdata , pl.ydata ,'ko');
    hold on
    
    % xlims
    xlim(1.1*[-1 1]*max(abs([pl.xdata;pl.ydata])));
    ylim(1.1*[-1 1]*max(abs([pl.xdata;pl.ydata])));
    % correlation
    [r.R,r.P]=corrcoef(pl.xdata,pl.ydata);
    text(mean(get(gca,'xlim')),max(get(gca,'ylim')),...
        sprintf('R=%1.3f p=%1.4f\n',r.R(2),r.P(2)), ...
        'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8);
    % add fitted line
    r.fit = polyfit(pl.xdata, pl.ydata,1);
    t.x=linspace(min(get(gca,'xlim')),max(get(gca,'xlim')),100);
    plot(t.x,r.fit(1)*t.x+r.fit(2),'-','Color',[0.6,0.6,0.6])
    
    t.pos=get(gca,'Position');
    
    % title and labels
    title(sprintf(...
        '%s | %1.1f to %1.1f Hz | %1.0f to %1.0f ms | %s X\n%s | %1.1f to %1.1f Hz | %1.0f to %1.0f ms | %s',...
        pl.parameters{i_pl,4}, pl.parameters{i_pl,1}, pl.parameters{i_pl,2}, vararg2str(pl.parameters{i_pl,3}),...
        pl.parameters{i_pl,8}, pl.parameters{i_pl,5}, pl.parameters{i_pl,6}, vararg2str(pl.parameters{i_pl,7})...
        ),'FontSize',8)
    xlabel(pl.parameters{i_pl,4},'FontSize',8)
    ylabel(pl.parameters{i_pl,8},'FontSize',8)
    
    grid on
    box on
    
    
end

%% calculate time-lagged correlations between signals
clear figs h


pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
    [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [1 0 0] 2 'vis alpha' 'noRESS';...
    [10 14] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [0 0.5 0] 2 'mot alpha' 'noRESS';...
    [18 30] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [0 0.5 0.9] 2 'mot beta' 'noRESS';...
    [15 23] {'CP1';'CPz';'CP2';'P1';'Pz';'P2'} 'data_induced_bc' [0.8 0 0.8] 2 'parietal beta' 'noRESS';...
    };

pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
    [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [255 151 0]./255 2 'vis alpha' 'noRESS';...
    [10 14] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [33 92 150]./255 2 'mot alpha' 'noRESS';...
    [18 30] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [170 99 57]./255 2 'mot beta' 'noRESS';...
    [15 23] {'CP1';'CPz';'CP2';'P1';'Pz';'P2'} 'data_induced_bc' [23 150 118]./255 2 'parietal beta' 'noRESS';...
    [2 4] {'C3';'CP3';'C1';'CP1'} 'data_induced_bc' [0.4 0.4 0.4] 2 'mot delta' 'noRESS';...
    };


% this one for the poster figures
% pl.parameters = {...
%     [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
%     [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [255 151 0]./255 2 'vis alpha' 'noRESS';...
%     [10 14] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [33 92 150]./255 2 'mot alpha' 'noRESS';...
%     [18 30] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [170 99 57]./255 2 'mot beta' 'noRESS';...
%     [15 23] {'CP1';'CPz';'CP2';'P1';'Pz';'P2'} 'data_induced_bc' [23 150 118]./255 2 'parietal beta' 'noRESS';...
%     };

pl.corrsigs = [1 2; 1 3; 2 3];
pl.corrsigs = [2 3; 3 4];
pl.corrtlim = [-1000 1000]; % time limits of shifted correlation analysis


pl.subs2use = 1:numel(F.Subjects2Use);

pl.xlim = [-3500 3500];

pl.permut_n = 10000;
% pl.permut_n = 100;


clear pl.tdata2 h.sp h.p h.pm h.pl
% do the analyisis
for i_fig = 1:size(pl.corrsigs,1)
    figs{i_fig}=figure;
    set(gcf,'Position',[100 100 600 950],'PaperPositionMode','auto')
    
    
    % first dataset
    % index frequencies
    [t.t t.ind1]=min(abs(TFA.frequency-pl.parameters{pl.corrsigs(i_fig,1),1}(1)));
    [t.t t.ind2]=min(abs(TFA.frequency-pl.parameters{pl.corrsigs(i_fig,1),1}(2)));
    
    % index time
    [t.t t.ind3]=min(abs(TFA.time-pl.xlim(1)));
    [t.t t.ind4]=min(abs(TFA.time-pl.xlim(2)));
    
    % index electrodes
    switch  pl.parameters{pl.corrsigs(i_fig,1),7}
        case 'noRESS'
            pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.parameters{pl.corrsigs(i_fig,1),2}, 'UniformOutput',false)),1));
            com=sprintf('pl.data1=squeeze(nanmean(nanmean(TFA.%s(t.ind1:t.ind2,:,pl.elec2plot_i,pl.subs2use),1),3));',pl.parameters{pl.corrsigs(i_fig,1),3});
        case 'RESS'
            pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes_RESS.labels},x), pl.parameters{pl.corrsigs(i_fig,1),2}, 'UniformOutput',false)),1));
            com=sprintf('pl.data1=squeeze(nanmean(TFA.%s(t.ind1:t.ind2,:,pl.subs2use),1));',pl.parameters{pl.corrsigs(i_fig,1),3});
        
    end
    eval(com)
    [tt.h1 tt.p1 tt.ci1 tt.stats1]=ttest(pl.data1(t.ind3:t.ind4,:)');
    
    % second dataset
    % index frequencies
    [t.t t.ind1]=min(abs(TFA.frequency-pl.parameters{pl.corrsigs(i_fig,2),1}(1)));
    [t.t t.ind2]=min(abs(TFA.frequency-pl.parameters{pl.corrsigs(i_fig,2),1}(2)));
    
    % index time
    [t.t t.ind3]=min(abs(TFA.time-pl.xlim(1)));
    [t.t t.ind4]=min(abs(TFA.time-pl.xlim(2)));
    
    % index electrodes
    switch  pl.parameters{pl.corrsigs(i_fig,2),7}
        case 'noRESS'
            pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.parameters{pl.corrsigs(i_fig,2),2}, 'UniformOutput',false)),1));
            com=sprintf('pl.data2=squeeze(nanmean(nanmean(TFA.%s(t.ind1:t.ind2,:,pl.elec2plot_i,pl.subs2use),1),3));',pl.parameters{pl.corrsigs(i_fig,2),3});
        case 'RESS'
            pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes_RESS.labels},x), pl.parameters{pl.corrsigs(i_fig,2),2}, 'UniformOutput',false)),1));
            com=sprintf('pl.data2=squeeze(nanmean(TFA.%s(t.ind1:t.ind2,:,pl.subs2use),1));',pl.parameters{pl.corrsigs(i_fig,2),3});
        
    end
    eval(com)
    [tt.h2 tt.p2 tt.ci2 tt.stats2]=ttest(pl.data2(t.ind3:t.ind4,:)');
    
    
    %% plot data
    
    % data1
    % actual plot with boundedline
    pl.mdata = mean(pl.data1(t.ind3:t.ind4,:),2)';
    pl.ddata = tt.ci1 - repmat(mean(pl.data1(t.ind3:t.ind4,:),2)',2,1);
    %pl.ddata = tt.ci;
    subplot(5,2,1)
    [h.l(1), h.p(1)] = boundedline(TFA.time(t.ind3:t.ind4),pl.mdata,pl.ddata(1,:));
    set(h.l(1),'Color',pl.parameters{pl.corrsigs(i_fig,1),4},'LineWidth',pl.parameters{pl.corrsigs(i_fig,1),5});
    set(h.p(1),'FaceColor',pl.parameters{pl.corrsigs(i_fig,1),4},'FaceAlpha',0.25)
    hold on
    
    % plot single subject data
%     pl.mdata = mean(pl.data,2)';
%     h.p(i_pl,:)=plot(TFA.time(t.ind3:t.ind4),pl.data,'Color',[pl.parameters{i_pl,4} 0.3],'LineWidth',0.2);
%     hold on
%     h.pm(i_pl)=plot(TFA.time(t.ind3:t.ind4),pl.mdata,'Color',[pl.parameters{i_pl,4} 1],'LineWidth',2);
    
    % significant effects?
    pl.tdata2=nan(size(pl.mdata));
    pl.tdata2(tt.p1<.05)=1;
    set(gca,'ylim',[-1.1 (1+0.15+0.05*4)]...
        .*max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false))))
    pl.maxdata = max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false)));
    pl.y = 1*pl.maxdata*(1+0.05+0.05*2);
    plot(TFA.time(t.ind3:t.ind4),pl.tdata2.*pl.y,'Color',[pl.parameters{pl.corrsigs(i_fig,1),4} 1],'LineWidth',6)
    
    
    % graphics
    xlim(pl.xlim)
    title(pl.parameters{pl.corrsigs(i_fig,1),6},'FontSize', 8)
    grid on
    box on
    hline(0,'k')
    vline(TFA.baseline,'k')
    xlabel('time in ms')
    ylabel('amplitude modulation in %')
    
    
    % data2
    % actual plot with boundedline
    pl.mdata = mean(pl.data2(t.ind3:t.ind4,:),2)';
    pl.ddata = tt.ci2 - repmat(mean(pl.data2(t.ind3:t.ind4,:),2)',2,1);
    %pl.ddata = tt.ci;
    subplot(5,2,2)
    [h.l(1), h.p(1)] = boundedline(TFA.time(t.ind3:t.ind4),pl.mdata,pl.ddata(1,:));
    set(h.l(1),'Color',pl.parameters{pl.corrsigs(i_fig,2),4},'LineWidth',pl.parameters{pl.corrsigs(i_fig,2),5});
    set(h.p(1),'FaceColor',pl.parameters{pl.corrsigs(i_fig,2),4},'FaceAlpha',0.25)
    hold on
    
    % plot single subject data
%     pl.mdata = mean(pl.data,2)';
%     h.p(i_pl,:)=plot(TFA.time(t.ind3:t.ind4),pl.data,'Color',[pl.parameters{i_pl,4} 0.3],'LineWidth',0.2);
%     hold on
%     h.pm(i_pl)=plot(TFA.time(t.ind3:t.ind4),pl.mdata,'Color',[pl.parameters{i_pl,4} 1],'LineWidth',2);
    
    % significant effects?
    pl.tdata2=nan(size(pl.mdata));
    pl.tdata2(tt.p2<.05)=1;
    set(gca,'ylim',[-1.1 (1+0.15+0.05*4)]...
        .*max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false))))
    pl.maxdata = max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false)));
    pl.y = 1*pl.maxdata*(1+0.05+0.05*2);
    plot(TFA.time(t.ind3:t.ind4),pl.tdata2.*pl.y,'Color',[pl.parameters{pl.corrsigs(i_fig,2),4} 1],'LineWidth',6)
    
    
    % graphics
    xlim(pl.xlim)
    title(pl.parameters{pl.corrsigs(i_fig,2),6},'FontSize', 8)
    grid on
    box on
    hline(0,'k')
    vline(TFA.baseline,'k')
    xlabel('time in ms')
    ylabel('amplitude modulation in %')
    
    %% calculate correlations
    corr.xdata = pl.data1(t.ind3:t.ind4,:);
    corr.xtimes = TFA.time(t.ind3:t.ind4);
    corr.ytimes = [flip(-median(diff(TFA.time)):-median(diff(TFA.time)):pl.corrtlim(1)) 0:median(diff(TFA.time)):pl.corrtlim(2)];
    corr.idx = (1:numel(corr.ytimes))-ceil((numel(corr.ytimes)/2));
    
    corr.r = nan(numel(corr.ytimes),numel(corr.xtimes));
    corr.p = corr.r;
    
    % loop across time-points
    fprintf('calculating %1.0f correlations, progress in %%: ', numel(corr.r))
    t.progress_n = round(linspace(1, size(corr.xdata,1), 11));
    t.progress_p = 0:10:100;
    for i_tp = 1:size(corr.xdata,1)
        if any(i_tp == t.progress_n)
            fprintf(' %3.0f', t.progress_p(i_tp == t.progress_n))
        end
        % index y-data
        t.idx1 = corr.idx+t.ind3+i_tp-1;
        t.idx2 = t.idx1>0;
        t.idx3 = find(t.idx2);
        
        % loop across lags: do separate correlations
        for i_corr = 1:numel(t.idx3)
            [t.R t.P] = corrcoef(...
                corr.xdata(i_tp,:),... % data of first signal at timepoint i_tp
                pl.data2(t.idx1(t.idx3(i_corr)),:)... % data of second signal for lag i_corr around timepoint i_tp
                );
            corr.r(t.idx3(i_corr),i_tp) = t.R(2);
            corr.p(t.idx3(i_corr),i_tp) = t.P(2);
        end
        
    end
    fprintf(' ...done!\n')
    
    %% plot r values
    subplot(5,2,[3 4 5 6])
    
    t.clims1=[-1 1]*max(max(abs(corr.r)));
    imagesc(corr.xtimes,corr.ytimes,corr.r,t.clims1)
    % colormap(gca,pl.colmap); % cbrewer('div','RdBu',256,'l')
    colormap(gca,flipud(cbrewer2('RdBu')))
    set(gca,'YDir','normal')
    title(sprintf('correlation R between %s and %s amplitudes of different lags',...
        pl.parameters{pl.corrsigs(i_fig,1),6},  pl.parameters{pl.corrsigs(i_fig,2),6}),'FontSize',8)
    h.cb1=colorbar;
    set(h.cb1,'FontSize',8)
    xlabel(sprintf('time in ms for %s', pl.parameters{pl.corrsigs(i_fig,1),6}))
    ylabel(sprintf('relative time lag in ms for %s', pl.parameters{pl.corrsigs(i_fig,2),6}))
    % freezeColors
    % h.cb1=cbfreeze(h.cb1);
    hline(0,'k'); vline(0,'k')
    set(gca,'FontSize',8)
  
    
    
    %% plot p-values


    subplot(5,2,[7 8 9 10])
    % ttest
    t.data = abs(log10(corr.p));
    t.clim = [0 max(t.data(:))];
    t.pcriterion = abs(log10(0.05));
    if max(t.data(:))<t.pcriterion
        % temp.colormap = repmat([0.5 0.5 0.5],100,1);
        t.colormap = repmat(linspace(1,0.3,1000)',1,3);
    else
        t.border = ceil((t.pcriterion / max(t.data(:)))*1000);
        % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
        t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
    end
    imagesc(corr.xtimes,corr.ytimes,t.data,t.clim)
    set(gca,'YDir','normal')
    title(sprintf('correlation p-values'...
        ),'FontSize',8)
    % freezeColors
    h.cb2=colorbar;
    xlabel(sprintf('time in ms for %s', pl.parameters{pl.corrsigs(i_fig,1),6}))
    ylabel(sprintf('relative time lag in ms for %s', pl.parameters{pl.corrsigs(i_fig,2),6}))
    t.yticks = get(h.cb2,'YTick');
    set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
        'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')),'FontSize',8)
    % h.cb2=cbfreeze(h.cb2);
    hline(0,'k'); vline(0,'k')
    set(gca,'FontSize',8)
    colormap(gca,t.colormap)
      
    
    
end


sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
% sav.filenames = {'Resp_AllSignals_Amp_Timecourse'};
% sav.filenames = {'Resp_AllSignals_Amp_Timecourse_v2'};
% sav.filenames = {'Resp_AllSignals_Amp_Timecourse_5'};
sav.filenames = {'Resp_AllSignals_Amp_Timecourse_sep_v2_5'};
% sav.filenames = {'Resp_AllSignals_Amp_Timecourse_v2_clustcorr_5b'};
for i_fig = 1:1
%     print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-dpng','-r300')
%     print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-djpeg','-r300')
%     saveas(figs{i_fig},fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'fig')
%     print(figs{i_fig},fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-depsc2', '-painters','-r300')
end


%% calculate cross-correlations between signals
clear figs h


pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
    [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [1 0 0] 2 'vis alpha' 'noRESS';...
    [10 14] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [0 0.5 0] 2 'mot alpha' 'noRESS';...
    [18 30] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [0 0.5 0.9] 2 'mot beta' 'noRESS';...
    [15 23] {'CP1';'CPz';'CP2';'P1';'Pz';'P2'} 'data_induced_bc' [0.8 0 0.8] 2 'parietal beta' 'noRESS';...
    };

pl.parameters = {...
    [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
    [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [255 151 0]./255 2 'vis alpha' 'noRESS';...
    [10 14] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [33 92 150]./255 2 'mot alpha' 'noRESS';...
    [18 30] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [170 99 57]./255 2 'mot beta' 'noRESS';...
    [15 23] {'CP1';'CPz';'CP2';'P1';'Pz';'P2'} 'data_induced_bc' [23 150 118]./255 2 'parietal beta' 'noRESS';...
    [2 4] {'C3';'CP3';'C1';'CP1'} 'data_induced_bc' [0.4 0.4 0.4] 2 'mot delta' 'noRESS';...
    };


% this one for the poster figures
% pl.parameters = {...
%     [14.16667 14.16667] {'RESS'} 'data_RESS_evoked_bc' [0 0 0] 2 sprintf('SSVEP') 'RESS';...
%     [8 12] {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'} 'data_induced_bc' [255 151 0]./255 2 'vis alpha' 'noRESS';...
%     [10 14] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [33 92 150]./255 2 'mot alpha' 'noRESS';...
%     [18 30] {'C3';'CP3';'C5';'CP5'} 'data_induced_bc' [170 99 57]./255 2 'mot beta' 'noRESS';...
%     [15 23] {'CP1';'CPz';'CP2';'P1';'Pz';'P2'} 'data_induced_bc' [23 150 118]./255 2 'parietal beta' 'noRESS';...
%     };

pl.corrsigs = [1 2; 1 3; 2 3];
pl.corrsigs = [3 2; 3 4];
pl.corrsigs = [4 5];
pl.corrsigs = [6 4; 6 3];
pl.corrtlim = [-1000 1000]; % time limits of shifted correlation analysis


pl.subs2use = 1:numel(F.Subjects2Use);

pl.xlim = [-3500 3500];

pl.permut_n = 10000;
% pl.permut_n = 100;


clear pl.tdata2 h.sp h.p h.pm h.pl
% do the analyisis
for i_fig = 1:size(pl.corrsigs,1)
    figs{i_fig}=figure;
    set(gcf,'Position',[100 100 500 950],'PaperPositionMode','auto')
    
    
    % first dataset
    % index frequencies
    [t.t t.ind1]=min(abs(TFA.frequency-pl.parameters{pl.corrsigs(i_fig,1),1}(1)));
    [t.t t.ind2]=min(abs(TFA.frequency-pl.parameters{pl.corrsigs(i_fig,1),1}(2)));
    
    % index time
    [t.t t.ind3]=min(abs(TFA.time-pl.xlim(1)));
    [t.t t.ind4]=min(abs(TFA.time-pl.xlim(2)));
    
    % index electrodes
    switch  pl.parameters{pl.corrsigs(i_fig,1),7}
        case 'noRESS'
            pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.parameters{pl.corrsigs(i_fig,1),2}, 'UniformOutput',false)),1));
            com=sprintf('pl.data1=squeeze(nanmean(nanmean(TFA.%s(t.ind1:t.ind2,:,pl.elec2plot_i,pl.subs2use),1),3));',pl.parameters{pl.corrsigs(i_fig,1),3});
        case 'RESS'
            pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes_RESS.labels},x), pl.parameters{pl.corrsigs(i_fig,1),2}, 'UniformOutput',false)),1));
            com=sprintf('pl.data1=squeeze(nanmean(TFA.%s(t.ind1:t.ind2,:,pl.subs2use),1));',pl.parameters{pl.corrsigs(i_fig,1),3});
        
    end
    eval(com)
    [tt.h1 tt.p1 tt.ci1 tt.stats1]=ttest(pl.data1(t.ind3:t.ind4,:)');
    
    % second dataset
    % index frequencies
    [t.t t.ind1]=min(abs(TFA.frequency-pl.parameters{pl.corrsigs(i_fig,2),1}(1)));
    [t.t t.ind2]=min(abs(TFA.frequency-pl.parameters{pl.corrsigs(i_fig,2),1}(2)));
    
    % index time
    [t.t t.ind3]=min(abs(TFA.time-pl.xlim(1)));
    [t.t t.ind4]=min(abs(TFA.time-pl.xlim(2)));
    
    % index electrodes
    switch  pl.parameters{pl.corrsigs(i_fig,2),7}
        case 'noRESS'
            pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.parameters{pl.corrsigs(i_fig,2),2}, 'UniformOutput',false)),1));
            com=sprintf('pl.data2=squeeze(nanmean(nanmean(TFA.%s(t.ind1:t.ind2,:,pl.elec2plot_i,pl.subs2use),1),3));',pl.parameters{pl.corrsigs(i_fig,2),3});
        case 'RESS'
            pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes_RESS.labels},x), pl.parameters{pl.corrsigs(i_fig,2),2}, 'UniformOutput',false)),1));
            com=sprintf('pl.data2=squeeze(nanmean(TFA.%s(t.ind1:t.ind2,:,pl.subs2use),1));',pl.parameters{pl.corrsigs(i_fig,2),3});
        
    end
    eval(com)
    [tt.h2 tt.p2 tt.ci2 tt.stats2]=ttest(pl.data2(t.ind3:t.ind4,:)');
    
    
    %% plot data
    
    % data1
    % actual plot with boundedline
    pl.mdata = mean(pl.data1(t.ind3:t.ind4,:),2)';
    pl.ddata = tt.ci1 - repmat(mean(pl.data1(t.ind3:t.ind4,:),2)',2,1);
    %pl.ddata = tt.ci;
    subplot(5,2,1)
    [h.l(1), h.p(1)] = boundedline(TFA.time(t.ind3:t.ind4),pl.mdata,pl.ddata(1,:));
    set(h.l(1),'Color',pl.parameters{pl.corrsigs(i_fig,1),4},'LineWidth',pl.parameters{pl.corrsigs(i_fig,1),5});
    set(h.p(1),'FaceColor',pl.parameters{pl.corrsigs(i_fig,1),4},'FaceAlpha',0.25)
    hold on
    
    % plot single subject data
%     pl.mdata = mean(pl.data,2)';
%     h.p(i_pl,:)=plot(TFA.time(t.ind3:t.ind4),pl.data,'Color',[pl.parameters{i_pl,4} 0.3],'LineWidth',0.2);
%     hold on
%     h.pm(i_pl)=plot(TFA.time(t.ind3:t.ind4),pl.mdata,'Color',[pl.parameters{i_pl,4} 1],'LineWidth',2);
    
    % significant effects?
    pl.tdata2=nan(size(pl.mdata));
    pl.tdata2(tt.p1<.05)=1;
    set(gca,'ylim',[-1.1 (1+0.15+0.05*4)]...
        .*max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false))))
    pl.maxdata = max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false)));
    pl.y = 1*pl.maxdata*(1+0.05+0.05*2);
    plot(TFA.time(t.ind3:t.ind4),pl.tdata2.*pl.y,'Color',[pl.parameters{pl.corrsigs(i_fig,1),4} 1],'LineWidth',6)
    
    
    % graphics
    xlim(pl.xlim)
    title(pl.parameters{pl.corrsigs(i_fig,1),6},'FontSize', 8)
    grid on
    box on
    hline(0,'k')
    vline(TFA.baseline,'k')
    xlabel('time in ms')
    ylabel('amplitude modulation in %')
    
    
    % data2
    % actual plot with boundedline
    pl.mdata = mean(pl.data2(t.ind3:t.ind4,:),2)';
    pl.ddata = tt.ci2 - repmat(mean(pl.data2(t.ind3:t.ind4,:),2)',2,1);
    %pl.ddata = tt.ci;
    subplot(5,2,2)
    [h.l(1), h.p(1)] = boundedline(TFA.time(t.ind3:t.ind4),pl.mdata,pl.ddata(1,:));
    set(h.l(1),'Color',pl.parameters{pl.corrsigs(i_fig,2),4},'LineWidth',pl.parameters{pl.corrsigs(i_fig,2),5});
    set(h.p(1),'FaceColor',pl.parameters{pl.corrsigs(i_fig,2),4},'FaceAlpha',0.25)
    hold on
    
    % plot single subject data
%     pl.mdata = mean(pl.data,2)';
%     h.p(i_pl,:)=plot(TFA.time(t.ind3:t.ind4),pl.data,'Color',[pl.parameters{i_pl,4} 0.3],'LineWidth',0.2);
%     hold on
%     h.pm(i_pl)=plot(TFA.time(t.ind3:t.ind4),pl.mdata,'Color',[pl.parameters{i_pl,4} 1],'LineWidth',2);
    
    % significant effects?
    pl.tdata2=nan(size(pl.mdata));
    pl.tdata2(tt.p2<.05)=1;
    set(gca,'ylim',[-1.1 (1+0.15+0.05*4)]...
        .*max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false))))
    pl.maxdata = max(cell2mat(cellfun(@(x) max(abs(double(x))), get(get(gca,'Children'),'YData'),'UniformOutput',false)));
    pl.y = 1*pl.maxdata*(1+0.05+0.05*2);
    plot(TFA.time(t.ind3:t.ind4),pl.tdata2.*pl.y,'Color',[pl.parameters{pl.corrsigs(i_fig,2),4} 1],'LineWidth',6)
    
    
    % graphics
    xlim(pl.xlim)
    title(pl.parameters{pl.corrsigs(i_fig,2),6},'FontSize', 8)
    grid on
    box on
    hline(0,'k')
    vline(TFA.baseline,'k')
    xlabel('time in ms')
    ylabel('amplitude modulation in %')
    
    %% calculate correlations
    corr.xdata = pl.data1(t.ind3:t.ind4,:);
    corr.xtimes = TFA.time(t.ind3:t.ind4);
    corr.ydata = pl.data2(t.ind3:t.ind4,:);
    corr.ytimes = TFA.time(t.ind3:t.ind4);
    
    
    corr.r = nan(numel(corr.ytimes),numel(corr.xtimes));
    corr.p = corr.r;
    
    % loop across time-points
    fprintf('calculating %1.0f correlations, progress in %%: ', numel(corr.r))
    t.progress_n = round(linspace(1, size(corr.xdata,1), 11));
    t.progress_p = 0:10:100;
    for i_tp = 1:size(corr.xdata,1)
        if any(i_tp == t.progress_n)
            fprintf(' %3.0f', t.progress_p(i_tp == t.progress_n))
        end
        
        % loop time of second signal
        for i_corr = 1:numel(corr.ytimes)
            [t.R t.P] = corrcoef(...
                corr.xdata(i_tp,:),... % data of first signal at timepoint i_tp
                corr.ydata(i_corr,:)... % data of second signal for lag i_corr around timepoint i_tp
                );
            corr.r(i_corr,i_tp) = t.R(2);
            corr.p(i_corr,i_tp) = t.P(2);
        end
        
    end
    fprintf(' ...done!\n')
    
    %% plot r values
    subplot(5,2,[3 4 5 6])
    
    t.clims1=[-1 1]*max(max(abs(corr.r)));
    imagesc(corr.xtimes,corr.ytimes,corr.r,t.clims1)
    % colormap(gca,pl.colmap); % cbrewer('div','RdBu',256,'l')
    colormap(gca,flipud(cbrewer2('RdBu')))
    set(gca,'YDir','normal')
    title(sprintf('cross time correlation R between %s and %s amplitudes',...
        pl.parameters{pl.corrsigs(i_fig,1),6},  pl.parameters{pl.corrsigs(i_fig,2),6}),'FontSize',8)
    h.cb1=colorbar;
    set(h.cb1,'FontSize',8)
    xlabel(sprintf('time in ms for %s', pl.parameters{pl.corrsigs(i_fig,1),6}))
    ylabel(sprintf('time in ms for %s', pl.parameters{pl.corrsigs(i_fig,2),6}))
    % freezeColors
    % h.cb1=cbfreeze(h.cb1);
    hline(0,'k'); vline(0,'k')
    set(gca,'FontSize',8)
  
    
    
    %% plot p-values


    subplot(5,2,[7 8 9 10])
    % ttest
    t.data = abs(log10(corr.p));
    t.data(t.data>abs(log10(0.00001))) = abs(log10(0.0000099));
    t.clim = [0 max(t.data(:))];
    t.pcriterion = abs(log10(0.05));
    if max(t.data(:))<t.pcriterion
        % temp.colormap = repmat([0.5 0.5 0.5],100,1);
        t.colormap = repmat(linspace(1,0.3,1000)',1,3);
    else
        t.border = ceil((t.pcriterion / max(t.data(:)))*1000);
        % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
        t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
    end
    imagesc(corr.xtimes,corr.ytimes,t.data,t.clim)
    set(gca,'YDir','normal')
    title(sprintf('correlation p-values'...
        ),'FontSize',8)
    % freezeColors
    h.cb2=colorbar;
    xlabel(sprintf('time in ms for %s', pl.parameters{pl.corrsigs(i_fig,1),6}))
    ylabel(sprintf('time in ms for %s', pl.parameters{pl.corrsigs(i_fig,2),6}))
    t.yticks = get(h.cb2,'YTick');
    set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
        'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')),'FontSize',8)
    % h.cb2=cbfreeze(h.cb2);
    hline(0,'k'); vline(0,'k')
    set(gca,'FontSize',8)
    colormap(gca,t.colormap)
      
    
    
end


sav.pathout = 'C:\Users\psy05cvd\Dropbox\work\matlab\AnalyzerUni\SSVEP_volmov\figures\';
% sav.filenames = {'Resp_AllSignals_Amp_Timecourse'};
% sav.filenames = {'Resp_AllSignals_Amp_Timecourse_v2'};
% sav.filenames = {'Resp_AllSignals_Amp_Timecourse_5'};
sav.filenames = {'Resp_AllSignals_Amp_Timecourse_sep_v2_5'};
% sav.filenames = {'Resp_AllSignals_Amp_Timecourse_v2_clustcorr_5b'};
for i_fig = 1:1
%     print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-dpng','-r300')
%     print(figs{i_fig}, fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-djpeg','-r300')
%     saveas(figs{i_fig},fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'fig')
%     print(figs{i_fig},fullfile(sav.pathout,sprintf('%s',sav.filenames{i_fig})),'-depsc2', '-painters','-r300')
end


