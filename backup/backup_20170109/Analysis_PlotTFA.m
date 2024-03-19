%% plot previously calculated TFA data


clear all
%% parameters
% F.PathIn                = 'F:\work\data\SSVEP_volmov\EEG\TFA';
F.PathIn                = 'F:\work\data\SSVEP_volmov\EEG\TFA_SCADS';
% F.PathIn                = 'F:\work\data\SSVEP_volmov\EEG\TFA_SCADS2';
F.Subjects2Use          = [1:20];

TFA.baseline            = [-3000 -2500];
% TFA.baseline            = [-2250 -2000];

pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];
%% loop for subjects
for i_sub=1:numel(F.Subjects2Use)
    
    % read in TFA data
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%02.0f_exp_tfa.mat ||\n',i_sub,numel(F.Subjects2Use),F.PathIn,F.Subjects2Use(i_sub))
    
    temp.tfa = open(sprintf('%s\\VP%02.0f_exp_tfa.mat',F.PathIn,F.Subjects2Use(i_sub)));
    
    % preallocate memory
    if i_sub == 1
        TFA.data_induced = nan([size(temp.tfa.TFA.data_induced) numel(F.Subjects2Use)]);
        TFA.data_evoked = TFA.data_induced;
        TFA.data_induced_bc = TFA.data_induced;
        TFA.data_evoked_bc = TFA.data_induced;
        TFA.time = temp.tfa.TFA.t;
        TFA.frequency = temp.tfa.TFA.f;
        TFA.electrodes = temp.tfa.TFA.electrodes;
    end
    
    % subject by subject
    TFA.data_induced(:,:,:,i_sub) = temp.tfa.TFA.data_induced; % induced data
    TFA.data_induced_bc(:,:,:,i_sub) = bsxfun(@minus, temp.tfa.TFA.data_induced, ...
        mean(temp.tfa.TFA.data_induced(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2));
    TFA.data_evoked(:,:,:,i_sub) = temp.tfa.TFA.data_evoked; % evoked data
    TFA.data_evoked_bc(:,:,:,i_sub) = bsxfun(@minus, temp.tfa.TFA.data_evoked, ...
        mean(temp.tfa.TFA.data_evoked(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:),2));
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
    
end

%% actual plotting
% figure; topoplot([],TFA.electrodes, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);


% plotting parameters
% pl.elec2plot = {'Oz'};
% pl.elec2plot = {'PO3';'POz';'PO4';'O1';'Oz';'O2';'Iz'}; % steady state I
% pl.elec2plot = {'POz';'O1';'Oz';'O2';'I1';'Iz';'I2'}; % steady state II
% pl.elec2plot = {'PO4';'O2';'PO8';'P8';'P10';'I2'}; % vis alpha
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; % vis alpha II
% pl.elec2plot = {'POz'}; % vis alpha II
% pl.elec2plot = {'POz'}; % vis alpha IV
% pl.elec2plot = {'C3';'CP3'}; % motor alpha/beta
% pl.elec2plot = {'C3';'CP3';'C5';'CP5'}; % motor alpha/beta II
pl.elec2plot = {TFA.electrodes(1:64).labels}';
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.xlims=[-2500 2500]; % index time 2 plot
% pl.xlims=[-3500 3500]; % index time 2 plot
[t.t t.ind1]=min(abs(TFA.time-pl.xlims(1)));
[t.t t.ind2]=min(abs(TFA.time-pl.xlims(2)));

pl.flims = TFA.frequency([1 end]); % index frequency 2 plot
% pl.flims = [5 40];
pl.flims = [8 40];
[t.t t.find1]=min(abs(TFA.frequency-pl.flims(1)));
[t.t t.find2]=min(abs(TFA.frequency-pl.flims(2)));

pl.subs2use = F.Subjects2Use;



%%%%%%%%%%%%%
% raw induced
pl.data_ind = squeeze(nanmean(nanmean(TFA.data_induced(:,:,pl.elec2plot_i,pl.subs2use),3),4));
pl.data_evo = squeeze(nanmean(nanmean(TFA.data_evoked(:,:,pl.elec2plot_i,pl.subs2use),3),4));
t.clims1=[0 1]*max(max(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2)));
t.clims2=[0 1]*max(max(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2)));
figure;
subplot(2,1,1)
imagesc(TFA.time,TFA.frequency,pl.data_ind,t.clims1)
set(gca,'YDir','normal')
title(sprintf('tfa, for channel %s; averaged in frequency domain', vararg2str(pl.elec2plot)))
colorbar
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
hline(14.16667,'m')

subplot(2,1,2)
imagesc(TFA.time,TFA.frequency,pl.data_evo,t.clims2)
set(gca,'YDir','normal')
title(sprintf('tfa, for channel %s; averaged in time domain', vararg2str(pl.elec2plot)))
colorbar
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
hline(14.16667,'m')

% normalized induced
pl.data_ind = squeeze(nanmean(nanmean(TFA.data_induced_bc(:,:,pl.elec2plot_i,pl.subs2use),3),4));
t.clims1=[-1 1]*max(max(abs(pl.data_ind(t.find1:t.find2,t.ind1:t.ind2))));
figure;
h.sp1=subplot(2,1,1);
colormap('default')
% colormap(cmap('R2'))
imagesc(TFA.time,TFA.frequency,pl.data_ind,t.clims1)
set(gca,'YDir','normal')
title(sprintf('normalized tfa, for channel %s; averaged in frequency domain; baseline [%1.0f %1.0f]ms', vararg2str(pl.elec2plot),TFA.baseline))
h.cb1=colorbar;
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
freezeColors
h.cb1=cbfreeze(h.cb1);
hline(14.16667,'m')

h.sp2=subplot(2,1,2);
for i_freq = 1:numel(TFA.frequency)
    [tt.h tt.p tt.ci tt.stats]=ttest(squeeze(nanmean(TFA.data_induced_bc(i_freq,:,pl.elec2plot_i,:),3))');
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
colormap(t.colormap)
imagesc(TFA.time,TFA.frequency,abs(log10(TFA.pvals_induced_bc)),t.clim)
set(gca,'YDir','normal')
title(sprintf('p-values, for channel %s; averaged in frequency domain; baseline [%1.0f %1.0f]ms', vararg2str(pl.elec2plot),TFA.baseline))
freezeColors
h.cb2=colorbar;
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
t.yticks = get(h.cb2,'YTick');
set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
    'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')))
h.cb2=cbfreeze(h.cb2);
hline(14.16667,'g')



% normalized evoked
pl.data_evo = squeeze(nanmean(nanmean(TFA.data_evoked_bc(:,:,pl.elec2plot_i,pl.subs2use),3),4));
t.clims2=[-1 1]*max(max(abs(pl.data_evo(t.find1:t.find2,t.ind1:t.ind2))));
figure
h.sp1=subplot(2,1,1);
colormap('default')
imagesc(TFA.time,TFA.frequency,pl.data_evo,t.clims2)
set(gca,'YDir','normal')
title(sprintf('normalized tfa, for channel %s; averaged in time domain; baseline [%1.0f %1.0f]ms', vararg2str(pl.elec2plot),TFA.baseline))
freezeColors

h.cb1=colorbar;
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)

h.cb1=cbfreeze(h.cb1);
hline(14.16667,'m')

h.sp2=subplot(2,1,2);
for i_freq = 1:numel(TFA.frequency)
    [tt.h tt.p tt.ci tt.stats]=ttest(squeeze(nanmean(TFA.data_evoked_bc(i_freq,:,pl.elec2plot_i,:),3))');
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
colormap(t.colormap)
imagesc(TFA.time,TFA.frequency,abs(log10(TFA.pvals_induced_bc)),t.clim)
set(gca,'YDir','normal')
title(sprintf('p-values, for channel %s; averaged in time domain; baseline [%1.0f %1.0f]ms', vararg2str(pl.elec2plot),TFA.baseline))
freezeColors
h.cb2=colorbar;
xlabel('time in ms')
ylabel('frequency in Hz')
xlim(pl.xlims)
ylim(pl.flims)
t.yticks = get(h.cb2,'YTick');
set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
    'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')))
h.cb2=cbfreeze(h.cb2);
hline(14.16667,'g')


%% plot timecourse of frequency
% pl.freq2plot=[14.16667 14.16667];
% pl.freq2plot=[10 14];
% pl.freq2plot=[15 25];
% pl.freq2plot=[20 30];
pl.freq2plot=[8 14];
% pl.freq2plot=[15 30];
[t.t t.ind3]=min(abs(TFA.frequency-pl.freq2plot(1)));
[t.t t.ind4]=min(abs(TFA.frequency-pl.freq2plot(end)));

% pl.elec2plot = {'Oz';'POz'};
% pl.elec2plot = {'PO3';'POz';'PO4';'O1';'Oz';'O2';'Iz'}; % steady state I
% pl.elec2plot = {'POz';'O1';'Oz';'O2';'I1';'Iz';'I2'}; % steady state II
% pl.elec2plot = {'PO4';'O2';'PO8';'P8';'P10';'I2'}; % vis alpha I
% pl.elec2plot = {'POz'}; % vis alpha II
pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; % vis alpha III
% pl.elec2plot = {'POz';'Pz';'P1';'P2'}; % vis alpha II
% pl.elec2plot = {'C3';'CP3'}; % motor alpha/beta
% pl.elec2plot = {'C3';'CP3';'C5';'CP5'}; % motor alpha/beta II
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.xlims=[-2750 3000]; % index time 2 plot
[t.t t.ind1]=min(abs(TFA.time-pl.xlims(1)));
[t.t t.ind2]=min(abs(TFA.time-pl.xlims(2)));

pl.subs2use = F.Subjects2Use;

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

pl.parameters = {...
    [14.16667 14.16667] {'POz';'O1';'Oz';'O2';'I1';'Iz';'I2'} 'data_evoked_bc' [0.2 0.2 0.2] 2 sprintf('vis evo\nSSVEP');...
    [8 12] {'PO4';'O2';'PO8';'P8';'P10';'I2'} 'data_induced_bc' [1 0 0] 1 'vis alpha';...
    [18 30] {'PO4';'O2';'PO8';'P8';'P10';'I2'} 'data_induced_bc' [1 0 1] 1 'vis beta';...
    [10 14] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0] 1 'mot alpha';...
    [18 30] {'C3';'CP3'} 'data_induced_bc' [0 0.5 0.9] 1 'mot beta';...
    };
pl.subs2use = F.Subjects2Use;
pl.xlim = [-3000 3000];

pl.toponum = [ceil(sqrt(size(pl.parameters,1))) round(sqrt(size(pl.parameters,1)))];
pl.plpos = repmat([1:pl.toponum(2)*3],pl.toponum(1),1)+repmat((0:pl.toponum(2)*4:pl.toponum(2)*4*(pl.toponum(1)-1))',1,pl.toponum(2)*3);
pl.topopos = repmat([pl.toponum(2)*3+1:pl.toponum(2)*4],pl.toponum(1),1)'+repmat((0:pl.toponum(2)*4:pl.toponum(2)*4*(pl.toponum(1)-1))',1,pl.toponum(2))';

figure;
subplot(pl.toponum(1),pl.toponum(2)*4,pl.plpos(:))
for i_pl = 1:size(pl.parameters)
    % index frequencies
    [t.t t.ind1]=min(abs(TFA.frequency-pl.parameters{i_pl,1}(1)));
    [t.t t.ind2]=min(abs(TFA.frequency-pl.parameters{i_pl,1}(2)));
    
    % index time
    [t.t t.ind3]=min(abs(TFA.time-pl.xlim(1)));
    [t.t t.ind4]=min(abs(TFA.time-pl.xlim(2)));
    
    % index electrodes
    pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.parameters{i_pl,2}, 'UniformOutput',false)),1));
    com=sprintf('pl.data=squeeze(nanmean(nanmean(nanmean(TFA.%s(t.ind1:t.ind2,t.ind3:t.ind4,pl.elec2plot_i,pl.subs2use),1),3),4));',pl.parameters{i_pl,3});
    eval(com)
    
    % actual plot
    plot(TFA.time(t.ind3:t.ind4),pl.data,'Color',pl.parameters{i_pl,4},'LineWidth',pl.parameters{i_pl,5})
    hold on
end
xlim(pl.xlim)
set(gca,'ylim',[-1.1 1.1].*max(cellfun(@(x) max(abs(x)), get(get(gca,'Children'),'YData'))))
legend(pl.parameters(:,6),'location','NorthWest','FontSize', 8)
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

%% topoplot
pl.time2plot=[500 2000];
% pl.time2plot=[-1500 500];
% pl.time2plot=[-500 300];
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

pl.titleadd = {'induced';'evoked'};

figure;
% loop for evoked/induced
for i_dat = 1:2
    colormap('default')
    h.s1=subplot(2,3,1+3*(i_dat-1));
    pl.data=squeeze(nanmean(pl.data_raw(:,:,i_dat),2));
    topoplot( pl.data, TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',[0 max(pl.data)]);
    title(sprintf('raw %s amp %1.2f to %1.2f Hz;\n%1.0f to %1.0f ms',pl.titleadd{i_dat},pl.freq2plot, pl.time2plot))
    freezeColors
    h.c1 = colorbar;
    h.c1=cbfreeze(h.c1);
    
    colormap('default')
    h.s2=subplot(2,3,2+3*(i_dat-1));
    pl.data=squeeze(nanmean(pl.data_bc(:,:,i_dat),2));
    topoplot( pl.data, TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits','absmax');
    title(sprintf('normalized %s amp %1.2f to %1.2f Hz;\n%1.0f to %1.0f ms',pl.titleadd{i_dat},pl.freq2plot, pl.time2plot))
    freezeColors
    h.c2 = colorbar;
    h.c2=cbfreeze(h.c2);
    
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
        'shading', 'interp', 'numcontour', 0, 'maplimits',t.clim);
    colormap(t.colormap)
    title(sprintf('p-values %s %1.2f to %1.2f Hz;\n%1.0f to %1.0f ms',pl.titleadd{i_dat},pl.freq2plot, pl.time2plot))
    freezeColors
    h.cb2 = colorbar;
    t.yticks = get(h.cb2,'YTick');
    set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
        'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')))
    h.cb2=cbfreeze(h.cb2);
    
end

%% plot topographies across time


% define parameters
t.p_time = [100 100 -2200 2500]; % width step min max
t.posScale = 1.1;

pl.freq2plot=[14.16667 14.16667];
% pl.freq2plot=[10 14];
% pl.freq2plot=[16 25];
% pl.freq2plot=[8 12];
% pl.freq2plot=[15 30];
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
t.col = round(sqrt((size(t.time,1)+2)/(1/(t.t(1)/t.t(2)))));
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
        colormap('default')
        topoplot( plotdata(i_fig,:,i_spl), TFA.electrodes(1:64), ...
            'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims,'electrodes','off');
        title(sprintf('[%1.0f %1.0f]',t.time(i_spl,1),t.time(i_spl,2)),'FontSize',8)
        t.pos = get(h.sp(i_spl),'Position');
        set(h.sp(i_spl),'Position',[t.pos(1:2)-(t.pos(3:4).*((t.posScale-1)/2)) t.pos(3:4).*t.posScale])
    end
    h.sp(i_spl)=subplot(t.row,t.col,i_spl+2);
    topoplot( plotdata(i_fig,:,i_spl), TFA.electrodes(1:64),  ...
        'style','blank');
    title(sprintf('%s\n[%1.2f %1.2f]Hz',t.conlabel{i_fig},pl.freq2plot),'FontSize',8)
    
    t.pos2 = get(h.sp(i_spl),'Position');
    t.pos3 = get(h.sp(i_spl),'OuterPosition');
    h.a1 = axes('position',[t.pos3(1) t.pos2(2) t.pos3(3) t.pos2(4)],'Visible','off');
    colormap('default')
    caxis(t.lims);
    h.c3 = colorbar('clim',t.lims);
    t.pos4 = get(h.c3,'Position');
    set(h.c3,'Position',[t.pos4(1)+0.03 t.pos4(2)+(t.pos4(4)*(1/6)) t.pos4(3)/2 t.pos4(4)*(2/3)])
end
axcopy(h.fig(i_fig))
