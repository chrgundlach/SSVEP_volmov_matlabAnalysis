%% script to process files
% matlabpool open 3



clear all
%% parameters
% F.Rootpath          = '/scr/curie2/christopher/PinchLearn_EEG/';
F.Rootpath          = 'F:\work\data\SSVEP_volmov\EEG\';
F.PathRAW           = [F.Rootpath 'RAW'];
F.PathSET           = [F.Rootpath 'SET'];
F.PathTFA           = [F.Rootpath 'TFA'];
F.PathICA           = [F.Rootpath 'ICA2'];
F.PathBehav         = [F.Rootpath 'PRESENTATION'];
% F.ChannFile         = ['/home/raid1/gundlach/Eigene Dateien/additional software/eeglab/eeglab9_0_5_6b/plugins/dipfit2.2/standard_BESA/standard-10-5-cap385.elp'];
F.ChannFile         = 'C:\Users\HP-User\matlab\Skripte\Auswertungsskripte\Analyzer_G\ChanLocs\BioSemi64_8_1020.epf';
F.Subjects          = ['01';];
F.Sub2Use           = [1];
F.Prefix            = ['VP'];
F.PresentName       = ['REFerence_lvl.txt'];
F.EEGChans          = [1:64];
F.Trigger           = {{'31'};{'26', '25'}};
F.Epoch             = {[-2.5 2.5];[0 2]};

F.Resample          = 250;
t={};

%% script
for i_sub = 1:numel(F.Sub2Use)
    %% read in raw data
    % index files
    h.time1 = tic;
    t{1}.dir = dir([F.PathRAW '\' F.Prefix F.Subjects(F.Sub2Use(i_sub),:) '_m*.bdf']);
    t{2}.dir = dir([F.PathRAW '\' F.Prefix F.Subjects(F.Sub2Use(i_sub),:) '_a*.bdf']);
    
    %%
    for i_file = 1:2
        % loop for potential files
        %     fprintf(1,'###\n%s\t subject %02.0f of %02.0f\t Reading in file: %s\n###\n',datestr(now),i_sub,numel(F.Sub2Use),t{i_sub}.dir.name)
        % load
        EEG=pop_biosig(sprintf('%s\\%s', F.PathRAW,t{i_file}.dir.name));

        if ~isempty(F.Resample)
            EEG=pop_resample(EEG,F.Resample);
        end
        
        EEG.setname=t{i_file}.dir.name(1:end-4);
        % chanlocs
        EEG.chanlocs=pop_chanedit(EEG.chanlocs,'load',{F.ChannFile,'filetype','besa (elp)'}); % Load Channel Locations
        % select data channels
        EEG = pop_select( EEG,'channel',F.EEGChans);
        % reref
        EEG = pop_reref( EEG, []);
        
        % epoch
        EEG = pop_epoch( EEG, F.Trigger{i_file}, F.Epoch{i_file}, 'epochinfo', 'yes');
        
        EEG = eeg_checkset( EEG );
        EEG = pop_rmbase( EEG, []);
        
        % discard trials with strong activity
        
        clean.param=8;
        for i_el=1:EEG.nbchan
            tt.data=squeeze(EEG.data(i_el,:,:));
            clean.m(i_el,:)=mean(tt.data(:));
            clean.std(i_el,:)=std(tt.data(:));
        end
        for i_tr = 1:EEG.trials
            tt.data=EEG.data(:,:,i_tr);
            clean.trials(i_tr)=~(any(any(tt.data<repmat(clean.m-clean.std*clean.param,1,EEG.pnts)))|...
                any(any(tt.data>repmat(clean.m+clean.std*clean.param,1,EEG.pnts))));
        end
        EEG = pop_select( EEG,'trial',find(clean.trials));
        
        
        %% run ica
        h.time2 = tic;
        %         [ EEG.icaweights, EEG.icasphere, mods ] = runamica12(EEG.data(:,:) );
        EEG = pop_runica(EEG, 'extended',1,'interupt','on');
        time(i_sub,1) = toc(h.time2);
        
        % select one epoch for saving
        EEG = pop_select( EEG,'trial',1);
        
        % save
        EEG.setname=t{i_file}.dir.name(1:end-4);
        EEG = pop_saveset( EEG, 'filename',[EEG.setname '.set'],'filepath',F.PathICA);
        
        EEG = [];
        time(i_sub,2) = toc(h.time1);
    end

    
end