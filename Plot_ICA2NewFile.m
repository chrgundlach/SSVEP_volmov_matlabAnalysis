%% script to check ICA decompositions
% loads set file with ICA-weights
% loads raw file, epochs, resamples and rereferences data
% get weights from first file


clear all
%% parameters
F.Rootpath          = 'F:\work\data\SSVEP_volmov\EEG\';
F.PathRAW           = [F.Rootpath 'RAW'];
F.PathSET           = [F.Rootpath 'SET'];
F.PathTFA           = [F.Rootpath 'TFA'];
F.PathICA           = [F.Rootpath 'ICA'];
F.PathBehav         = [F.Rootpath 'PRESENTATION'];
% F.ChannFile         = ['/home/raid1/gundlach/Eigene Dateien/additional software/eeglab/eeglab9_0_5_6b/plugins/dipfit2.2/standard_BESA/standard-10-5-cap385.elp'];
F.ChannFile         = 'C:\Users\HP-User\matlab\Skripte\Auswertungsskripte\Analyzer_G\ChanLocs\BioSemi64_8_1020.epf';
F.Subjects          = ['01';];
F.Sub2Use           = [1];
F.Prefix            = ['VP'];
F.Suffix            = {'_motor';'_alpha'};
F.PresentName       = ['REFerence_lvl.txt'];
F.EEGChans          = [1:64];
F.Trigger           = {{'31'};{'26', '25'}};
F.Epoch             = {[-2.5 2.5];[0 2]};

F.Resample          = 250;

pl.file             = 2; % 1=motor; 2=visual

%% load ICA file
ALLEEG = []; EEG = [];
EEG = pop_loadset('filename',[F.Prefix F.Subjects(F.Sub2Use,:) F.Suffix{pl.file} '.set'],'filepath',F.PathICA );
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

%% load RAW file
t.dir{1} = dir([F.PathRAW '\' F.Prefix F.Subjects(F.Sub2Use,:) '_m*.bdf']);
t.dir{2} = dir([F.PathRAW '\' F.Prefix F.Subjects(F.Sub2Use,:) '_a*.bdf']);
EEG=pop_biosig(sprintf('%s\\%s', F.PathRAW,t.dir{pl.file}.name));
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

% resample
if ~isempty(F.Resample)
    EEG=pop_resample(EEG,F.Resample);
end

EEG.setname=t.dir{pl.file}.name(1:end-4);
        

% chanlocs
EEG.chanlocs=pop_chanedit(EEG.chanlocs,'load',{F.ChannFile,'filetype','besa (elp)'}); % Load Channel Locations
% select data channels
EEG = pop_select( EEG,'channel',F.EEGChans);
% reref
EEG = pop_reref( EEG, []);

% epoch
EEG = pop_epoch( EEG, F.Trigger{pl.file}, F.Epoch{pl.file}, 'epochinfo', 'yes');
EEG = eeg_checkset( EEG );
EEG = pop_rmbase( EEG, []);

% get ICA information
EEG = pop_editset(EEG, 'icachansind', 'ALLEEG(1).icachansind', 'icaweights', 'ALLEEG(1).icaweights', 'icasphere', 'ALLEEG(1).icasphere');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% show ICA components
pop_selectcomps(EEG, [1:64] );
