function SEEG_preproc(cfg)
%% SEEG_data_preprocessing
% The current script includes: 1) Getting triggers, 2) Resample if necessary, 
%            3)Re-reference to average, 4)Epoch data

%% Load parameters
clc,clearvars('-except','cfg')
v2struct(cfg);

%% Load data
% Getting file path for storing data, by choose the .edf file manually.
[~, filepath]=uigetfile('*.edf'); 
% Getting data files
files = dir(fullfile(filepath,'*.edf'));
% Folder for storing preprocessed data
if ~exist(writdir,'dir')
    mkdir(writdir)
end
% Define the number of participants for preprocessing
if strcmp(single_subj,'yes')
    sub_num = 1;
    dataset= [filepath single_name];
else
    sub_num = length(files); 
end

%% Start preprocessing
for id = 1:sub_num
    fprintf('Running subject %i of %i for preprocessing...\n',id,sub_num)
    %clear previous dataset
    clear EEG
    %load dataset when preprocessing multiple datasets
    if strcmp(single_subj,'no')
        dataset= [filepath files(id).name];
        f_name = files(id).name(1:end-4);
    else
        f_name = single_name(1:end-4);
    end
    sub  = str2num([f_name(1:2)]);
    EEG = pop_biosig(dataset);
    
%% Get triggers
    % Get channel lables
    chan_labels = transpose({EEG.chanlocs(:).labels});
    % Load trigger function
    EEG = creat_events1(EEG,dc_chans,chan_labels,cfg); 

%% Resampling
    EEG = pop_resample(EEG,res_rate);
%% Re-reference
    EEG = pop_select(EEG,'nochannel',irrelevant_channels);
    EEG = pop_reref(EEG,[],'refstate','averef');
%% Filtering, and filter 50, 100, and 150 Hz
    fltord = 2; 
    if  strcmp(filt_type,'fieldtrip')
        cfg.bsfilter      = 'yes' ;%'no' or 'yes' bandstop filter (default = 'no')
        cfg.bsfreq        = [49 51];
        cfg.bsfiltord     =  2;
        cfg.bsfilttype    =  'but';
        data= eeglab2fieldtrip(EEG,'raw','none');
        data = ft_preprocessing(cfg,data);
        cfg.bsfreq        = [99 101];
        data = ft_preprocessing(cfg,data);
        cfg.bsfreq        = [149 151];
        data = ft_preprocessing(cfg,data);
        
        EEG.data = data.trial{1};
        EEG = pop_basicfilter(EEG,1:EEG.nbchan,'Filter','bandpass','Design','butter','Cutoff',bandcut,'Order',fltord); 
        
    elseif strcmp(filt_type,'eeglab')
        EEG  = pop_basicfilter( EEG,  1:EEG.nbchan, 'Cutoff',  50, 'Design', 'notch', 'Filter', 'PMnotch', 'Order',  180 );
        EEG  = pop_basicfilter( EEG,  1:EEG.nbchan, 'Cutoff',  100, 'Design', 'notch', 'Filter', 'PMnotch', 'Order',  180 );
        EEG  = pop_basicfilter( EEG,  1:EEG.nbchan, 'Cutoff',  150, 'Design', 'notch', 'Filter', 'PMnotch', 'Order',  180 );
        EEG  = pop_basicfilter( EEG,  1:EEG.nbchan, 'Filter','highpass','Design','butter','Cutoff',bandcut,'Order', fltord); 
    end
%% EEG = pop_eegfiltnew(EEG,highcut,0); Backup for filtering
%% Epoch
    EEG = pop_epoch(EEG, [triggers{:,2}],[epochtime]);
    %% Baseline
    EEG = pop_rmbase(EEG,[],[]);      
    save([writdir filesep f_name '.mat'],'EEG','-v7.3') 
end
end
