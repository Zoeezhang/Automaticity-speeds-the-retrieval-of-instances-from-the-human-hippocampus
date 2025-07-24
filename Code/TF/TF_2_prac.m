%% SEEG_data,Time-frequency analysis
% revised by BW 11.30.2021
clc,clear all;

%% Only set those parameters when you need to run single subject!
single_subj = 'no'; % If 'yes', it will run for one subject only; If 'no', it will run for all subjects.
single_name = 'pp_subject02_pra.mat'; % If single_subj = yes, put the name here for running.

%% Add TF toolbox dir
tf_dir = ['..\code\core' filesep 'tfdecomp-master']; % Folder for tf decomp
addpath(tf_dir)
%% Load data
% Getting file path for storing data, by choose the .mat file manually.
[~, filepath]=uigetfile('*.mat'); 
% Getting data files
files = dir(fullfile(filepath,'*.mat'));
% Folder for storing preprocessed data
write_clean_dir = ['..\data\250Hz\load2pra\TF\TF_power\subjects\'];
if ~exist(write_clean_dir,'dir')
    mkdir(write_clean_dir)
end
% Define the number of participants for preprocessing
if strcmp(single_subj,'yes')
    sub_num = 1;
    dataset= [filepath single_name];
else
    sub_num = length(files);
end

%% Start TF processing
for id = 56%sub_num
    clear tf_pow tf_phase dim eegdat tmprej; 
    fprintf('Running subject %i of %i for TF processing...\n',id,sub_num)
    %clear previous dataset
    clear EEG
    %load dataset when preprocessing multiple datasets
    if strcmp(single_subj,'no')
        dataset= [filepath files(id).name];
    end
    load(dataset);
    if EEG.trials == 180
        %% Parameters for TF processing  
        for i = 1:EEG.nbchan
            fprintf('Running channel %i of %i for TF processing...\n',i,EEG.nbchan)
            EEGb = pop_select(EEG,'channel',[i]);
            cfg = [];
            cfg.channels     = 1:EEGb.nbchan;  %根据被试改导数channel
            cfg.chanlocs     = EEGb.chanlocs;  
            cfg.frequencies  = [1:30 35:5:120]; % from min to max in nsteps 
            cfg.cycles       = [4 20]; % min max number of cycles used for min max frequency 
            cfg.scale        = 'lin'; % whether the above frequency settings will be logarithmically (recommended) or linearly scaled ('lin') 
            cfg.basetime     = [-500 -300];% pre-stim baseline
            cfg.baselinetype = 'conavg'; % Baseline correction type: 'conavg' or 'conspec'
            cfg.times2save   = -1500:4:7700;% Time window for saving data
            cfg.erpsubtract  = false; % if true, non-phase-locked (i.e. "induced") power will be computed by subtracting the ERP from each single trial
            cfg.matchtrialn  = false; % if true, conditions will be equated in terms of trial count (so SNR is comparable across conditions)
            cfg.srate        = EEGb.srate; % sample rate
            cfg.eegtime      = EEGb.times;  %epoch time
            cfg.singletrial  = true; 
            
            cfg.report_progress = true;
            cfg.save_output = true;
            cfg.overwrite = false;
    
            cfg.writdir = write_clean_dir;
            if strcmp(single_subj,'yes')
                f_name = single_name;
            else
                f_name = files(id).name;
            end
    
            cfg.filename = ['tf_' f_name];          
            eegdat{1} = EEGb.data;
           
            [tf_pow_ele,tf_phase_ele,dim_ele] = tfdecomp_raw_power(cfg,eegdat); 
            
            if cfg.singletrial
                tf_pow(i,:,:,:)= tf_pow_ele{1};
                tf_phase = NaN;
                name = 'trial-';
            else 
                tf_pow(:,i,:,:)=tf_pow_ele;
                tf_phase(:,i,:,:)=tf_phase_ele;
                name = 'block';
            end
            
            if i ==1
                dim = dim_ele;
            end
            dim.chans(i,1)=dim_ele.chans;
            clear tf_pow_ele tf_phase_ele dim_ele;  
        end
        save([write_clean_dir filesep [name f_name]], 'tf_pow', 'dim','-v7.3')
    else
        fprintf('Subject %i of %i trials are not enough...\n',id,sub_num)
    end
end

