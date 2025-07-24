%% Add EEGLAB packages only when you first analyze the data, and remember to save the path.
% eeglab_path   = ['D:' filesep 'Toolbox' filesep 'eeglab2019']; % Folder for eeglab
% run([eeglab_path filesep 'eeglab.m']);
% close all;

%% Add fieldtrip packages only when you need it! better not to save the path for fieldtrip.
% ft_path     = ['D:' filesep 'Toolbox' filesep 'fieldtrip-20170704'];% Folder for fieldtrip
% addpath(ft_path)
% ft_defaults;

%% Parameters for preprocessing
clc,clear
cfg = [];
cfg.writdir     = '..\data\After_prepro\'; % Folder name for preprocessed data
cfg.dc_chans    = {'POL DC09';'POL DC10';'POL DC11'};% Names for channels sending triggers
cfg.cor_10      = 'no';% 5 ms correction for triggers if necessary, dont change it before you understand it
cfg.res_rate    = 250; % resample rate
cfg.filt_type   = 'eeglab'; % 'fieldtrip' or 'eeglab' for bandstop filtering.
cfg.bandcut     = [0.3];  % low-pass 350, high-pass 0.3（高于0.3通过）
cfg.triggers    = {'display' {'1'};}; % Triggers
cfg.epochtime   = [-1.5 7.7];  % Time window for single epoch
cfg.irrelevant_channels = {'POL DC01';'POL DC02';'POL DC03';'POL DC04';'POL DC05';'POL DC06';'POL DC07';'POL DC08';
                        'POL DC09';'POL DC10';'POL DC11';'POL DC12';'POL DC13';'POL DC14';'POL DC15';'POL DC16';
                        'POL $A11';'POL $A12';'POL $A''11';'POL $A''12';'POL Pulse';'POL SpO2';
                        'POL BP1';'POL BP2'; 'POL BP3'; 'POL BP4';'POL $B''11';'POL $B''12';
                        'POL $C''11';'POL $C''12';'POL CO2Wave';'POL $J''11';'POL $J''12';
                        'POL';'POL E';'POL EKG';'POL EKG1';'POL EKG2';'EDF Annotations';'POL EtCO2';
                        'POL LDELT';'POL LDELT1';'POL LDELT2';'POL LEXT1';'POL LEXT2';'POL LFLEX1';'POL LFLEX2';
                        'POL LTOE1';'POL LTOE2';'POL LTHE1';'POL LTHE2';'POL LTA1';'POL LTA2';                        
                        'POL RDELT';'POL REXT1';'POL REXT2';'POL RFLEX1';'POL RFLEX2';'POL RTA1';'POL RTA2';
                        'POL RTOE1';'POL RTOE2';'POL RTHE1';'POL RTHE2'; 'POL RDELT1'; 'POL RDELT2';'POL ECG-1';'POL ECG-2'};

%% Only set those parameters when you need to run single subject!
cfg.single_subj = 'no'; % If 'yes', it will run for one subject only; If 'no', it will run for all subjects.
cfg.single_name = 'subject21_pra.edf'; % If single_subj = yes, put the name here for running.

%% Load core code and run preprocessing
SEEG_preproc_250(cfg);