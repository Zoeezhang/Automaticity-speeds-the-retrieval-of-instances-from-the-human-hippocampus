clc,clear
filepath = '..\data\After_prepro\'; %eeg data
fname = dir(fullfile(filepath,'*.mat'));
w = 2.3;

sub = 1;
for id = 1:length(fname)
    clear eegdat EEG 
    fprintf('Running subject %i of %i for removing ouliners...\n',id,length(fname))
    load([filepath fname(id).name]);
    subjID = fname(id).name(8:9);
    clear Q1 Q3 threshold ele_list tr_list
    Q1 = quantile(squeeze(mean(EEG.data,2)),0.25,2);
    Q3 = quantile(squeeze(mean(EEG.data,2)),0.75,2);
    threshold = Q3+w*(Q3-Q1);
    outline = zeros(EEG.nbchan,EEG.trials);
    for e = 1:EEG.nbchan
        clear data1 index
        data1 = squeeze(mean(EEG.data(e,:,:),2));
        index = find(data1 > threshold(e));
        outline(e,index) = 1;
    end
    ele_list = find(sum(outline,2)~=0);
    tr_list = find(sum(outline,1)~=0);
    ele_out{sub} = ele_list;
    tr_out{sub} = tr_list;
    ele_num(sub) = EEG.nbchan;
    per_ele(sub) = length(ele_out{sub})/EEG.nbchan;
    per_tr(sub)  = length(tr_out{sub})/EEG.trials;
    sub = sub+1;
end

save(['..\data\250Hz\load2pra\remove_trials.mat'],'tr_out','ele_out','ele_list','outline','ele_num','per_ele','per_tr');

%% kick out fake ripple
clc,clear
Fs = '1000';
rataAcq = 1000/str2double(Fs);
ele_conds = {'HPC','PFC','MTG'};
remove_time = -500:rataAcq:5200;

for a = 1:3
    filepath = ['..\data\After_prepro' Fs filesep]; %eeg data
    fname = dir(fullfile(filepath,'*.mat'));
    if strcmp(ele_conds{a},'HPC')
        load('..\electrodes\electrodes_hippo_greymatter(321).mat');
        ele_name = hippo_name;
    else
        load('..\electrodes\brain_allpp_index1130.mat');
        bb = find(strcmp(brain_name,ele_conds{a}));
        ele_name = brain_ele(bb,1:end);
    end
    clear tm_index out_dex EEG
    for id = 1:length(fname)
        clear eegdat EEG zdata count_out number ele_list tr_list EEGb EEGc ele_index all_ele
        fprintf('Running subject %i of %i for TF processing...\n',id,length(fname))
        load([filepath fname(id).name]);
        rataAcq = 1000/round(EEG.srate);

        all_ele = find(ismember({EEG.chanlocs.labels},ele_name));
        ele_index = setdiff(all_ele,ele_out{id});
        EEG = pop_select(EEG,'channel',ele_index);  
        hpFilt = designfilt('highpassfir','PassbandFrequency', 250,'SampleRate', EEG.srate,'StopbandFrequency',200);
        gradient_1st = diff(EEG.data,1,2);
        zdata = zscore(gradient_1st,1,2); 
        gra_list = zeros(EEGc.nbchan,EEGc.trials);
        high_list = zeros(EEGc.nbchan,EEGc.trials);
        
        tm_index = dsearchn(EEG.times',remove_time');
                  
        for e = 1:length(ele_index)
            for tr = 1:EEGc.trials
                clear out_gra out_high gra_list high_list all_index
                out_dex{id}{e,tr} = [];
                EEGa = filtfilt(hpFilt,double(squeeze(EEG.data(e,:,tr))));
                amp = zscore(EEGa,1,2);
                out_gra = find(zdata(e,tm_index,tr) > 5);  
                out_high = find(amp(tm_index) > 5);
                gra_list = find(sum(out_gra~=0));
                high_list= find(sum(out_high~=0));
                all_index = union(out_gra,out_high);
                out_dex{id}{e,tr} = all_index;
            end
        end   
    end
    save(['..\data\ripple\remove_ripple_' ele_conds{a} '_' Fs 'Hz.mat'],'all_out','ele_out','out_dex','tm_index','remove_time','hpFilt');
end

