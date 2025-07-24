clc,clear
standard = 'voltage';
ele_conds = {'HPC','PFC','MTG'};
save_names = {'Z_correction_blk123','Z_correction_blk456'};

Fs = 1000; fre = [70 180]; 
filepath = ['..\data\After_prepro' num2str(Fs) filesep]; %eeg data
fname = dir(fullfile(filepath,'*.mat'));
rateAcq = 1000/Fs;
time = 3200:rateAcq:5200;
duration_time = [15 25]; %ms
nsample = duration_time/rateAcq;
min_dura = [num2str(duration_time(1)) 'ms'];
band = [num2str(fre(1)) '-' num2str(fre(end)) 'Hz']; 
ri_time = [num2str(time(1)) '-' num2str(time(end))];
order = 20;
re_cond = '-500-5200';
load('..\data\250Hz\load2pra\remove_trials.mat');

%% Start ripple
for a = 1:3
    clear out_dex edges
    load(['..\data\ripple\remove_ripple_' ele_conds{a} '_' num2str(Fs) 'Hz.mat'],'out_dex','remove_time'); % attention to time!!!
    if strcmp(ele_conds{a},'HPC')
        load('..\electrodes\electrodes_hippo_greymatter(321).mat','hippo_name');
        ele_name = hippo_name;
    else
        load('..\electrodes\brain_allpp_index1130.mat','brain_name','brain_ele');
        bb = find(strcmp(brain_name,ele_conds{a}));
        ele_name = brain_ele(bb,1:end);
    end
    ids = length(fname)-sum((cellfun(@isempty,ele_name)));
    for b = 1:2
        clear draw tm
        root = ['..\data\MVPA\' ele_conds{a}  filesep save_names{b} filesep];
        load([root '\per\noInf\' save_names{b} '-4-8Hz-tt-smooth.mat'],'tm','draw');
        edges{b} = [tm(draw(1)) tm(draw(end))];
    end

    SD = 2.5;
    ripples = cell(2,ids);      max_ripple = cell(2,ids);
    ele_rate_1st = cell(2,ids); ele_rate_med = cell(2,ids);
    rate_med = zeros(2,ids);    rate_1st = zeros(2,ids);
    out_ri=cell(1,ids);         subjID = {ids};
    sub = 0; th = [SD 3];
    tic 
    for id = 1:length(fname)
        if ~isempty(ele_name{id})
            clear EEG ele_index all_ele
            fprintf('Running subject %i of %i for TF processing...\n',id,length(fname))
            sub = sub+1;
            if id ~= 13  
                load([filepath fname(id).name]);
            else
                load([filepath '\..\After_prepro500\' fname(id).name]);
            end
            subjID{sub} = fname(id).name(8:9);
            Fs = EEG.srate;
            rateAcq = 1000/Fs;
            time = 3200:rateAcq:5200;   
            epoch_index = dsearchn(EEG.times',time'); % Define epoch
            all_ele = find(ismember({EEG.chanlocs.labels},ele_name{id}));
            ele_index = setdiff(all_ele,ele_out{id});  
            lpFilt = designfilt('bandpassfir','CutoffFrequency1', fre(1), 'CutoffFrequency2', fre(2),'SampleRate',EEG.srate,'FilterOrder',order);

            for b = 1:2
                clear EEGb EEGc EEGd tr_index
                tr_index = setdiff(1+(b-1)*90:90+(b-1)*90,tr_out{id});                                          
                                            
                for e = 1:length(ele_index)
                    clear hil_envelope data_2std data_3std 
                    EEGa = squeeze(double(EEG.data(ele_index(e),epoch_index,tr_index))); 
                    EEGb = filtfilt(lpFilt,EEGa);
                    hil_envelope = EEGb;
                    data_2std = mean(hil_envelope,1) + std(hil_envelope,0,1)*th(1);
                    data_3std = mean(hil_envelope,1) + std(hil_envelope,0,1)*th(2);
                    data_2std_neg = mean(hil_envelope,1) - std(hil_envelope,0,1)*th(1);
                    data_3std_neg = mean(hil_envelope,1) - std(hil_envelope,0,1)*th(2);
                    hil_envelope(hil_envelope>data_2std_neg & hil_envelope<data_2std)=0;

                    out_ri{id}{b,e} = []; count_r1=[]; count_r2=[];tr_ri_med=[];tr_ri_1st=[];
                    for tr = 1:size(hil_envelope,2)
                        clear amp_tr islands join_ripple toMerge iri clustNs dura_index 
                        amp_tr = hil_envelope(:,tr);
                        islands = bwconncomp(logical(amp_tr));
                        if ~isempty(islands.PixelIdxList)
                            for i = 1:length(islands.PixelIdxList)
                                join_ripple(i,1) = islands.PixelIdxList{i}(1);
                                join_ripple(i,2) = islands.PixelIdxList{i}(end);
                            end
                            % Merge ripples if inter-ripple period is too short (unless this would yield too long a ripple)
                            iri = join_ripple(2:end,1) - join_ripple(1:end-1,2);
                            toMerge = iri<=nsample(1); %% changed here
                            while any(toMerge)
                                % Get indices of first ripples in pairs to be merged
                                rippleStart = strfind([0 toMerge'],[0 1])';
                                % Incorporate second ripple into first in all pairs
                                rippleEnd = rippleStart+1;
                                join_ripple(rippleStart,2) = join_ripple(rippleEnd,2);
                                % Remove second ripples and loop
                                join_ripple(rippleEnd,:) = [];
                                iri = join_ripple(2:end,1) - join_ripple(1:end-1,2);
                                toMerge = iri<=nsample(1) ; %% changed here
                            end

                            clustNs = join_ripple(:,2)-join_ripple(:,1)+1;
                            dura_index = find(clustNs<nsample(2));
                            join_ripple(dura_index,:)=[];

                            fit_ripple =[]; 
                            for j = 1:size(join_ripple,1)
                                clear power
                                power_max = max(amp_tr(join_ripple(j,1):join_ripple(j,2)));
                                power_min = min(amp_tr(join_ripple(j,1):join_ripple(j,2)));
                                if power_max >= data_3std(tr) || power_min <= data_3std_neg(tr)
                                    fit_ripple = [fit_ripple;join_ripple(j,:)];
                                end
                            end
    
                            n = 1;                             
                            outline_time = []; last_ripple = [];clear artifacts
                            tr_com = intersect(1:90,tr_out{id});  
                            artifacts = out_dex{id}{e,tr+(90-length(tr_com))*(b-1)}; 
                            if ~isempty(artifacts) 
                                
                                for q = 1:length(artifacts)
                                    outline_time = [outline_time,artifacts(q)-100/rateAcq:artifacts(q)+100/rateAcq;];
                                end 
                                for j = 1:size(fit_ripple,1)
                                    clear rr_dex                           
                                    rr_dex = ismember(fit_ripple(j,1):fit_ripple(j,2),outline_time);
                                    if sum(rr_dex) == 0
                                        last_ripple = [last_ripple;fit_ripple(j,:)];
                                    else
                                        out_ri{id}{b,e}(n) = tr;
                                        n = n+1;
                                    end 
                                end
                            else
                                last_ripple = fit_ripple;
                            end

                            tr_ri_med = [];tr_ri_1st = [];
                            if ~isempty(last_ripple)
                                for j = 1:size(last_ripple,1)
                                    clear pp index
                                    [pp,index]= max(amp_tr(last_ripple(j,1):last_ripple(j,2)));
                                    max_ripple{b,sub}{e,tr}(j) = time(last_ripple(j,1)+index-1);
                                end
                                ripples{b,sub}{e,tr} = time(last_ripple); 
                                tr_ri_med = (time(last_ripple(:,1))+time(last_ripple(:,2)))/2;
                                tr_ri_1st = time(last_ripple(:,1));
                            end
                        end
                        count_r1 = [count_r1,tr_ri_med];
                        count_r2 = [count_r2,tr_ri_1st];
                    end
                    ele_rate_med{b,sub}(e) = histcounts(count_r1,edges{b})/(tr*(edges{b}(2)-edges{b}(1)));
                    ele_rate_1st{b,sub}(e) = histcounts(count_r2,edges{b})/(tr*(edges{b}(2)-edges{b}(1)));
                end
                rate_med(b,sub) = mean(ele_rate_med{b,sub});
                rate_1st(b,sub) = mean(ele_rate_1st{b,sub});
            end
        end
    end
    root1 = ['..\data\ripple\'];
    if ~exist(root1,'dir')
        mkdir(root1)
    end
    save([root1 '\ripple_' ele_conds{a} '_' ri_time '_' band '_2secs'],...
         'ripples','subjID','fre','duration_time','ele_name','time','max_ripple','amp_tr','hil_envelope','out_ri','nsample',...
         'min_dura','ele_rate_med','ele_rate_1st','rate_med','rate_1st','order','Fs','rateAcq','th','edges','re_cond');
    toc
end

