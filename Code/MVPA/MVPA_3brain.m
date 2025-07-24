%% MVPA decoding for orientation in load 2 condition.
clc,clear
filepath = '..\data\behavior\all\'; % behavior
fname=dir(fullfile(filepath,'*.mat'));

%% paramater
fre = [4 8]; 
band = '4-8Hz';
ele_cond = 'PFC';
time_cond = '3200-5200';
tm = 3200:4:5200;
tt = length(tm);
svmECOC.nBlocks = 3;
svmECOC.nBins = 2; % # of direction bin 
svmECOC.nIter = 10; % # of iterations 
svmECOC.frequencies = fre; % frequency bands to analyze  / low pass filter 
svmECOC.time = tm; % time points for re-sampling  
svmECOC.Fs = 250; % sampling rate of pre-processed data 
svmECOC.window = 1000/svmECOC.Fs; % 1 data point per 4 ms  
save_names = {'Z_correction_blk123','Z_correction_blk456'};
sub = {'01','03','05','07','08','10','11','12','13','15','16','17','18','19','20','21','22','23'};

%% load all subjects' tf_pow
dirs  = 'Z_meanfre_2blk-baseline_';
load(['..\load2pra\MVPA\' time_cond '\data\' dirs band '(' time_cond ').mat']);
load('..\data\250Hz\load2pra\remove_trials.mat');

%% define electrodes
if strcmp(ele_cond,'HPC')          
    load('..\electrodes\electrodes_hippo_greymatter(321).mat');
    ele_name = hippo_name;
else
    load('..\electrodes\brain_allpp_index1130.mat');
    bb = find(strcmp(brain_name,ele_cond));
    ele_name = brain_ele(bb,1:end);
end

%% training SVM and decoding
for id = 1:length(fname)
    fprintf('Subject:\t%d\n',id)
    aa = fname(id).name(7:8);
    if aa == sub{id}
        clear findele_index all_ele stim
        % load behavioral data
        load([filepath fname(id).name]);
        cue = reshape(permute(stim.cue,[2,1]),[1,180]);  % Bin numbers
        all_ele = find(ismember({subs_dim{id}.chans.labels},ele_name{id}));
        findele_index = setdiff(all_ele,ele_out{id});
        
        if ~isempty(findele_index) 
            for cond = 1:2
                clear pow eeg_data n_cue tr_cue
                tr_index = setdiff(1+(cond-1)*90:90+(cond-1)*90,tr_out{id})-(cond-1)*90;
                pow = Z_correction_meantrial_2blk{id,cond}(findele_index,:,tr_index);
                trial_num = size(pow,3);
                nBlocks = svmECOC.nBlocks;
                nBins = svmECOC.nBins;
                nElectrodes = size(pow,1);
                stimBins = 1:nBins;
                nIter = svmECOC.nIter;
                svmECOC.nElectrodes = findele_index;
                eeg_data = permute(pow,[3 1 2]);
                n_cue = cue(tr_index+(cond-1)*90);

                % Define the labels for different conditions                 
                tic
                svm_predict = nan(nIter,tt,nBins,2);
                tst_target = nan(nIter,tt,nBins,2);
                for iter = 1:nIter
                    fprintf('iter:\t%d\n',iter)
                    clear random_index binCnt minCnt nPerBin shuffBin random_block x tr_cue 
                    for bin = 1:nBins
                        binCnt(bin) = sum(n_cue == bin); 
                    end
                    minCnt = min(binCnt); % # of trials for position bin with fewest trials
                    nPerBin = floor(minCnt/nBlocks); % max # of trials such that the # of trials for each bin can be equated within each block
                    random_index = randperm(trial_num); % shuffle trials
                    shuffBin = n_cue(random_index); % shuffle trial order
                    shuffBlocks = zeros(1,trial_num);
                    for bin = 1:nBins 
                        idx = find(shuffBin == bin); % get index for trials belonging to the current bin
                        idx = idx(1:nPerBin*nBlocks); % drop excess trials
                        x = repmat((1:nBlocks)',nPerBin,1);
                        shuffBlocks(idx) = x; % assign randomly order trials to blocks
                    end
                    
                    random_block(random_index) = shuffBlocks; % unshuffle block assignment
                    tr_cue = n_cue;
                    tr_cue(random_block==0)= [];
                    random_block(random_block==0)=[];
                    
                    blockDat_filtData = nan(nBins*nBlocks,nElectrodes,tt);  % averaged & filtered EEG data
                    labels = nan(nBins*nBlocks,1);                              % bin labels for averaged & filtered EEG data
                    blockNum = nan(nBins*nBlocks,1);                            % block numbers for averaged & filtered EEG data 
                    bCnt = 1;
                    for ii = 1:nBins
                        for iii = 1:nBlocks
                            blockDat_filtData(bCnt,:,:) = squeeze(mean(eeg_data(tr_cue==stimBins(ii) & random_block==iii,:,:),1)); % 这里按照block和类别先把trial提取出来进行平均，然后在进行训练
                            labels(bCnt) = ii;
                            blockNum(bCnt) = iii;
                            bCnt = bCnt+1;
                        end
                    end

                    % decoding over time
                    for t = 1:tt
                        dataAtTimeT = squeeze(mean(blockDat_filtData(:,:,t),3));  
                        for i = 1:nBlocks % loop through blocks, holding each out as the test set
                            trnl = labels(blockNum~=i); % training labels
                            tstl = labels(blockNum==i); % test labels
                            trnD = dataAtTimeT(blockNum~=i,:);    % training data
                            tstD = dataAtTimeT(blockNum==i,:);    % test data
                            % here define how to perform the learning
                            mdl = fitcecoc(trnD,trnl, 'Coding','onevsall','Learners','SVM' );   %train support vector mahcine
                            LabelPredicted = predict(mdl, tstD);  % predict class of new data
                            svm_predict(iter,t,i,:) = LabelPredicted;  % save predicted labels
                            tst_target(iter,t,i,:) = tstl;             % save true target labels                  
                        end
                    end  
                 end
                
                svmECOC.targets = tst_target;   
                svmECOC.modelPredict = svm_predict;
                
                toc           
                root = ['..\data\MVPA\' ele_cond filesep save_names{cond} filesep];
                if ~exist(root,'dir')
                    mkdir(root)
                end
                data_root = [root '\pp' aa '-SVM.mat']; 
                save(data_root,'svmECOC','nBlocks','n_cue','tr_index','tr_cue','-v7.3');
            end
        end
    else
        fprintf('Attention:\t%d\n',id)   
    end
end

