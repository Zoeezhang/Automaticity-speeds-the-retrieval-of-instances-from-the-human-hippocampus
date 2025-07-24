clc,clear
tm = 3200:4:5200;
tt = length(tm);
fre = [4 8];
band = '4-8Hz';
ele_cond ='PFC';
time_cond = '3200-5200';
Nitr = 10; 
NBins = 2;
nBlocks = 3;
simulations = [1000,1];
save_names = {'Z_correction_blk123','Z_correction_blk456'};

%% each trial contrast with prediction
 for conds = 1
    clear AverageAccuracy simulationT draw
    filepath = ['..\data\MVPA\' ele_cond filesep save_names{conds} filesep ]; 
    fname=dir(fullfile(filepath,'*.mat'));
    nSubjects = length(fname);
    AverageAccuracy = zeros(nSubjects,tt);
    for s = 1:2
        for siml = 1:simulations(s)
            fprintf('permutation:\t%d\n',siml)
            for id = 1:nSubjects      
                load([filepath fname(id).name]);
                DecodingAccuracy = zeros(tt,nBlocks,Nitr);
                svm_predict = svmECOC.modelPredict; 
                svm_target  = svmECOC.targets;
                for tp = 1:tt
                    for block = 1:nBlocks
                        for iter = 1:Nitr
                            RandomAnswer = shuffle(1:2); 
                            prediction = squeeze(svm_predict(iter,tp,block,:)); % this is predictions from models
                            target = squeeze(svm_target(iter,tp,block,:));
                            if s == 1
                                value = RandomAnswer';
                            else
                                value = target;
                            end
                            Dist = value - prediction;
                            ACC = mean(Dist==0); 
                            DecodingAccuracy(tp,block,iter,:) = ACC; 
                        end
                    end
                end
                %% now compute grand average for a participant
                grandAvg = squeeze(mean(mean(DecodingAccuracy,3),2)); 
                smoothedAvg = nan(1,tt);
                for tAvg = 1:tt
                     if tAvg ==1
                       smoothedAvg(tAvg) = mean(grandAvg((tAvg):(tAvg+2)));
                     elseif tAvg == 2
                       smoothedAvg(tAvg) = mean(grandAvg((tAvg-1):(tAvg+2)));
                     elseif tAvg == (tt-1)
                       smoothedAvg(tAvg) = mean(grandAvg((tAvg-2):(tAvg+1)));
                     elseif tAvg == tt
                       smoothedAvg(tAvg) = mean(grandAvg((tAvg-2):(tAvg)));
                     else
                       smoothedAvg(tAvg) = mean(grandAvg((tAvg-2):(tAvg+2)));  
                     end
                end
                AverageAccuracy(id,:) =smoothedAvg; % average across iteration and block
            end   
             %% now compute average accuracy across participants
             subAverage = squeeze(mean(AverageAccuracy,1)); 
             seAverage = squeeze(std(AverageAccuracy,1))/sqrt(nSubjects);   
            %% do cluster mass analyses
            Ps = nan(2,tt);
            for tp = 1:tt % make sure this time range is correct
                [H,P,CI,STATS] =  ttest(AverageAccuracy(:,tp),0.5,'tail','right'); 
                Ps(1,tp) = STATS.tstat; 
                Ps(2,tp) = P; 
            end
            % find significant points
            candid = Ps(2,:) <= .05 & Ps(2,:)~= 0;
            candid_marked = zeros(1,length(candid));
            candid_marked(1,1) = candid(1,1);
            candid_marked(1,length(candid)) = candid(1,length(candid));
            %remove orphan time points
            for i = 2:tt-1
                if candid(1,i-1) == 0 && candid(1,i) ==1 && candid(1,i+1) ==0 
                    candid_marked(1,i) = 0; 
                else
                    candid_marked(1,i) = candid(1,i);     
                end
            end
            % combine whole time range with relevent time & significant information
            clusters = zeros(tt,1); 
            clusterT = zeros(tt,1);
            clusters(:,1) = candid_marked;
            clusterT(:,1) = Ps(1,:);
            clusterTsum = sum(Ps(1,logical(candid_marked)));
            %find how many clusters are there, and compute summed T of each cluster
            tmp = zeros(10,25);
            cl = 1;
            member = 1;
            for i = 2:length(clusters)-1
                if clusters(i-1) ==0 && clusters(i) == 1 && clusters(i+1) == 1 
                    cl = cl+1;
                    member = member +1;
                    tmp(cl,member) = i;           
                elseif clusters(i-1) ==1 && clusters(i) == 1 && clusters(i+1) == 0 
                    member = member +1;  
                    tmp(cl,member) = i;    
                    member = 0;  
                elseif clusters(i-1) ==1 && clusters(i) == 1 && clusters(i+1) == 1             
                    member = member +1;  
                    tmp(cl,member) = i;           
                else                    
                end
            end        
            HowManyClusters = cl;
            a = tmp(1:cl,:);
            eachCluster = a(:,logical(sum(a,1)~=0));      
            %% now, compute summed T of each cluster 
            dat_clusterSumT = zeros(HowManyClusters,1);

            for c = 1:HowManyClusters
               dat_clusterSumT(c,1) = sum(clusterT(eachCluster(c,eachCluster(c,:) ~=0)));
            end 
            if s == 1
                if size(dat_clusterSumT,1) > 0
                    simulationT(1,siml) = max(dat_clusterSumT);
                else
                    simulationT(1,siml) = 0;     
                end  
            end
        end % end of simulation
        if s == 1
            simulationT = sort(simulationT);
            root = [filepath '\per\noInf\'];
            if ~exist(root,'dir')
                mkdir(root)
            end
            save([root save_names{conds} '-' band '-per-smooth.mat'],'simulationT','tm','fre','nSubjects');
        else
            load([root save_names{conds} '-' band '-per-smooth.mat']);
            iteration = 1000;
            cutOff = iteration - iteration * 0.05; %one tailed
            % sortedTvlaues = sort(EmpclusterTvalue,2);
            critT = simulationT(cutOff); % 2 tailed iteration * 0.025
            sigCluster = dat_clusterSumT > critT;
            draw = eachCluster(sigCluster,:);
            draw = sort(reshape(draw,1,size(draw,1)*size(draw,2)));
            draw = draw(draw>0);
            
            for si = 1 :size(dat_clusterSumT,1)  
                [~,where] = min(abs(simulationT - dat_clusterSumT(si)));
                pv = 1 - where/iteration;
            end
            save([root save_names{conds} '-' band '-tt-smooth.mat'],'draw','seAverage','subAverage','tm','fre','nSubjects');
        end
    end
 end

