% This function does ERP-based decoding for orientation reported in Experiment 1.

% Structure of analysis and names of some variables used in this script

% were borrowed from the analysis script used in Foster et al. (2016).

% Note: Randomization routine included in this analysis (e.g., random permutation)

% can produce results slightly different from results reported in the paper. 

% Gi-Yeul Bae 2017.10.3.

function SVM_ECOC_ERP_Decoding(subs)
% delete(gcp)
% parpool
addpath("G:\My Drive\MudrikLab020818\Experiments_new\Jonathan\erp-decoding\software\eeglab2020_0")
if nargin ==0

    subs = [505 506 507 508 509 510 512 514 516 517 519 520 521 523 524 525];       
end

nSubs = length(subs);
% parameters to set

svmECOC.nChans = 16; % # of channels

svmECOC.nBins = svmECOC.nChans; % # of stimulus bins

svmECOC.nIter = 10; % # of iterations

svmECOC.nBlocks = 3; % # of blocks for cross-validation

svmECOC.frequencies = [0 6]; % low pass filter

svmECOC.time = -500:20:1496; % time points of interest in the analysis

svmECOC.window = 4; % 1 data point per 4 ms in the preprocessed data

svmECOC.Fs = 250; % samplring rate of in the preprocessed data for filtering

ReleventChan = sort([2,3,4,18,19, 5,6,20, 7,8,21, 9,10,11,12,13,14, 22,23,24,25,26,27, 15,16,1,17]); %electrodes

svmECOC.nElectrodes = length(ReleventChan); % # of electrode included in the analysis

% for brevity in analysis

nChans = svmECOC.nChans;

nBins = svmECOC.nBins;

nIter = svmECOC.nIter;

nBlocks = svmECOC.nBlocks;

freqs = svmECOC.frequencies;

times = svmECOC.time;
 
nElectrodes = svmECOC.nElectrodes;

nSamps = length(svmECOC.time);

Fs = svmECOC.Fs;

t = templateSVM('Verbose', 1)
%% Loop through participants
for cond = 1:1
    % 1: orientation
    
%% Loop through participants
for s = 1:nSubs
    sn = subs(s);

    fprintf('Subject:\t%d\n',sn)

    % load data
    currentSub = num2str(sn);
    dataLocation = strcat("G:\My Drive\MudrikLab020818\Experiments_new\Jonathan\erp-decoding\Exp1_Data");% set directory of data set
    loadThis = strcat(dataLocation,'\Decoding_DE_',currentSub,'.mat');
    load(loadThis);
    % data.eeg = data.eeg(data.channel == 1 | data.channel == 10);
    % where to save decoding output
    saveLocation = pwd; % set directory for decoding results.

    
    % set up locaiton bin of each trial

    channel = data.channel;
    %channel = data.channel(data.channel == 1 | data.channel == 10); % tip location of sample teardrop

    svmECOC.posBin = channel'; % add to fm structure so it's saved

    posBin = svmECOC.posBin;
    
    % grab EEG data
    eegs = data.eeg(:,ReleventChan,:); 
    
    % set up time points
    tois = ismember(data.time.pre:4:data.time.post,svmECOC.time);
    nTimes = length(tois);
    
    % # of trials
    svmECOC.nTrials = length(posBin); nTrials = svmECOC.nTrials;

    % Preallocate Matrices
    svm_predict = nan(nIter,nSamps,nBlocks,nChans); % a matrix to save prediction from SVM
    tst_target = nan(nIter,nSamps,nBlocks,nChans);  % a matrix to save true target values

    svmECOC.blocks = nan(nTrials,nIter);  % create svmECOC.block to save block assignments

    % low-pass filtering
    filtData = nan(nTrials,nElectrodes,nTimes);


 %   for c = 1:nElectrodes
%          filtData(:,c,:) = squeeze(eegs(:,c,:)); % low pass filter
  %  end

    size(filtData)
    parfor c = 1:nElectrodes
       filtData(:,c,:) = eegfilt(squeeze(eegs(:,c,:)),Fs,freqs(1,1),freqs(1,2)); % low pass filter
    end

    % Loop through each iteration
    tic % start timing iteration loop

    for iter = 1:nIter

        % preallocate arrays

        blocks = nan(size(posBin));

        shuffBlocks = nan(size(posBin));

        % count number of trials within each position bin

        clear binCnt

        for bin = 1:nBins

            binCnt(bin) = sum(posBin == bin); 

        end

        minCnt = min(binCnt); % # of trials for position bin with fewest trials

        nPerBin = floor(minCnt/nBlocks); % max # of trials such that the # of trials for each bin can be equated within each block

        % shuffle trials

        shuffInd = randperm(nTrials)'; % create shuffle index

        shuffBin = posBin(shuffInd); % shuffle trial order

        % take the 1st nPerBin x nBlocks trials for each position bin.

        for bin = 1:nBins;   

            idx = find(shuffBin == bin); % get index for trials belonging to the current bin

            idx = idx(1:nPerBin*nBlocks); % drop excess trials

            x = repmat((1:nBlocks)',nPerBin,1); shuffBlocks(idx) = x; % assign randomly order trials to blocks - actually assign block IDs to data points

        end


        % unshuffle block assignment

        blocks(shuffInd) = shuffBlocks;

        % save block assignment

        svmECOC.blocks(:,iter) = blocks; % block assignment

        svmECOC.nTrialsPerBlock = length(blocks(blocks == 1)); % # of trials per block

        % Average data for each position bin across blocks   

        posBins = 1:nBins;

        blockDat_filtData = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged & filtered EEG data with resampling at 50 Hz

        labels = nan(nBins*nBlocks,1);  % bin labels for averaged & filtered EEG data

        blockNum = nan(nBins*nBlocks,1); % block numbers for averaged & filtered EEG data

        bCnt = 1;

        for ii = 1:nBins

            for iii = 1:nBlocks

                blockDat_filtData(bCnt,:,:) = squeeze(mean(filtData(posBin==posBins(ii) & blocks==iii,:,tois),1));

                labels(bCnt) = ii;

                blockNum(bCnt) = iii;

                bCnt = bCnt+1;

            end

        end

        % Do SVM_ECOC at each time point
        parfor t = 1:nSamps

            % grab data for timepoint t

            toi = ismember(times,times(t)-svmECOC.window/2:times(t)+svmECOC.window/2);
            
            % average across time window of interest

            dataAtTimeT = squeeze(mean(blockDat_filtData(:,:,toi),3));  

            % Do SVM_ECOC for each block
            for i=1:nBlocks % loop through blocks, holding each out as the test set

                trnl = labels(blockNum~=i); % training labels

                tstl = labels(blockNum==i); % test labels

                trnD = dataAtTimeT(blockNum~=i,:);    % training data

                tstD = dataAtTimeT(blockNum==i,:);    % test data

                % SVM_ECOC
                mdl = fitcecoc(trnD,trnl, 'Coding','onevsall','Learners','SVM');   %train support vector mahcine
                LabelPredicted = predict(mdl, tstD);       % predict classes for new data
                
                svm_predict(iter,t,i,:) = LabelPredicted;  % save predicted labels
               
                tst_target(iter,t,i,:) = tstl;             % save true target labels

            end % end of block

        end % end of time points

    end % end of iteration

    toc % stop timing the iteration loop


    OutputfName = strcat(saveLocation,'\Orientation_Results_ERPbased_',currentSub,'.mat');
    
    svmECOC.targets = tst_target;
    svmECOC.modelPredict = svm_predict; 
    svmECOC.testResults =  svmECOC.targets == svmECOC.modelPredict;
    svmECOC.successRate =  mean(svmECOC.testResults, "all") * 100;
    svmECOC.nBlocks = nBlocks;

    save(OutputfName,'svmECOC','-v7.3');

end % end of subject loop
end % end of condition loop


for s = 1:numel(subs)
sub = subs(s);
varName = strcat('sub_',num2str(sub))
resultsFile = strcat(saveLocation,'\Orientation_Results_ERPbased_',num2str(sub),'.mat');
tmp =  load(resultsFile);
results.(varName) = tmp.svmECOC;
results.(varName).successRatesTime = mean(results.(varName).testResults, [1 3 4])
allResults(s) = results.(varName).successRate;
end

end
