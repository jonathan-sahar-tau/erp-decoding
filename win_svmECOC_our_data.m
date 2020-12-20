% This function does ERP-based decoding for orientation reported in Experiment 1.

% Structure of analysis and names of some variables used in this script

% were borrowed from the analysis script used in Foster et al. (2016).

% Note: Randomization routine included in this analysis (e.g., random permutation)

% can produce results slightly different from results reported in the paper. 

% Gi-Yeul Bae 2017.10.3.

function svmECOC_our_data()
delete(gcp)
parpool
addpath("G:\My Drive\MudrikLab020818\Experiments_new\Jonathan\erp-decoding\software\eeglab2020_0")

if nargin == 0
    subjects = [11 13 14 15 16 17 19 20 21];
end
% subjects = [11 13] %14 15 16];
subjects
nSubs = length(subjects);

% parameters to set

svmECOC.conditions = ["cong_int", "cong_scr", "inc_int"];

svmECOC.nChans = numel(svmECOC.conditions); % # of channels

svmECOC.nBins = svmECOC.nChans; % # of stimulus bins

svmECOC.nIter = 10; % # of iterations

svmECOC.nBlocks = 5; % # of blocks for cross-validation

svmECOC.frequencies = [0 30]; % low pass filter

svmECOC.time = -99:2:898; % time points of interest in the analysis

svmECOC.window = 4; % 1 data point per 4 ms in the preprocessed data

svmECOC.Fs = 512 ; % sampling rate of the preprocessed data for filtering


svmECOC.nElectrodes = 64; % # of electrode included in the analysis

% for brevity in analysis

nChans = svmECOC.nChans;

nBins = svmECOC.nBins;

nIter = svmECOC.nIter;

nBlocks = svmECOC.nBlocks;
nAveragedTrials = nBlocks;

freqs = svmECOC.frequencies;

times = svmECOC.time;
 
nElectrodes = svmECOC.nElectrodes;

nSamps = length(svmECOC.time);

Fs = svmECOC.Fs;

conditions = svmECOC.conditions;

t = templateSVM('Verbose', 1);

fieldsToDelete = {'filename', 'filepath', 'subject', 'group', ...
                'condition', 'setname', 'session', 'srate', ...
                'xmin', 'xmax', 'icaact', 'icawinv', ...
                'icasphere', 'icaweights', 'icachansind', 'chanlocs', ...
                'urchanlocs', 'chaninfo', 'ref', 'event', ...
                'urevent', 'eventdescription', 'epoch', 'epochdescription', ...
                'reject', 'stats', 'specdata', 'specicaact', ...
                'splinefile', 'icasplinefile', 'dipfit', 'history', ...
                'saved', 'etc', 'run'};

addpath("G:\My Drive\MudrikLab020818\Experiments_new\Jonathan\erp-decoding\software\eeglab2020_0")

baseDir = "G:\My Drive\MudrikLab020818\Experiments_new\Jonathan\erp-decoding\";
dataLocation = strcat(baseDir, "exported data\");
exportDataLocation = strcat(baseDir, "experiment_data\");
outputDir = strcat(baseDir, "output\");

outputFile = strcat(outputDir, datestr(now,'mm-dd-yyyy-HH-MM'), '_run-log')

nCVBlocks = 3;
averagedSuccessRates = nan(numel(subjects), 1);
allSuccessRates = nan(nCVBlocks, numel(subjects));
%% Loop through participants
    for subjectIdx = 1:numel(subjects)
        subject = subjects(subjectIdx);
        subjectName = num2str(subject, '%03.f');
        fprintf('Subject:\t%d\n',subject);
        
        % assign an integer value to each condition, so we can later put all the data in
        % one matrix and reference it with this number.
        enumVal = 1;

        for condition = conditions
            exportFileName = strcat(subjectName, '_', condition, '_bc', '.mat');
            exportLocation = strcat(exportDataLocation, '\', exportFileName);
%           currentFilename= strcat(subjectName, '_', condition, '_bc', '.vhdr');
%           tempData = (pop_loadbv(dataLocation, currentFilename, [], 1:64));
%           save(exportLocation, 'tempData', '-v7.3')
 
            % remove uneeded fields from the struct
            load(exportLocation); %into variable tempData
            tempData = rmfield(tempData, fieldsToDelete);

            % swap the 1st and 2nd dimensions, effectively transposing the
            % "inner" (electrode) matrices to be nTrials x nSamples -
            % every row is the time course of one trial
            tempData.data = permute(tempData.data, [3, 1, 2]);
            EEG.(condition) = tempData;
            EEG.(condition).name = condition;
            EEG.(condition).enumVal = enumVal;

            % get the minimal number of trials across conditions
            trialNums(enumVal) = EEG.(condition).trials;
            enumVal = enumVal + 1;

        end

        % find the actual # of trials that can accomodate our desired # of averagedTrials,
        nTrialsPerAveragedTrial = floor(min(trialNums)/nAveragedTrials);
        nTrialsPerCondition = nTrialsPerAveragedTrial * nAveragedTrials;

        % normalize all conditions to have the same number of time points in each trial recording
        nSamples = EEG.(conditions(1)).pnts;
        for condition = conditions
            enumVal = EEG.(condition).enumVal;
            conditionData{enumVal} =  EEG.(condition).data(1:nTrialsPerCondition,:,:);
            conditionLabels{enumVal} = repmat(enumVal, nTrialsPerCondition, 1); 
        end
        
        allData = cat(1, conditionData{:});
        allLabels = cat(1, conditionLabels{:});


        data.eeg = allData;
        data.time.pre = -99;
        data.time.post = 898;

        % set up locaiton bin of each trial

        channel = allLabels;

        svmECOC.posBin = channel; % add to fm structure so it's saved

        posBin = svmECOC.posBin;

        % grab EEG data
        eegs = data.eeg;

        % set up time points
        tois = ismember(data.time.pre:2:data.time.post,svmECOC.time);
        nTimes = length(tois);

        % # of trials
        svmECOC.nTrials = length(posBin); nTrials = svmECOC.nTrials;

        % Preallocate Matrices
        svm_predict = nan(nIter,nSamps,nBlocks,nChans); % a matrix to save prediction from SVM
        tst_target = nan(nIter,nSamps,nBlocks,nChans);  % a matrix to save true target values

        svmECOC.blocks = nan(nTrials,nIter);  % create svmECOC.block to save block assignments

        % low-pass filtering
        filtData = nan(nTrials,nElectrodes,nTimes);

    % 
    %     for c = 1:nElectrodes
    %           filtData(:,c,:) = squeeze(eegs(:,c,:)); % low pass filter
    %     end

        size(filtData)

        tic
        parfor c = 1:nElectrodes
              tmp = eegfilt(squeeze(eegs(:,c,:)),Fs,freqs(1,1),freqs(1,2));
              filtData(:,c,:) = tmp(:,1:nTimes); % low pass filter
        end
        toc
    
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
            for t = 1:nSamps

                % grab data for timepoint t

                toi = ismember(times,times(t)-svmECOC.window/2:times(t)+svmECOC.window/2);

                % average across time window of interest

                dataAtTimeT = squeeze(mean(blockDat_filtData(:,:,toi),3));  

                % Do SVM_ECOC for each block
                fprintf("training...\n");
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


    OutputfName = strcat(outputDir, 'congruency_and_intactness_results_',subjectName,'.mat');
    
    svmECOC.nBlocks = nBlocks;
    svmECOC.targets = tst_target;
    svmECOC.modelPredict = svm_predict; 
    svmECOC.testResults =  svmECOC.targets == svmECOC.modelPredict;
    svmECOC.overallSuccessRatePcnt =  mean(svmECOC.testResults, "all") * 100;
    svmECOC.successRatesTime = mean(svmECOC.testResults, [1 3 4]);
    tmp = sort(svmECOC.successRatesTime, 'descend');
    svmECOC.topSuccessRates = tmp(1:10) 
%     allSubjectSuccessRatesTime{subjectName} = svmECOC.successRatesTime;

    save(OutputfName,'svmECOC','-v7.3');
    end % end of subject loop
   
    OutputfName = strcat(outputDir, 'congruency_and_intactness_results_all.mat');
    save(OutputfName,'allSubjectSuccessRatesTime','-v7.3');

    
% subjects = [11 13 14 15 16 17 19 20 21];    
% baseDir = "G:\My Drive\MudrikLab020818\Experiments_new\Jonathan\erp-decoding\";
% dataLocation = strcat(baseDir, "exported data\");
% exportDataLocation = strcat(baseDir, "experiment_data\");
% outputDir = strcat(baseDir, "output\");

% for s = 1:numel(subjects)
% sub = subjects(s);
% subjectName = num2str(sub, '%03.f');
% resultsFile = strcat(outputDir, 'congruency_and_intactness_results_',subjectName,'.mat');
% tmp =  load(resultsFile);
% results{s} = tmp.svmECOC.successRatesTime;
% end
end
