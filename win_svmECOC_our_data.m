% This function does ERP-based decoding for orientation reported in Experiment 1.

% Structure of analysis and names of some variables used in this script

% were borrowed from the analysis script used in Foster et al. (2016).

% Note: Randomization routine included in this analysis (e.g., random permutation)

% can produce results slightly different from results reported in the paper. 

% Gi-Yeul Bae 2017.10.3.


function svmECOC_our_data(subjects)
% delete(gcp)
% parpool

if nargin == 0
%         subjects = [102 104:106 108 109]
    %     subjects = [201] %:206]
%         subjects = [11 13:17 19:21];
%     subjects = [22 24:27 29 30 35:42]
                subjects = [102 104:106 108:112 114:116 118:120 122]


end

    C = Constants();
    subjects = C.subjects;
    nSubjects = length(subjects);
    debugAddBump = 0
    debugApplyLowpassFilter = 1
    % parameters to set

    allElectrodesInOrder = C.allElectrodesInOrder;

    frontalElectrodes = C.frontalElectrodes;

    centralElectrodes = C.centralElectrodes;

    occipitalElectrodes = C.occipitalElectrodes;

    combinedElectrodes = union(centralElectrodes, occipitalElectrodes);
    
    conditionDesc = "Congruency";

    conditions = ["ConInt", "IncInt"] %, "ConScr", "IncScr"]
    labels = [1, 2] %, 1, 2]; % intact = 1, scrambled = 2


    
    nCVBlocks = 3; % # of blocks for cross-validation
    
    frequencies = [0 30]; % low pass filter
    
    window = 4; % 1 data point per 4 ms in the preprocessed data
    
    Fs = 512 ; % sampling rate of the preprocessed data for filtering
    
    
    % for brevity in analysis
    
    nIter = C.nIter;

    nCVBlocks = C.nCVBlocks;
    nAveragedTrials = C.nAveragedTrials;
    
    freqs = C.frequencies;
    
    Fs = C.Fs;
    
    svmTemplate = C.svmTemplate;
    
    fieldsToDelete = C.fieldsToDelete;

    baseDir = C.baseDir;
    eeglabPath = C.eeglabPath;
    dataLocation = C.dataLocation;
    outputDir = C.outputDir;
    resultsDir = C.resultsDir;
    figuresDir = C.figuresDir;
    
    averagedSuccessRates = nan(numel(subjects), 1);
    allSuccessRates = nan(nCVBlocks, numel(subjects));
    
    % set the conditions and the most relevant electrode for those conditions
    % ----------------------------

    relevantElectrodes = centralElectrodes; 

    nConditions = numel(conditions); % # of channels
    nBins = numel(unique(labels)); % # of stimulus bins

    nElectrodes = numel(relevantElectrodes); % # of electrode included in the analysis
    conditions = conditions;
    nChans = nConditions;
    nBins = nBins;
    
    analysisTic = tic;
    %% Loop through participants
    for subjectIdx = 1:numel(subjects)
        subjectTic = tic;
        subject = subjects(subjectIdx);
        subjectName = num2str(subject, '%03.f');
        fprintf('Subject:\t%d\n',subject);
        
        
        for conditionIdx  = 1:numel(conditions)
            condition = conditions(conditionIdx);
            fileName = strcat(subjectName, '_', condition, '_bc', '.mat')
            % fileName = strcat('S', subjectName, '_N300Exp1', '.set');
            filePath = strcat(dataLocation, '\', fileName);
            
            % remove unneeded fields from the struct
            tempData = load(filePath); % into variable tempData
            tempData = tempData.EEGData;
            x = tempData.chanlocs.labels;
            % tempData = rmfield(tempData, fieldsToDelete);
            
            % swap the the matrix dimensions so that they're
            % nTrials x nElectrodes x nSamples
            tempData.data = permute(tempData.data, [3, 1, 2]);
            if debugAddBump == 1
                if conditionIdx == 2
                    bump = zeros(size(tempData.data));
                    electrodeMax = 100 * max(squeeze(mean(tempData.data, 1)), [], 2);
                    bumpTimes = (250 < tempData.times  &  tempData.times < 300);
                    bump(:,:,bumpTimes) = repmat(electrodeMax', size(bump, 1), 1, sum(bumpTimes));
                    tempData.data =  tempData.data + bump;
                end
                
                if conditionIdx == 3
                    bump = zeros(size(tempData.data));
                    electrodeMax = 100 * max(squeeze(mean(tempData.data, 1)), [], 2);
                    bumpTimes = (450 < tempData.times  &  tempData.times < 500);
                    bump(:,:,bumpTimes) = repmat(electrodeMax', size(bump, 1), 1, sum(bumpTimes));
                    tempData.data =  tempData.data - bump;
                end
            end
            EEG.(condition) = tempData;
            EEG.(condition).name = condition;

            % get the minimal number of trials across conditions
            trialNums(conditionIdx) = EEG.(condition).trials;
        end
        
        % find the actual # of trials that can accomodate our desired # of averagedTrials,
        nTrialsPerAveragedTrial = floor(min(trialNums)/nAveragedTrials);
        nTrialsPerCondition = nTrialsPerAveragedTrial * nAveragedTrials;
        % normalize all conditions to have the same number of time points in each trial recording
        for conditionIdx  = 1:numel(conditions)
            condition = conditions(conditionIdx);
            electrodeIdx = ismember(allElectrodesInOrder, relevantElectrodes);
            conditionData{conditionIdx} =  ...
                EEG.(condition).data(1:nTrialsPerCondition,...
                electrodeIdx, ...
                :);
            conditionLabels{conditionIdx} = repmat(labels(conditionIdx), nTrialsPerCondition, 1);
        end
        
        clear labels;
        
        allData = cat(1, conditionData{:});
        allLabels = cat(1, conditionLabels{:});
        
        data.eeg = allData;
        data.labels = allLabels;
        data.times = floor(EEG.(conditions(1)).times);
        data.time.pre = data.times(1);
        data.time.post = data.times(end);
  
        if debugAddBump == 1
            figure
            subplot(2, 2, 1);
            plot(data.times, squeeze(allData(nTrialsPerCondition+1,1,:)))
            subplot(2, 2, 2);
            plot(data.times, squeeze(conditionData{2}(1,1,:)))
        end

        resamplingRatio = 1000/250; % resample the data during analysis to 250 Hz
        downsampledTimes = data.time.pre:resamplingRatio:data.time.post; % time points of interest in the analysis
        nSamps = length(downsampledTimes);
        
        % set up locaiton bin of each trial
        posBin = data.labels;
        posBin = posBin;
        
        eegs = data.eeg;
        
        % set up time points
        tois = ismember(data.time.pre:2:data.time.post,downsampledTimes);
        nTimes = length(tois);
        
        % # of trials
        nTrials = length(posBin);
        nTrials = nTrials;
        
        % Preallocate Matrices
        svm_predict = nan(nIter,nSamps,nCVBlocks,nBins); % a matrix to save prediction from SVM
        tst_target = nan(nIter,nSamps,nCVBlocks,nBins);  % a matrix to save true target values
        
        blocksAssign = nan(nTrials,nIter);  % create block to save block assignments
        
        % low-pass filtering
        filtData = nan(nTrials,nElectrodes,nTimes);

        if debugApplyLowpassFilter == 1
            tic
            for c = 1:nElectrodes
                    tmp = eegfilt(squeeze(eegs(:,c,:)),Fs,freqs(1,1),freqs(1,2)); % low pass filter
                    filtData(:,c,:) = tmp(:,1:nTimes);
            end
            toc
        else
            filtData = eegs(:,:,1:nTimes);
        end

        
        % Loop through each iteration
        for iter = 1:nIter
            fprintf("training iteration %d\n", iter);
            iterTic = tic; % start timing iteration loop
            
            % preallocate arrays
            blocks = nan(size(posBin));
            
            shuffBlocks = nan(size(posBin));
            
            % count number of trials within each position bin
            
            clear binCnt
            
            for bin = 1:nBins
                
                binCnt(bin) = sum(posBin == bin);
                
            end
            
            minCnt = min(binCnt); % # of trials for position bin with fewest trials
            
            nPerBin = floor(minCnt/nCVBlocks); % max # of trials such that the # of trials for each bin can be equated within each block
            
            % shuffle trials
            
            shuffInd = randperm(nTrials)'; % create shuffle index
            
            shuffBin = posBin(shuffInd); % shuffle trial order
            
            % take the 1st nPerBin x nCVBlocks trials for each position bin.
            
            for bin = 1:nBins;
                
                idx = find(shuffBin == bin); % get index for trials belonging to the current bin
                
                idx = idx(1:nPerBin*nCVBlocks); % drop excess trials
                
                x = repmat((1:nCVBlocks)',nPerBin,1); shuffBlocks(idx) = x; % assign randomly order trials to blocks - actually assign block IDs to data points
                
            end
            
            
            % unshuffle block assignment
            
            blocks(shuffInd) = shuffBlocks;
            
            % save block assignment
            
            blocksAssign(:,iter) = blocks; % block assignment
            
            nTrialsPerBlock = length(blocks(blocks == 1)); % # of trials per block
            
            % Average data for each position bin across blocks
            
            posBins = 1:nBins;
            
            blockDat_filtData = nan(nBins*nCVBlocks,nElectrodes,nSamps);  % averaged & filtered EEG data with resampling at 50 Hz
            
            labels = nan(nBins*nCVBlocks,1);  % bin labels for averaged & filtered EEG data
            
            blockNum = nan(nBins*nCVBlocks,1); % block numbers for averaged & filtered EEG data
            
            bCnt = 1;
            
            for ii = 1:nBins
                
                for iii = 1:nCVBlocks
                    
                    blockDat_filtData(bCnt,:,:) = squeeze(mean(filtData(posBin==posBins(ii) & blocks==iii,:,tois),1)); %downsample and average trials into blocks
                    labels(bCnt) = ii;
                    
                    blockNum(bCnt) = iii;
                    
                    bCnt = bCnt+1;
                    
                end
                
            end
            
            slidingWindowAverage = movmean(blockDat_filtData ,3, 3);
            % Do SVM_ECOC at each time point
            parfor t = 1:nSamps
                
                % grab data for timepoint t
                toi = ismember(downsampledTimes,downsampledTimes(t)-window/2:downsampledTimes(t)+window/2);
                % average across time window of interest
                dataAtTimeT = squeeze(mean(blockDat_filtData(:,:,toi),3));
                
                %actually, take a sliding window average across the tois (which are downsampled), and sample from
                %it
                dataAtTimeT = slidingWindowAverage(:,:,t);
                % Do SVM_ECOC for each block
                for i=1:nCVBlocks % loop through blocks, holding each out as the test set
                    
                    trnl = labels(blockNum~=i); % training labels
                    
                    tstl = labels(blockNum==i); % test labels
                    
                    trnD = dataAtTimeT(blockNum~=i,:);    % training data
                    
                    tstD = dataAtTimeT(blockNum==i,:);    % test data
                    
                    % SVM_ECOC
                    mdl = fitcecoc(trnD,trnl, ...
                        'Coding','ternarycomplete', ...
                        ... 'Coding','onevsall', ...
                        'Learners', svmTemplate);   %train support vector mahcine
                    LabelPredicted = predict(mdl, tstD);       % predict classes for new data
                    
                    svm_predict(iter,t,i,:) = LabelPredicted;  % save predicted labels
                    
                    tst_target(iter,t,i,:) = tstl;             % save true target labels
                    
                end % end of block
                
            end % end of time points
            
        end % end of iteration
        toc(iterTic)
        
        
        OutputfName = strcat(resultsDir, ...
            subjectName, '_', ...
            strjoin(conditions, '_'), ...
            '_results_central_electrodes', ...
            '.mat');
        
        svmECOC.data = data;
        svmECOC.data.electrodes = relevantElectrodes;
        svmECOC.data.nElectrodes = numel(relevantElectrodes);
        svmECOC.downsampledTimes = downsampledTimes;
        svmECOC.nCVBlocks = nCVBlocks;
        svmECOC.targets = tst_target;
        svmECOC.modelPredict = svm_predict;
        svmECOC.testResults =  svmECOC.targets == svmECOC.modelPredict;
        svmECOC.overallSuccessRatePcnt =  mean(svmECOC.testResults, "all") * 100;
        svmECOC.successRatesTime = mean(svmECOC.testResults, [1 3 4]);
        
        save(OutputfName,'svmECOC','-v7.3');
        toc(subjectTic)
        
    end % end of subject loop

    
    for i = 1:numel(subjects)
        sub = subjects(i);
        subjectName = num2str(sub, '%03.f')
        
        resultsFile = strcat(resultsDir, ...
            subjectName, '_', ...
            strjoin(conditions, '_'), ...
            '_results_central_electrodes',  ...
            '.mat');
        tmp =  load(resultsFile);
        svmECOC = tmp.svmECOC;
        results{i} = svmECOC.successRatesTime;
    end
    
    allResults = cat(1, results{:});
    times = svmECOC.downsampledTimes;
    
    numPlotsPerFigure = 6;
    plotIdx = repmat(1:numPlotsPerFigure, 1, floor(numel(subjects)/numPlotsPerFigure + 1));
    for i = 1:numel(subjects)
        if  mod(i, numPlotsPerFigure) == 1
            figure('units','normalized', 'WindowState', 'maximized')
            titleString = sprintf('Decoding conditions: %s', conditionDesc);
            % titleString = sprintf('Decoding conditions: %s %s %s',strrep(string(conditions), '_', '-'));
            sgtitle(titleString)
        end
        subplot(2, numPlotsPerFigure/2, plotIdx(i));
        plot(times, allResults(i, :)*100);
        hold on
        plot(times, repmat((1/numel(unique(labels)) * 100), 1, numel(times)), 'b--');
        plot(times, repmat((1/numel(unique(labels)) * 100 * 2), 1, numel(times)), 'b--');
        hold off
        
        ylim([0 120])
        xlim([times(1) times(end)]);
        xlabel('Time')
        ylabel('success rate %')
        titleString = sprintf('Sucess rate, subject %d', subjects(i));
        title(titleString)
        if  mod(i, numPlotsPerFigure) == 0 || i == numel(subjects)
            if  mod(i, numPlotsPerFigure) == 0
                firstSubIdx = i - numPlotsPerFigure +1;
            else %i == numel(subjects)
                firstSubIdx = i - mod(i, numPlotsPerFigure) +1
            end
            figureFileName = sprintf('%d-%d-%s',subjects(firstSubIdx), subjects(i), conditionDesc);
            figureFileName = strcat(figuresDir, 'latest\', figureFileName);
            print(gcf, figureFileName, '-djpeg',  '-r0');
        end
    end
    toc(analysisTic)
end