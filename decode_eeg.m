% This function does ERP-based decoding for orientation reported in Experiment 1.

% Structure of analysis and names of some variables used in this script

% were borrowed from the analysis script used in Foster et al. (2016).

% Note: Randomization routine included in this analysis (e.g., random permutation)

% can produce results slightly different from results reported in the paper. 

% Gi-Yeul Bae 2017.10.3.


function decode_eeg(C)
% delete(gcp)
% parpool
    if nargin ==  0
        C = Constants();
    end

    subjects = C.subjects;
    nSubjects = C.nSubjects;
    bApplyLowpassFilter = 0;
    bPlotResults = 0;
    
    combinedElectrodes = union(C.centralElectrodes, C.occipitalElectrodes);
    
    averagedSuccessRates = nan(nSubjects, 1);
    allSuccessRates = nan(C.nCVBlocks, nSubjects);

    relevantElectrodes = C.relevantElectrodes;
    nElectrodes = numel(relevantElectrodes);
    
    nCVBlocks = C.nCVBlocks; % otherwise the parfor loop won't work
    nClasses = C.nUniqueLables;
    
    analysisTic = tic;
    
    outputFileNameResults = strcat(C.resultsDir, ...
        C.conditionDesc, ...
        '_results', ...
        C.data_suffix, ...
        C.result_suffix, ...
        '.mat');

    
    %% Loop through participants
    for subjectIdx = 1:nSubjects
        subjectTic = tic;
        subject = subjects(subjectIdx);
        subjectName = num2str(subject, '%03.f');
        fprintf('Subject:\t%d\n',subject);

        outputFileNameDecoder = strcat(C.resultsDir, ...
                             'decoder-params-output/', ...
                             subjectName, '_', ...
                             C.conditionDesc, ...
                             '_decoder_params', ...
                             C.data_suffix, ...
                             C.result_suffix, ...
                             '.mat');


        inputFileName = strcat(C.resultsDir, ...
            subjectName, '_', ...
            'data_all_conditions', ...
            C.data_suffix, ...
            '.mat');

        load(inputFileName); % into "data" variable
                             % perform the actual translation of the labels
        if C.bTranslateLabels == 1
            for i = C.relevantConditionsIdx
                data.labels(data.labels == C.labels(i)) = C.labelTranslation(i);
            end
        end

        % filter out the irrelevant electrodes and trials based on stimulus (condition) label
        electrodeIdx = ismember(C.allElectrodesInOrder, relevantElectrodes);
        labelIdx = ismember(data.labels, C.labels);
        data.eeg = data.eeg(labelIdx, electrodeIdx, :);
        data.labels = data.labels(labelIdx);

        % create an array of time points belonging to the current analysis
        % - downsampled by a factor of 2
        downsampleFactor = 2;
        origSamplingInterval = fix(1000/C.origSamplingFreq);
        resamplingInterval = origSamplingInterval * downsampleFactor;
        downsampledTimes = data.time.pre:resamplingInterval:data.time.post;

        % create a boolean array with length equal to the number of timepoints in
        % the original data, where "1" at an index i denotes that timepoint i should
        % be used for this analysis - this will be used for the downsampling itself
        nOrigDataPoints = fix((data.time.post - data.time.pre)/1000 * C.origSamplingFreq);
        t1 = 1:nOrigDataPoints;
        t2 = 1:2:nOrigDataPoints
        relevantTimes = ismember(t1, t2);

        nSamps = length(downsampledTimes);
        nTimes = length(relevantTimes);
        
        % overall # of trials (all conditions)
        nTrials = numel(data.labels);

        % preallocate Matrices
        svm_predict = nan(C.nIter,nSamps,C.nCVBlocks,nClasses); % a matrix for saving prediction from SVM
        iterTestLabels = nan(C.nIter,nSamps,C.nCVBlocks,nClasses);  % a matrix for saving true target values
        blocksAssign = nan(nTrials,C.nIter);  % a matrix for saving block assignments

        filtData = data.eeg(:,:,1:nTimes);

        % bootstrapping iterations
        for iter = 1:C.nIter
            fprintf("training iteration %d\n", iter);
            iterTic = tic; % start timing iteration loop
            
            % block assignments for every trial
            blocks = nan(nTrials,1);
            
            shuffBlocks = nan(nTrials, 1);
            
            % count number of trials within each position bin
            clear binCnt
            for class = 1:nClasses
                calssCnt(class) = sum(data.labels == class);
            end
            
            minTrialsPerClass = min(calssCnt); % # of trials for position bin with fewest trials
            nPerClass = floor(minTrialsPerClass/C.nCVBlocks); % max # of trials such that the # of trials for each bin can be equated within each block
            
            % shuffle trials
            shuffInd = randperm(nTrials)'; % create shuffle index
            shuffTrials = data.labels(shuffInd); % shuffle trial order
            
            % take the 1st nPerClass x C.nCVBlocks trials for each position bin.
            for class = 1:nClasses;
                idx = find(shuffTrials == class); % get index for trials belonging to the current class
                idx = idx(1:nPerClass*C.nCVBlocks); % drop excess trials
                x = repmat((1:C.nCVBlocks)',nPerClass,1); shuffBlocks(idx) = x; % assign randomly order trials to blocks - actually assign block IDs to data points
            end
            
            
            % unshuffle block assignment
            blocks(shuffInd) = shuffBlocks;
            
            % save block assignment
            blocksAssign(:,iter) = blocks; % block assignment
            nTrialsPerBlock = length(blocks(blocks == 1)); % # of trials per block
            
            % Average data for each position bin across blocks
            blockData = nan(nClasses*C.nCVBlocks,nElectrodes,nSamps);  % averaged EEG data with resampling
            labels = nan(nClasses*C.nCVBlocks,1);  % bin labels for averaged EEG data
            blockNum = nan(nClasses*C.nCVBlocks,1); % block numbers for averaged EEG data
            
            blockCounter = 1;
            for class = 1:nClasses
                for block = 1:C.nCVBlocks
                    blockData(blockCounter,:,:) = squeeze(mean(filtData(data.labels==C.uniqueLables(class) & blocks == block,:,relevantTimes),1)); %downsample and average trials into blocks
                    labels(blockCounter) = class;
                    blockNum(blockCounter) = block;
                    blockCounter = blockCounter+1;
                end
            end
            
            % take a sliding window average across relevantTimes (which are downsampled)
            slidingWindowAverage = movmean(blockData ,3, 3);

           % train a classifier at each time point
           parfor t = 1:nSamps
               % sample from the sliding window average
                dataAtTimeT = slidingWindowAverage(:,:,t);

                % iterate over cross validation blocks, holding each out as the test set
                for i=1:nCVBlocks
                    trnl = labels(blockNum~=i); % training labels
                    tstl = labels(blockNum==i); % test labels
                    trnD = dataAtTimeT(blockNum~=i,:);    % training data
                    tstD = dataAtTimeT(blockNum==i,:);    % test data
                    
                    % train the SVM
                    mdl = fitcecoc(trnD,trnl, ...
                        'Coding','ternarycomplete', ...
                        ... 'Coding','onevsall', ...
                        'Learners', C.svmTemplate);   %train support vector mahcine

                    LabelPredicted = predict(mdl, tstD); % predict classes for new data
                    
                    svm_predict(iter,t,i,:) = LabelPredicted;  % save predicted labels
                    
                    iterTestLabels(iter,t,i,:) = tstl; % save true target labels
                end % end of block
            end % end of time points
        end % end of iteration

        toc(iterTic)
        
        % save decoding parametersa and results for later analysis
        decoder.data = data;
        decoder.data.electrodes = relevantElectrodes;
        decoder.data.nElectrodes = numel(relevantElectrodes);
        decoder.downsampledTimes = downsampledTimes;
        decoder.nCVBlocks = C.nCVBlocks;
        decoder.targets = iterTestLabels;
        decoder.modelPredict = svm_predict;
        decoder.testResults =  decoder.targets == decoder.modelPredict;
        decoder.overallSuccessRatePcnt =  mean(decoder.testResults, "all") * 100;
        decoder.successRatesTime = mean(decoder.testResults, [1 3 4]);
        decoder.params = C;
        
        results{subjectIdx} = decoder.successRatesTime;
        
        save(outputFileNameDecoder,'decoder','-v7.3');
        toc(subjectTic)
        
    end % end of subject loop

    
    allResults = 100 * cat(1, results{:});
    decodingResults.successRates = allResults;
    decodingResults.downsampledTimes = decoder.downsampledTimes;
    decodingResults.nClasses = nClasses;
    save(outputFileNameResults,'decodingResults','-v7.3');

    if bPlotResults == 1
        plot_results(C);
    end

    toc(analysisTic)

end
