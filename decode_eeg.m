% This function does ERP-based decoding for orientation reported in Experiment 1.

% Structure of analysis and names of some variables used in this script

% were borrowed from the analysis script used in Foster et al. (2016).

% Note: Randomization routine included in this analysis (e.g., random permutation)

% can produce results slightly different from results reported in the paper. 

% Gi-Yeul Bae 2017.10.3.


function decode_eeg(subjects)
% delete(gcp)
% parpool

    C = Constants();
    subjects = C.subjects;
    nSubjects = C.nSubjects;
    bApplyLowpassFilter = 0
    bPlotResults = 0
    
    combinedElectrodes = union(C.centralElectrodes, C.occipitalElectrodes);
    
    averagedSuccessRates = nan(nSubjects, 1);
    allSuccessRates = nan(C.nCVBlocks, nSubjects);

    relevantElectrodes = C.allElectrodesInOrder(3:end);
    nElectrodes = numel(relevantElectrodes); % # of electrode included in the analysis
    
    nCVBlocks    = C.nCVBlocks; % otherwise the parfor loop won't work
    % global values that we may want to override
    % ------------------------------
    relevantConditionsIdx = 1:C.nConditions;
    relevantLabels = C.origLabels;
    nClasses = C.nUniqueLables;
    
    % overriding global values here
    % ------------------------------
%     relevantConditionsIdx = [1, 2] %, 3, 4];
%     relevantLabels = C.origLabels(relevantConditionsIdx);
%     nClasses = numel(relevantLabels);
%     
    analysisTic = tic;
    
    %% Loop through participants
    for subjectIdx = 1:nSubjects
        subjectTic = tic;
        subject = subjects(subjectIdx);
        subjectName = num2str(subject, '%03.f');
        fprintf('Subject:\t%d\n',subject);

        outputfName = strcat(C.resultsDir, ...
                             subjectName, '_', ...
                             C.conditionDesc, ...
                             '_results', ...
                             C.data_suffix, ...
                             C.result_suffix, ...
                             '.mat');


        inputFileName = strcat(C.resultsDir, ...
            subjectName, '_', ...
            'data_all_conditions', ...
            C.data_suffix, ...
            '.mat');

        load(inputFileName); % into "data" variable

        if C.translateLabels == 1
            for i = relevantConditionsIdx
                data.labels(data.labels == C.conditionDescriptors{i}.labelTranslation) = C.conditionDescriptors{i}.labelTranslation;
            end
        end

        % filter out the irrelevant electrodes and trials based on simulus (condition) label
        electrodeIdx = ismember(C.allElectrodesInOrder, relevantElectrodes);
        labelIdx = ismember(data.labels, relevantLabels);
        data.eeg = data.eeg(labelIdx, electrodeIdx, :);
        data.labels = data.labels(labelIdx);

        resamplingRatio = 1000/250; % resample the data during analysis to 250 Hz
        downsampledTimes = data.time.pre:resamplingRatio:data.time.post; % time points of interest in the analysis
        nSamps = length(downsampledTimes);
        
        % set up time points
        tois = ismember(data.time.pre:2:data.time.post,downsampledTimes);
        nTimes = length(tois);
        
        % # of trials
        nTrials = numel(data.labels);

        % Preallocate Matrices
        svm_predict = nan(C.nIter,nSamps,C.nCVBlocks,nClasses); % a matrix to save prediction from SVM
        iterTestLabels = nan(C.nIter,nSamps,C.nCVBlocks,nClasses);  % a matrix to save true target values
        
        blocksAssign = nan(nTrials,C.nIter);  % create block to save block assignments
        
        % low-pass filtering
        filtData = nan(nTrials,nElectrodes,nTimes);

        if bApplyLowpassFilter == 1
            tic
            for c = 1:nElectrodes
                    tmp = eegfilt(squeeze(data.eeg(:,c,:)),C.Fs,C.frequencies(1,1),C.frequencies(1,2)); % low pass filter
                    filtData(:,c,:) = tmp(:,1:nTimes);
            end
            toc
        else
            filtData = data.eeg(:,:,1:nTimes);
        end

        
        % Loop through each iteration
        for iter = 1:C.nIter
            fprintf("training iteration %d\n", iter);
            iterTic = tic; % start timing iteration loop
            
            % preallocate arrays
            blocks = nan(nTrials,1);
            
            shuffBlocks = nan(nTrials, 1);
            
            % count number of trials within each position bin
            
            clear binCnt
            
            for bin = 1:nClasses
                
                binCnt(bin) = sum(data.labels == bin);
                
            end
            
            minCnt = min(binCnt); % # of trials for position bin with fewest trials
            
            nPerBin = floor(minCnt/C.nCVBlocks); % max # of trials such that the # of trials for each bin can be equated within each block
            
            % shuffle trials
            
            shuffInd = randperm(nTrials)'; % create shuffle index
            
            shuffBin = data.labels(shuffInd); % shuffle trial order
            
            % take the 1st nPerBin x C.nCVBlocks trials for each position bin.
            
            for bin = 1:nClasses;
                
                idx = find(shuffBin == bin); % get index for trials belonging to the current bin
                
                idx = idx(1:nPerBin*C.nCVBlocks); % drop excess trials
                
                x = repmat((1:C.nCVBlocks)',nPerBin,1); shuffBlocks(idx) = x; % assign randomly order trials to blocks - actually assign block IDs to data points
                
            end
            
            
            % unshuffle block assignment
            
            blocks(shuffInd) = shuffBlocks;
            
            % save block assignment
            
            blocksAssign(:,iter) = blocks; % block assignment
            
            nTrialsPerBlock = length(blocks(blocks == 1)); % # of trials per block
            
            % Average data for each position bin across blocks

            blockDat_filtData = nan(nClasses*C.nCVBlocks,nElectrodes,nSamps);  % averaged & filtered EEG data with resampling at 50 Hz
            
            labels = nan(nClasses*C.nCVBlocks,1);  % bin labels for averaged & filtered EEG data
            
            blockNum = nan(nClasses*C.nCVBlocks,1); % block numbers for averaged & filtered EEG data
            
            bCnt = 1;
            
            for ii = 1:nClasses
                
                for iii = 1:C.nCVBlocks
                    
                    blockDat_filtData(bCnt,:,:) = squeeze(mean(filtData(data.labels==C.uniqueLables(ii) & blocks==iii,:,tois),1)); %downsample and average trials into blocks
                    labels(bCnt) = ii;
                    
                    blockNum(bCnt) = iii;
                    
                    bCnt = bCnt+1;
                    
                end
                
            end
            
            slidingWindowAverage = movmean(blockDat_filtData ,3, 3);
            % Do SVM_ECOC at each time point
           parfor t = 1:nSamps
                
                % grab data for timepoint t
                toi = ismember(downsampledTimes,downsampledTimes(t)-C.window/2:downsampledTimes(t)+C.window/2);
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
                        'Learners', C.svmTemplate);   %train support vector mahcine
                    LabelPredicted = predict(mdl, tstD);       % predict classes for new data
                    
                    svm_predict(iter,t,i,:) = LabelPredicted;  % save predicted labels
                    
                    iterTestLabels(iter,t,i,:) = tstl;             % save true target labels
                    
                end % end of block
                
            end % end of time points
            
        end % end of iteration
        toc(iterTic)
        
        
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
        
        save(outputfName,'decoder','-v7.3');
        toc(subjectTic)
        
    end % end of subject loop

    if bPlotResults == 1
    % load saved results files
    % ----------------------
    for i = 1:nSubjects
        sub = subjects(i);
        subjectName = num2str(sub, '%03.f')
        
        resultsFile = strcat(C.resultsDir, ...
            subjectName, '_', ...
            strjoin(C.conditions, '_'), ...
            '_results',  ...
            C.results_suffix,  ...
            '.mat');
        load(resultsFile); % into "decoder"
        results{i} = decoder.successRatesTime;
    end
    
    allResults = cat(1, results{:});
    times = decoder.downsampledTimes;


    % plot results
    % ------------
    numPlotsPerFigure = 6;
    plotIdx = repmat(1:numPlotsPerFigure, 1, floor(nSubjects/numPlotsPerFigure + 1));
    for i = 1:nSubjects
        if  mod(i, numPlotsPerFigure) == 1
            figure('units','normalized', 'WindowState', 'maximized')
            titleString = sprintf('Decoding conditions: %s%s', C.conditionDesc, C.results_suffix);
            % titleString = sprintf('Decoding conditions: %s %s %s',strrep(string(conditions), '_', '-'));
            sgtitle(titleString)
        end
        subplot(2, numPlotsPerFigure/2, plotIdx(i));
        plot(times, allResults(i, :)*100);
        hold on
        plot(times, repmat(((1/nClasses) * 100), 1, numel(times)), 'b--');
        plot(times, repmat((1/nClasses * 100 * 2), 1, numel(times)), 'b--');
        hold off
        
        ylim([0 120])
        xlim([times(1) times(end)]);
        xlabel('Time')
        ylabel('success rate %')
        titleString = sprintf('Sucess rate, subject %d ', subjects(i));
        title(titleString)
        if  mod(i, numPlotsPerFigure) == 0 || i == nSubjects
            if  mod(i, numPlotsPerFigure) == 0
                firstSubIdx = i - numPlotsPerFigure +1;
            else %i == nSubjects
                firstSubIdx = i - mod(i, numPlotsPerFigure) +1
            end
            figureFileName = sprintf('%d-%d-%s-%s',subjects(firstSubIdx), subjects(i), C.conditionDesc, C.results_suffix);
            figureFileName = strcat(C.figuresDir, 'latest\', figureFileName);
            print(gcf, figureFileName, '-djpeg',  '-r0');
        end
    end
    toc(analysisTic)
    end
end
