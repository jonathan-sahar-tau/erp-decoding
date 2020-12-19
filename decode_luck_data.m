% Based on work by Gi-Yeul Bae 2017.10.3.

function decode_luck_data(subjects)
    delete(gcp('nocreate'))
    parpool;
    eeglab;

    if nargin == 0
        subjects = [505, 506, 507, 508, 509, 510];
    end
    subjects = [505]; % 14 15 16];
        subjects
% parameters to set

    nIter = 1; % # of iterations
    nAveragedTrials = 5; % # of averagedTrials for cross-validation
    % nAveragedTrials = 10; % # of averagedTrials for cross-validation
    lowPassBorders = [0 30]; % low pass filter
    windowSizeMs = 4; % 1 data point per 4 ms in the preprocessed data
    samplingFreq = 512; % samplring rate of the preprocessed data for filtering
    electrodes = 1:27;
    nElectrodes = length(electrodes); % # of electrode included in the analysis
    nSamples = 0;

    dataLocation = strcat("/home/jonathan/google_drive/Msc neuroscience/lab", ...
                            "/data/data_luck_2018/Exp1_Data"); % set directory of data set
    outputDir = strcat(dataLocation, '/'); %, 'output');

    outputFile = strcat(outputDir, '/', datestr(now,'mm-dd-yyyy-HH-MM'), '_run-log');
    fprintf("output file: %s", outputFile);
    fileId = fopen(outputFile, 'w');

    nCVBlocks = 3;
    averagedSuccessRates = nan(numel(subjects), 1);
    allSuccessRates = nan(nCVBlocks, numel(subjects));
    %% Loop through participants
    for subjectIdx = 1:numel(subjects)
        subject = subjects(subjectIdx);
        subjectName = num2str(subject, '%03.f');
        fprintf(fileId, 'Subject:\t%d\n',subject);
        % conditions = ["cong_int", "inc_int", "cong_scr", "inc_scr"]
        minTrialnum = intmax;

        % assign an integer value to each condition, so we can later put all the data in
        % one matrix and reference it with this number.
        enumVal = 1;

        currentFilename= strcat("Decoding_DE_", subjectName, '.mat');
        fprintf("fetching data from file: %s\n", strcat(dataLocation, '/', currentFilename));
        fprintf(fileId, "fetching data from file: %s\n", strcat(dataLocation, '/', currentFilename));

        %load the curren subject's data into "data" variable
        load(strcat(dataLocation, "/", currentFilename));
        conditions = unique(data.channel);
        nConditions = size(conditions, 2);

        fmt = ['Running decoder, subject no. %d,  conditions: [', repmat('%s, ', 1, numel(conditions)-1), '%s]\n'];
        fprintf(fileId, fmt, subject, conditions);
        fprintf(fileId, "===================================================\n");
        % swap the 2nd and 3rd dimensions, effectively transposing the
        % "inner" (electrode) matrices to be nTrials x nSamples -
        % every row is the time course of one trial
        data.eeg = permute(data.eeg, [2,1,3]);
        subjectEEGs = data
        for condition = conditions

            subjectEEGs.allConditions{condition}.name = condition;
            subjectEEGs.allConditions{condition}.enumVal = condition;
            enumVal = enumVal + 1;
            subjectEEGs.allConditions{condition}.data = data.eeg(:, data.channel == condition, :);
            subjectEEGs.allConditions{condition}.trials = size(subjectEEGs.allConditions{condition}.data, 2) ;
            % get the minimal number of trials acrosso conditions
            if (subjectEEGs.allConditions{condition}.trials < minTrialnum) minTrialnum = subjectEEGs.allConditions{condition}.trials;, end
            % set directory for decoding results.
        end

        for condition = conditions
            subjectEEGs.allConditions{condition}.nTrialsPerAveragedTrial = floor(minTrialnum/nAveragedTrials);

            fprintf(fileId, "nTrialsPerAveragedTrial, condition %s: %d trials\n", condition, subjectEEGs.allConditions{condition}.nTrialsPerAveragedTrial);
            % the actrual # of trials that can accomodate our desired # of averagedTrials,
            % such that all averagedTrials contain the same # of trials.
            origNTrials = subjectEEGs.allConditions{condition}.trials;
            subjectEEGs.allConditions{condition}.nTrials = ...
                subjectEEGs.allConditions{condition}.nTrialsPerAveragedTrial * nAveragedTrials;
            nTrials = subjectEEGs.allConditions{condition}.nTrials; %global variable
            fprintf(fileId, "decreasing nTrials, condition %s: from %d to %d trials\n", condition, origNTrials, nTrials);

            % create a mapping for shuffling the trials: the value of index[i]
            % gives the index of trial that was shuffled into the i-th
            % position for this iteration
            subjectEEGs.allConditions{condition}.averagedTrialIndex = repmat((1:nAveragedTrials)',subjectEEGs.allConditions{condition}.nTrialsPerAveragedTrial, 1);
        end % conditions

        % we assume here that all conditions are normalized to have the same
        % number of time points in each trial recording
        nSamples = size(subjectEEGs.allConditions{conditions(1)}.data, 3);

        for condition = conditions % low-pass filtering
            conditionEnum = subjectEEGs.allConditions{condition}.enumVal;
            subjectEEGs.allConditions{condition}.filteredData = ...
                                    nan(nElectrodes, ...
                                        nTrials, ...
                                        nSamples);
            filteredTemp = nan(size(subjectEEGs.allConditions{condition}.filteredData));
            parfor electrode = 1:nElectrodes
            % subjectEEGs.allConditions{condition}.filteredData(electrode,:,:) = ...
            %     subjectEEGs.allConditions{condition}.data(electrode,1:nTrials,:);

            data = squeeze(subjectEEGs.allConditions{condition}.data(electrode,1:nTrials,:));
            % EEGs(electrode, :, :, conditionEnum) = ...
            %     subjectEEGs.allConditions{condition}.data(electrode,1:nTrials,:);
            fprintf("filtering...")
            filteredTemp(electrode,:,:) = ...
            eegfilt(...
                data,...
                samplingFreq,...
                lowPassBorders(1, 1),...
                lowPassBorders(1, 2));
            end
            subjectEEGs.allConditions{condition}.filteredData = filteredTemp;
        end % Loop through each iteration

        testLabels = nan(nSamples * nConditions, nCVBlocks, nIter);
        testPredictions = nan(nSamples * nConditions, nCVBlocks , nIter);
        % for iter = 1:nIter
        iter = 1;
            for condition = conditions
                averagedTrialsAssignment = nan(nTrials, 1);
                % we're shuffling the indices of the trials - we'll then use
                % these to randomly designate an averagedTrial ID for each trial.
                subjectEEGs.allConditions{condition}.trialIndex = randperm(nTrials)';
                % unshuffle averagedTrial assignment - remember that averagedTrialIndex is
                % just a cyclyc list of assignments: e.g. 123123123123
                averagedTrialsAssignment(subjectEEGs.allConditions{condition}.trialIndex) = ...
                                    subjectEEGs.allConditions{condition}.averagedTrialIndex;

                % averaged & filtered eeg data
                averagedTrials = nan(nElectrodes,nSamples, nAveragedTrials);
                filteredData = subjectEEGs.allConditions{condition}.filteredData;

                for averagedTrial = 1:nAveragedTrials
                    % averagedTrials == averagedTrialAssignment means: get me all the trial data that
                    % has been assigned to averagedTrial "averagedTrial". this results in a
                    % matrix of dimensions nElectrodes x nSamples per each averagedTrial
                    averagedTrials(:,:, averagedTrial) = ...
                    squeeze(mean(filteredData(:, averagedTrialsAssignment == averagedTrial, :), 2));

                    % smooth out the ERP signal with a moving average
                    averagedTrials(:, :, averagedTrial) = ...
                    movmean(squeeze(averagedTrials(:, :, averagedTrial)), 4, ...
                    2);
                end % averagedTrials

                subjectEEGs.allConditions{condition}.averagedTrials = averagedTrials;
            end %condition

                % subjectEEGs.(conditions(1)).averagedTrials(:,:,:) = veragedTrialsAssignment()-1 ;
                % subjectEEGs.(conditions(2)).averallllgedTrials(:,:,:) = 1;
            averagedTrialIds = 1:nAveragedTrials;
            cvBlocks = randperm(nAveragedTrials);
            cvSuccessRates = nan(1,nCVBlocks);
            cvBlocks = cvBlocks(1:nCVBlocks)
            trainingStart = tic; % start timing iteration loop
            parfor cvBlockId = 1:nCVBlocks % cross validation across averagedTrials
                cvBlock = cvBlocks(cvBlockId);
                blockStart = tic; % start timing iteration loop
                % labels for this iteration's train set - all but one averagedTrial
                trainLabels = nan(nSamples * (nAveragedTrials - 1), nConditions);
                % labels for this iteration's test set
                iterTestLabels = nan(nSamples, nConditions);
                trainData = nan(nElectrodes, nSamples * (nAveragedTrials - 1), nConditions);
                testData = nan(nElectrodes, nSamples, nConditions);
                size(blockStart)
                size(trainLabels)
                size(iterTestLabels)
                size(trainData)
                size(testData)
                % collect data for this division of CV - DON'T PARALLEL!
                for condition = conditions
                    conditionData = subjectEEGs.allConditions{condition};
                    conditionAveragedTrials = conditionData.averagedTrials;
                    conditionEnum = conditionData.enumVal;

                    % fprintf(fileId, 'Collecting data, averagedTrial: %d, condition: %s\n', ...
                    %         cvBlock, condition)

                    trainLabels(:, conditionEnum) = repmat(conditionEnum,...
                                        size(trainLabels, 1), 1);

                    % tile the averagedTrials into a 2-D matrix of the data for this condition and
                    % this CV round
                    trainData(:, :,conditionEnum) = ...
                        reshape(conditionAveragedTrials(:, :, averagedTrialIds ~= cvBlock), ...
                                                                    nElectrodes, []);

                    iterTestLabels(: , conditionEnum) = repmat(conditionEnum,...
                                        size(iterTestLabels, 1), 1);

                    testData(:, :,conditionEnum) = ...
                                squeeze(conditionAveragedTrials (:, :, averagedTrialIds == cvBlock));
                    end

            trainLabels = reshape(trainLabels, 1, [])';

            testLabels(:, cvBlockId, iter) = reshape(iterTestLabels,  1, [])';

            trainData = reshape(trainData, size(trainData, 1), []);

            testData = reshape(testData, size(testData, 1), []);

            fprintf('Training SVM, subject %d, holding out CV block %d for testing\n', subject, cvBlock)
            % fprintf(fileId, 'Training SVM, holding out averagedTrial %d for testing\n',cvBlock);
            % SVMModel = fitcsvm(trainData', trainLabels, 'ShrinkagePeriod', 1000, "BoxConstraint", 111, "KernelScale", 13);
            t = templateSVM('KernelFunction','gaussian', ...
                            'KernelScale', 'auto', ...
                            'ShrinkagePeriod',1000, ...
                            'CacheSize', 14000, ...
                            'Verbose', 1)
            t = templateSVM('Verbose', 1)
            SVMModel = fitcecoc(trainData', trainLabels, ...
                                'Coding', 'onevsall', ...
                                'Learners', t);

            fprintf('Done training SVM');
            % SVM_ECOC model = fitcecoc(trainData,trainLabels,
            % 'Coding','onevsall','Learners','SVM' ); %train support
            % vector mahcine
            classLoss = kfoldLoss(crossval(SVMModel))
            %store test results for this <averagedTrial held out, iter>
            testPredictions(:, cvBlockId, iter) = predict(SVMModel, testData');
            cvSuccessRates(cvBlockId) =  mean(squeeze(testPredictions(:, cvBlockId, iter) == testLabels(:, cvBlockId, iter))) * 100;
            % fprintf(fileId, "success rate, block %d: %.1f%%\n", cvBlock, cvSuccessRates(cvBlockId));
            elapsed = toc(blockStart);
            elapsed = datevec(elapsed./(60*60*24))
            % fprintf(fileId, "time elapsed for this CV block: %d minutes, %.0f seconds\n", elapsed(end-1), elapsed(end));
            end % end of cross validation

        % end % end of iteration

        totalElapsed = toc(trainingStart); % stop timing the iteration loop
        totalElapsed = datevec(totalElapsed./(60*60*24));
        cvSuccessRates
        allSuccessRates(:, subjectIdx) =  cvSuccessRates;
        averagedSuccessRates(subjectIdx) =  mean(cvSuccessRates);

        % fileID = fopen('output.org','a');
        fprintf("overall success rate, subject %d: %.1f%%\n", subject, averagedSuccessRates(subjectIdx))
        fprintf("total time elapsed for subject: %d minutes, %.0f seconds\n", totalElapsed(end-1), totalElapsed(end));

        fprintf(fileId, "overall success rate, subject %d: %.1f%%\n", subject, averagedSuccessRates(subjectIdx));
        fprintf(fileId, "total time elapsed for subject: %d minutes, %.0f seconds\n", totalElapsed(end-1), totalElapsed(end));
        % OutputfName = strcat(outputDir,'/Orientation_Results_ERPbased_', subjectName,'.mat');
        % save(OutputfName,'decoder','-v7.3');;
    end % subjects

    fprintf("overall average success rate: %.1f%%\n" ,mean(averagedSuccessRates));
    fprintf(fileId, "overall average success rate: %.1f%%\n" ,mean(averagedSuccessRates));
    allSuccessRates
    averagedSuccessRates
    fclose(fileId); % close the log file descriptor
    save(strcat(outputDir, "/", subjects.mat), 'subjects');
    save(strcat(outputDir, "/", allSuccessRates.mat), 'allSuccessRates');
end
