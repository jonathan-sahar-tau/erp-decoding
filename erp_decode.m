% Based on work by Gi-Yeul Bae 2017.10.3.

function erp_decode(subjects)
    % delete(gcp);
    % parpool;
    % eeglab;

    if nargin == 0
        subjects = [11 13 14 15 16]
    end

% parameters to set

    nClasses = 2; % # of possible stimuli we want to differentiate between
    nIter = 1; % # of iterations
    nAveragedTrials = 10; % # of averagedTrials for cross-validation
    lowPassBorders = [0 6]; % low pass filter
    allTimePoints = -100:2:900; % time points of interest in the analysis
    windowSizeMs = 4; % 1 data point per 4 ms in the preprocessed data
    samplingFreq = 512; % samplring rate of the preprocessed data for filtering
    electrodes = 1:64;
    nElectrodes = length(electrodes); % # of electrode included in the analysis
    nSamples = 0;

    nCVBlocks = 3;
    allSuccessRates = nan(nCVBlocks, numel(subjects));
    conditions = ["cong_int", "cong_scr"];
    nConditions = size(conditions, 2);
    fieldsToDelete = {'filename', 'filepath', 'subject', 'group', ...
                    'condition', 'setname', 'session', 'srate', ...
                    'xmin', 'xmax', 'icaact', 'icawinv', ...
                    'icasphere', 'icaweights', 'icachansind', 'chanlocs', ...
                    'urchanlocs', 'chaninfo', 'ref', 'event', ...
                    'urevent', 'eventdescription', 'epoch', 'epochdescription', ...
                    'reject', 'stats', 'specdata', 'specicaact', ...
                    'splinefile', 'icasplinefile', 'dipfit', 'history', ...
                    'saved', 'etc', 'run'};
    %% Loop through participants
    for subjectIdx = 1:numel(subjects)
        subject = subjects(subjectIdx);
        subjectName = num2str(subject, '%03.f');
        fprintf('Subject:\t%d\n',subject);
        dataLocation = strcat("/home/jonathan/google_drive/Msc neuroscience/lab", ...
                              "/data/experiment_data"); % set directory of data set

        outputDir = strcat(dataLocation, '/', 'output');
        % conditions = ["cong_int", "inc_int", "cong_scr", "inc_scr"]
        minTrialnum = intmax;

        % assign an integer value to each condition, so we can later put all the data in
        % one matrix and reference it with this number.
        enumVal = 1;

        for condition = conditions
            currentFilename= strcat(subjectName, '_', condition, '_bc', '.vhdr');
            fprintf("fetching data from file: %s\n", strcat(dataLocation, '/', currentFilename))
            tempData = (pop_loadbv(dataLocation, currentFilename, [], 1:64));
            % remove uneeded fields from the struct
            eegDataAndMetadata = rmfield(tempData, fieldsToDelete);

            % swap the 2nd and 3rd dimensions, effectively transposing the
            % "inner" (electrode) matrices to be nTrials x nSamples -
            % every row is the time course of one trial
            eegDataAndMetadata.data = permute(eegDataAndMetadata.data, [1,3,2]);
            subjectEEGs.(condition) = eegDataAndMetadata;
            subjectEEGs.(condition).name = condition;
            subjectEEGs.(condition).enumVal = enumVal;
            enumVal = enumVal + 1;

            % get the minimal number of trials acrosso conditions
            if (subjectEEGs.(condition).trials < minTrialnum) minTrialnum = subjectEEGs.(condition).trials;, end
            % set directory for decoding results.

            subjectEEGs.(condition).nTrialsPerAveragedTrial = floor(minTrialnum/nAveragedTrials);

            % the actrual # of trials that can accomodate our desired # of averagedTrials,
            % such that all averagedTrials contain the same # of trials.
            subjectEEGs.(condition).nTrials = subjectEEGs.(condition).nTrialsPerAveragedTrial * nAveragedTrials;
            nTrials = subjectEEGs.(condition).nTrials;

            % create a mapping for shuffling the trials: the value of index[i]
            % gives the index of trial that was shuffled into the i-th
            % position for this iteration
            subjectEEGs.(condition).averagedTrialIndex = repmat((1:nAveragedTrials)',subjectEEGs.(condition).nTrialsPerAveragedTrial, 1);
        end % conditions

        % we assume here that all conditions are normalized to have the same
        % number of time points in each trial recording
        nSamples = subjectEEGs.(conditions(1)).pnts;

        EEGs = nan(nElectrodes, nTrials, nSamples, nConditions);
        for condition = conditions % low-pass filtering
            conditionEnum = subjectEEGs.(condition).enumVal;
            subjectEEGs.(condition).filteredData = ...
                                    nan(nElectrodes, ...
                                        nTrials, ...
                                        subjectEEGs.(condition).pnts);

            for electrode = 1:nElectrodes
            % TODO: fixme
            subjectEEGs.(condition).filteredData(electrode,:,:) = ...
                subjectEEGs.(condition).data(electrode,1:nTrials,:);

            EEGs(electrode, :, :, conditionEnum) = ...
                subjectEEGs.(condition).data(electrode,1:nTrials,:);
            % subjectEEGs.(condition).filteredData(electrode,:,:) = ...
            % eegfilt(...
            %     squeeze(subjectEEGs.(condition).data(electrode,1:nTrials,:)),...
            %     samplingFreq,...
            %     lowPassBorders(1, 1),...
            %     lowPassBorders(1, 2));
            end
        end % Loop through each iteration

        testLabels = nan(nSamples * nConditions, nAveragedTrials, nIter);
        testPredictions = nan(nSamples * nConditions, nAveragedTrials, nIter);
        % for iter = 1:nIter
        iter = 1
            for condition = conditions
                % iter = 1

                % we're shuffling the indices of the trials - we'll then use
                % these to  randomly assign a averagedTrial to each trial.
                subjectEEGs.(condition).trialIndex = randperm(nTrials)';
                % unshuffle averagedTrial assignment - remember that averagedTrialIndex is
                % just a cyclyc list of assignments: e.g. 123123123123
                averagedTrialsAssignment(subjectEEGs.(condition).trialIndex) = ...
                                    subjectEEGs.(condition).averagedTrialIndex;

                % averaged & filtered eeg data
                averagedTrials = nan(nElectrodes,nSamples, nAveragedTrials);
                filteredData = subjectEEGs.(condition).filteredData;

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

                subjectEEGs.(condition).averagedTrials = averagedTrials;
            end %condition

                % subjectEEGs.(conditions(1)).averagedTrials(:,:,:) = -1 ;
                % subjectEEGs.(conditions(2)).averagedTrials(:,:,:) = 1;
            averagedTrialIds = 1:nAveragedTrials;
            percentCorrect = nan(1,nCVBlocks);
            indices = randperm(nAveragedTrials);
            trainingStart = tic; % start timing iteration loop
            for cvBlock = indices(1:nCVBlocks) % cross validation across averagedTrials
                blockStart = tic; % start timing iteration loop
                trainLabels = nan(nSamples * (nAveragedTrials - 1), nConditions);
                iterTestLabels = nan(nSamples, nConditions);
                trainData = nan(nElectrodes, nSamples * (nAveragedTrials - 1), nConditions);
                testData = nan(nElectrodes, nSamples, nConditions);
                % collect data for this division of CV - DON'T PARALLEL!
                for condition = conditions
                    conditionData = subjectEEGs.(condition);
                    conditionAveragedTrials = conditionData.averagedTrials;
                    conditionEnum = conditionData.enumVal;

                    % fprintf('Collecting data, averagedTrial: %d, condition: %s\n', ...
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

            testLabels(:, cvBlock, iter) = reshape(iterTestLabels,  1, [])';

            trainData = reshape(trainData, size(trainData, 1), []);

            testData = reshape(testData, size(testData, 1), []);

            fprintf('Training SVM, holding out averagedTrial %d for testing\n',cvBlock)
            SVMModel = fitcsvm(trainData', trainLabels, 'ShrinkagePeriod', 1000);
            % SVM_ECOC model = fitcecoc(trainData,trainLabels,
            % 'Coding','onevsall','Learners','SVM' ); %train support
            % vector mahcine
            classLoss = kfoldLoss(crossval(SVMModel))
            %store test results for this <averagedTrial,iter>
            testPredictions(:, cvBlock, iter) = predict(SVMModel, testData');
            percentCorrect(1, cvBlock) =  mean(squeeze(testPredictions(:, cvBlock, iter) == testLabels(:, cvBlock, iter))) * 100;
            fprintf("success rate, block %d: %.1f%%\n", cvBlock, percentCorrect(cvBlock))
            elapsed = toc(blockStart);
            elapsed = datevec(elapsed./(60*60*24));
            fprintf("time elapsed for this CV block: %d minutes, %.0f seconds\n", elapsed(end-1), elapsed(end))
            end % end of cross validation

        % end % end of iteration

        totalElapsed = toc(trainingStart); % stop timing the iteration loop
        totalElapsed = datevec(totalElapsed./(60*60*24));
        percentCorrect
        allSuccessRates(:, subjectIdx) =  percentCorrect;
        averagedPercentCorrect(subjectIdx) =  mean(percentCorrect);
        fprintf("overall success rate, iteration %d: %.1f%%\n", iter, averagedPercentCorrect(subjectIdx))
        fprintf("total time elapsed: %d minutes, %.0f seconds\n", totalElapsed(end-1), totalElapsed(end))
        % allResults(subject) = percentCorrect;
        % OutputfName = strcat(outputDir,'/Orientation_Results_ERPbased_', subjectName,'.mat');
        % fprintf("percent success overall: %d%%\n" ,mean(allResults))
        % save(OutputfName,'decoder','-v7.3');
    end % subjects
end
