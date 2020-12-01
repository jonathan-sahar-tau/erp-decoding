% Based on work by Gi-Yeul Bae 2017.10.3.

function erp_decode(subjects)
    % delete(gcp);
    % parpool;
    % eeglab;

    if nargin == 0 subjects = [11]; end

% parameters to set
    subjects = [13 14 15 16];
    subjects = [13];

    nClasses = 2; % # of possible stimuli we want to differentiate between
    nIter = 1; % # of iterations
    nBlocks = 10; % # of blocks for cross-validation
    lowPassBorders = [0 6]; % low pass filter
    allTimePoints = -100:2:900; % time points of interest in the analysis
    windowSizeMs = 4; % 1 data point per 4 ms in the preprocessed data
    samplingFreq = 512; % samplring rate of the preprocessed data for filtering
    electrodes = 1:64;
    nElectrodes = length(electrodes) % # of electrode included in the analysis
    nSamples = 0;

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
    for subject = subjects

        subjectName = num2str(subject, '%03.f');
        fprintf('Subject:\t%d\n',subject);

        % load data subjectName = num2str(subject,'%03.f');
        dataLocation = strcat("/home/jonathan/google_drive/Msc neuroscience/lab", ...
                              "/data/exported_data"); % set directory of data set

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

            subjectEEGs.(condition).nTrialsPerBlock = floor(minTrialnum/nBlocks);

            % the actrual # of trials that can accomodate our desired # of blocks,
            % such that all blocks contain the same # of trials.
            subjectEEGs.(condition).nTrials = subjectEEGs.(condition).nTrialsPerBlock * nBlocks;
            nTrials = subjectEEGs.(condition).nTrials;

            % create a mapping for shuffling the trials: the value of index[i]
            % gives the index of trial that was shuffled into the i-th
            % position for this iteration
            subjectEEGs.(condition).blockIndex = repmat((1:nBlocks)',subjectEEGs.(condition).nTrialsPerBlock, 1);
        end % conditions

        % we assume here that all conditions are normalized to have the same
        % number of time points in each trial recording
        nSamples = subjectEEGs.(conditions(1)).pnts;

        lowPassStart = tic;
        for condition = conditions % low-pass filtering
            subjectEEGs.(condition).filteredData = ...
                                    nan(nElectrodes, ...
                                        nTrials, ...
                                        subjectEEGs.(condition).pnts);

            for electrode = 1:nElectrodes
            % TODO: fixme
            subjectEEGs.(condition).filteredData(electrode,:,:) = ...
                subjectEEGs.(condition).data(electrode,1:nTrials,:);

            % subjectEEGs.(condition).filteredData(electrode,:,:) = ...
            % eegfilt(...
            %     squeeze(subjectEEGs.(condition).data(electrode,1:nTrials,:)),...
            %     samplingFreq,...
            %     lowPassBorders(1, 1),...
            %     lowPassBorders(1, 2));
            end
        end % Loop through each iteration
        toc(lowPassStart)
        % end
        % lowPassElapsed = toc(lowPassStart)
        % fprintf('elapsed time: %d', lowPassElapsed)

        testLabels = nan(nSamples * nConditions, nBlocks, nIter);
        testPredictions = nan(nSamples * nConditions, nBlocks, nIter);
        for iter = 1:nIter
            iterStart = tic; % start timing iteration loop
            for condition = conditions
                % iter = 1

                % we're shuffling the indices of the trials - we'll then use
                % these to  randomly assign a block to each trial.
                subjectEEGs.(condition).trialIndex = randperm(nTrials)';
                % unshuffle block assignment - remember that blockIndex is
                % just a cyclyc list of assignments: e.g. 123123123123
                blocks(subjectEEGs.(condition).trialIndex) = ...
                                    subjectEEGs.(condition).blockIndex;

                % averaged & filtered eeg data
                averagedERPBlocks = nan(nElectrodes,nSamples, nBlocks);
                filteredData = subjectEEGs.(condition).filteredData;

                for block = 1:nBlocks
                    % blocks == block means: get me all the trial data that
                    % has been assigned to block "block". this results in a
                    % matrix of dimensions nElectrodes x nSamples per each block
                    averagedERPBlocks(:,:, block) = ...
                    squeeze(mean(filteredData(:, blocks == block, :), 2));

                    % smooth out the ERP signal with a moving average
                    averagedERPBlocks(:, :, block) = ...
                    movmean(squeeze(averagedERPBlocks(:, :, block)), 4, ...
                    2);
                end % blocks

                subjectEEGs.(condition).averagedERPBlocks = averagedERPBlocks;
            end %condition

            blockIds = 1:nBlocks;
            for cvBlock = 1:nBlocks % cross validation across blocks
                trainLabels = nan(nSamples * (nBlocks - 1), nConditions);
                iterTestLabels = nan(nSamples, nConditions);
                trainData = nan(nElectrodes, nSamples * (nBlocks - 1), nConditions);
                testData = nan(nElectrodes, nSamples, nConditions);
                % collect data for this division of CV - DON'T PARALLEL!
                for condition = conditions
                    conditionData = subjectEEGs.(condition);
                    conditionBlocks = conditionData.averagedERPBlocks;
                    conditionEnum = conditionData.enumVal;

                    fprintf('Collecting data, block: %d, condition: %s\n', ...
                            cvBlock, condition)

                    trainLabels(:, conditionEnum) = repmat(conditionEnum,...
                                        size(trainLabels, 1), 1);

                    % tile the blocks into a 2-D matrix of the data for this condition and
                    % this CV round
                    trainData(:, :,conditionEnum) = ...
                        reshape(conditionBlocks(:, :, blockIds ~= cvBlock), ...
                                                                    nElectrodes, []);

                    iterTestLabels(: , conditionEnum) = repmat(conditionEnum,...
                                        size(iterTestLabels, 1), 1);

                    testData(:, :,conditionEnum) = ...
                                squeeze(conditionBlocks (:, :, blockIds == cvBlock));
                    end

            trainLabels = reshape(trainLabels, 1, []);

            testLabels(:, cvBlock, iter) = reshape(iterTestLabels,  1, [])';

            trainData = reshape(trainData, size(trainData, 1), []);

            testData = reshape(testData, size(testData, 1), []);

            fprintf('Training SVM, block: %d\n',cvBlock)

                    model = fitcsvm(trainData', trainLabels');
                    % SVM_ECOC model = fitcecoc(trainData,trainLabels,
                    % 'Coding','onevsall','Learners','SVM' ); %train support
                    % vector mahcine

                    %store test results for this <block,iter>
                    testPredictions(:, cvBlock, iter) = predict(model, testData');
            end % end of block

        end % end of iteration


    toc(iterStart) % stop timing the iteration loop
        percentCorrect =  mean2(reshape(testPredictions, [], 1) == reshape(testLabels, [], 1)) * 100;
        fprintf("percent success, iteration %d: %d%%\n", iter, percentCorrect)
        allResults(subject) = percentCorrect;
        % OutputfName = strcat(outputDir,'/Orientation_Results_ERPbased_', subjectName,'.mat');
        fprintf("percent success overall: %d%%\n" ,mean(allResults))
        % save(OutputfName,'decoder','-v7.3');

    end % subjects


end

