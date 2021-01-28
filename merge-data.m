function svmECOC_our_data(subjects)

    C = Constants();
    subjects = C.subjects;
    nSubjects = C.nSubjects;
    debugAddBump = 0
    debugApplyLowpassFilter = 1
    % parameters to set

    combinedElectrodes = union(C.centralElectrodes, C.occipitalElectrodes);


    conditions = ["ConInt"] %, "IncInt", "ConScr", "IncScr"]
    nConditions = numel(conditions);

    origLabels = [1, 2, 3, 4]; % must be same length as conditions array


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
            filePath = strcat(C.dataLocation, '\', fileName);

            tempData = load(filePath);
            tempData = tempData.EEGData;
            % tempData = rmfield(tempData, C.fieldsToDelete);

            % swap the the matrix dimensions so that they're
            % nTrials x nElectrodes x nSamples
            tempData.data = permute(tempData.data, [3, 1, 2]);
            EEG.(condition) = tempData;
            EEG.(condition).name = condition;

            % get the minimal number of trials across conditions
            trialNums(conditionIdx) = EEG.(condition).trials;
        end

        % find the actual # of trials that can accomodate our desired # of averagedTrials,
        nTrialsPerAveragedTrial = floor(min(trialNums)/C.nAveragedTrials);
        nTrialsPerCondition = nTrialsPerAveragedTrial * C.nAveragedTrials;

        % normalize all conditions to have the same number of time points in each trial recording
        for conditionIdx  = 1:numel(conditions)
            condition = conditions(conditionIdx);
            electrodeIdx = ismember(C.allElectrodesInOrder, relevantElectrodes);
            conditionData{conditionIdx} =  ...
                EEG.(condition).data(1:nTrialsPerCondition,...
                :, :);
            conditionLabels{conditionIdx} = repmat(origLabels(conditionIdx), nTrialsPerCondition, 1);
        end

        allData = cat(1, conditionData{:});
        allLabels = cat(1, conditionLabels{:});

        data.eeg = allData;
        data.labels = allLabels;
        data.times = floor(EEG.(conditions(1)).times);
        data.time.pre = data.times(1);
        data.time.post = data.times(end);


        outputfName = strcat(C.resultsDir, ...
            subjectName, '_', ...
            'data_all_conditions', ...
            suffix, ...
            '.mat');

        save(outputfName,'data','-v7.3');
    end
end
