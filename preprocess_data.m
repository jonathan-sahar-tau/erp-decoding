function preprocess_data(subjects)
    % This function goes over .mat files (possibly imported by import_data.m),
    % that contain the trial data from different conditions of the same subject, and creates
    % a single .mat file which contains all trials from all conditions, as well as a vector
    % of labels (describing the condition) that match the trials

    % In addition, it can also apply a low pass filter to the data.

    C = Constants();
    bApplyLowpassFilter = 1
    
    %% Loop through participants
    for subjectIdx = 1:C.nSubjects
        subjectTic = tic;
        subject = C.subjects(subjectIdx);
        subjectName = num2str(subject, '%03.f');
        fprintf('Subject:\t%d\n',subject);


        for conditionIdx  = 1:C.nConditions
            condition = C.conditions(conditionIdx)
            fileName = strcat(subjectName, '_', condition, '.mat')
            filePath = strcat(C.dataLocation, '\', fileName);

            load(filePath); %into EEGData
            % swap the the matrix dimensions so that they're
            % nTrials x nElectrodes x nSamples
            EEGData.data = permute(EEGData.data, [3, 1, 2]);
            EEG.(condition) = EEGData;
            EEG.(condition).name = condition;

            % get the minimal number of trials across C.conditions
            trialNums(conditionIdx) = EEG.(condition).trials;
        end

        % find the actual # of trials that can accomodate our desired # of averagedTrials,
        nTrialsPerCondition = min(trialNums);

        % normalize all C.conditions to have the same number of time points in each trial recording
        for conditionIdx  = 1:C.nConditions
            condition = C.conditions(conditionIdx);
            conditionData{conditionIdx} =  ...
                EEG.(condition).data(1:nTrialsPerCondition,...
                :, :);
            conditionLabels{conditionIdx} = repmat(C.origLabels(conditionIdx), nTrialsPerCondition, 1);
        end

        allData = cat(1, conditionData{:});
        allLabels = cat(1, conditionLabels{:});
        
             
        data.times = floor(EEG.(C.conditions(1)).times(1:size(allData, 3)));
        data.time.pre = data.times(1);
        data.time.post = data.times(end);
        

        % create a boolean array with length equal to the number of timepoints in
        % the original data, where "1" at an index i denotes that timepoint i should
        % be used for this analysis - this will be used for the downsampling itself
        % then create an array of time points belonging to the current analysis
        % - downsampled by a factor of 2
        downsampleFactor = 2;

        nOrigDataPoints = numel(data.times);
        t1 = 1:nOrigDataPoints;
        t2 = 1:downsampleFactor:nOrigDataPoints;
        relevantPoints = ismember(t1, t2);
        nTimes = length(relevantPoints);
        downsampledTimes = data.times(relevantPoints);

        nTotalTrials = numel(allLabels);

        nElectrodes = numel(C.allElectrodesInOrder);
        filteredData = nan(nTotalTrials,nElectrodes,nTimes);
       if bApplyLowpassFilter == 1
            parfor electrodeIdx = 1:nElectrodes
                    tmp = eegfilt(squeeze(allData(:,electrodeIdx,:)),C.origSamplingFreq,C.lowpassFrequencies(1,1),C.lowpassFrequencies(1,2)); % low pass filter
                    filteredData(:,electrodeIdx,:) = tmp(:,1:nTimes);
            end
            allData = filteredData;
            data.downsampledTimes = downsampledTimes;
       else
       
       end

        
        data.eeg = allData;
        data.labels = allLabels;


        outputfName = strcat(C.resultsDir, ...
            subjectName, '_', ...
            'data_all_conditions', ...
            C.data_suffix, ...
            '.mat');

        save(outputfName,'data','-v7.3');
    end
end
