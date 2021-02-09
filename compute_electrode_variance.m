% This function does ERP-based decoding for orientation reported in Experiment 1.

% Structure of analysis and names of some variables used in this script

% were borrowed from the analysis script used in Foster et al. (2016).

% Note: Randomization routine included in this analysis (e.g., random permutation)

% can produce results slightly different from results reported in the paper. 

% Gi-Yeul Bae 2017.10.3.


function compute_electrode_variance(subjects)
% delete(gcp)
% parpool

    C = Constants();
    subjects = C.subjects;
    nSubjects = C.nSubjects;
    bApplyLowpassFilter = 0;
    bPlotResults = 1;
    
    combinedElectrodes = union(C.centralElectrodes, C.occipitalElectrodes);
    relevantElectrodes = C.allElectrodes;
    nElectrodes = numel(relevantElectrodes); % # of electrode included in the analysis
    
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

        outputfName = strcat(C.resultsDir, ...
                             C.conditionDesc, ...
                             '_electrode_variance', ...
                             C.data_suffix, ...
                             C.result_suffix, ...
                             '.mat');

    allSubjectStds = nan(nSubjects, C.nConditions, nElectrodes);
    allSubjectData = nan(nSubjects, C.nConditions, nElectrodes);
    %% Loop through participants
    
    for subjectIdx = 1:nSubjects
        subjectTic = tic;
        subject = subjects(subjectIdx);
        subjectName = num2str(subject, '%03.f');
        fprintf('Subject:\t%d\n',subject);



        inputFileName = strcat(C.resultsDir, ...
            subjectName, '_', ...
            'data_all_conditions', ...
            C.data_suffix, ...
            '.mat');

        load(inputFileName); % into "data" variable

        % filter out the irrelevant electrodes and trials based on simulus (condition) label
        electrodeIdx = ismember(C.allElectrodesInOrder, relevantElectrodes);
        labelIdx = ismember(data.labels, relevantLabels);
        data.eeg = data.eeg(labelIdx, electrodeIdx, :);
        data.labels = data.labels(labelIdx);

        allCondStds = nan(C.nConditions, nElectrodes);
        allCondData = nan(C.nConditions, nElectrodes);

        for condIdx = 1:C.nConditions
            labelIdx = (data.labels == relevantLabels(condIdx));
            conditionData = squeeze(data.eeg(labelIdx, :, :));
            conditionStd = squeeze(std(conditionData, 1)); %across trials
            conditionMeanStd = squeeze(mean(conditionStd, 2)); %across time points, should be 1 x nElectrodes
            allCondStds(condIdx, :) = conditionMeanStd;
            allCondData(condIdx,:, :) = squeeze(mean(mean(conditionData, 1), 3)); % average data across trials, then time
        end


       allSubjectStds(subjectIdx, :, :) = allCondStds;
       allSubjectData(subjectIdx, :, :) = allCondData;

    end % end of subject loop

    allConditions_averageStd = squeeze(mean(allSubjectStds, 1)); % across conditions
    SEs = squeeze(std(allConditions_averageStd)/sqrt(size(allConditions_averageStd, 1))); % standard error of the across-conditions std
    allSubjects_averageStd = squeeze(mean(allConditions_averageStd, 1)); % across subjects    
    averagedData = mean(squeeze(mean(allSubjectData, 2)), 1);
    
    variance.allSubjectStds = allSubjectStds;
    variance.averagedData = averagedData;
    variance.subjectIds = subjects;
    
    save(outputfName,'variance','-v7.3');
    
    %  plot results 
    % ----------------------
    if bPlotResults == 1
        figure('units','normalized', 'WindowState', 'maximized')
        titleString = sprintf("Mean standard deviation across subjects and conditions, and SE thereof");
        title(titleString,'Interpreter','none')

        hold on
        xticks(1:nElectrodes);
        ylim([0,20]);
        set(gca,'XTickLabel', relevantElectrodes);
        xtickangle(45)
        bar(allSubjects_averageStd, "faceColor", rand(1,3));
        plot(1:nElectrodes, repmat(mean(allSubjects_averageStd),1, nElectrodes));
        er = errorbar(1:62, allSubjects_averageStd, SEs);
        er.Color = [0 0 0];
        er.LineStyle = 'none';
       
        bar(abs(averagedData),  "faceColor", rand(1,3));
        
        figureFileName = sprintf('%s-%s-%s-%s','all subjects', 'all-conditions', 'electrode_SD', C.data_suffix);
        figureFileName = strcat(C.figuresDir, 'latest\', figureFileName);
        print(gcf, figureFileName, '-djpeg',  '-r0');
    end
end
