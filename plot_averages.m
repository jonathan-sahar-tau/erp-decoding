function plot_averages(subjects)
   

    if nargin == 0
                subjects = [102 104:106 108:112 114:116 118:120 122]
%         subjects = [114:116 118:120 122]
    end
    
    nSubs = length(subjects);

    % parameters to set
    conditions = ["ConInt", "ConScr", "IncInt", "IncScr"];
%     conditions = ["Cong_Int", "Cong_Scr", "Inc_Int"];
    nTimes = 512;
    % allElectrodesInOrder = ["Fp1", "AF7", "AF3", "F1", "F3", "F5", "F7", "FT7", "FC5", "FC3", "FC1", "C1", "C3", "C5", "T7", "TP7", "CP5", "CP3", "CP1", "P1", "P3", "P5", "P7", "P9", "PO7", "PO3", "O1", "Iz", "Oz", "POz", "Pz", "CPz", "Fpz", "Fp2", "AF8", "AF4", "AFz", "Fz", "F2", "F4", "F6", "F8", "FT8", "FC6", "FC4", "FC2", "FCz", "Cz", "C2", "C4", "C6", "T8", "TP8", "CP6", "CP4", "CP2", "P2", "P4", "P6", "P8", "P10", "PO8", "PO4", "O2"]

    allElectrodesInOrder = ["Vertical", "Horizontal", "Fp1", "AF7", "AF3", "F1", "F3", "F5", "F7", "FT7", "FC5", "FC3", "FC1", "C1", "C3", "C5", "T7", "TP7", "CP5", "CP3", "CP1", "P1", "P3", "P5", "P7", "P9", "PO7", "PO3", "O1", "Iz", "Oz", "POz", "Pz", "CPz", "Fpz", "Fp2", "AF8", "AF4", "AFz", "Fz", "F2", "F4", "F6", "F8", "FT8", "FC6", "FC4", "FC2", "FCz", "Cz", "C2", "C4", "C6", "T8", "TP8", "CP6", "CP4", "CP2", "P2", "P4", "P6", "P8", "P10", "PO8"]

    leftFrontalElectrodes = ["Fp1", "AF7", "AF3", "F3", "F5", "F7"];
    middleFrontalElectrodes = ["FpZ", "AFZ", "F1" "Fz", "F2"];
    rightFrontalElectrodes = ["Fp2", "AF8", "AF4", "F4", "F6", "F8"];

    leftCentralElectrodes = ["FT7", "FC5", "FC3", "C3", "C5", "T7", "TP7", "CP5", "CP3"];
    middleCentralElectrodes = ["FC1", "FCz", "FC2", "C1", "C2", "Cz", "CP2", "CP1", "CPZ"];
    rightCentralElectrodes = ["FC4", "FC6", "FT8", "C4", "C6", "T8", "CP4", "CP6", "TP8"];
    
    leftOccipitalElectrodes = ["P3", "P5", "P7", "P9", "PO3", "PO7", "O1"];
    middleOcccipitalElectrodes = ["Pz", "P1", "P2", "POz", "Oz", "Iz"];
    rightOccipitalElectrodes = ["P4", "P6", "P8", "P10", "PO4", "PO8", "O2"];

    frontalElectrodes = union(union(leftFrontalElectrodes, middleOcccipitalElectrodes), rightFrontalElectrodes);
    centralElectrodes = union(union(leftCentralElectrodes, middleOcccipitalElectrodes), rightCentralElectrodes);
    occipitalElectrodes = union(union(leftOccipitalElectrodes, middleOcccipitalElectrodes), rightOccipitalElectrodes);



    regions{1} = leftFrontalElectrodes
    regions{2} = middleFrontalElectrodes
    regions{3} = rightFrontalElectrodes
    regions{4} = leftCentralElectrodes
    regions{5} = middleCentralElectrodes
    regions{6} = rightCentralElectrodes
    regions{7} = leftOccipitalElectrodes
    regions{8} = middleOcccipitalElectrodes
    regions{9} = rightOccipitalElectrodes


    regionNames = ["Left Frontal", "Middle Frontal", "Right Frontal", "Left Central", "Middle Central", "Right Central", "Left Occipital", "Middle Occcipital", "Right Occipital"]

    % baseDir("G:\My Drive\MudrikLab020818\Experiments_new\Jonathan\")
    baseDir = "C:\Users\Jonathan\Google Drive\Msc neuroscience\lab\";
    eeglabPath = strcat(baseDir, "software\eeglab2020_0")
    addpath(eeglabPath);

    dataLocation = strcat(baseDir, "analysis_scripts\erp-decoding\experiment_data\conscious_conditions_long_exposure\");
    outputDir = strcat(baseDir, "analysis_scripts\erp-decoding\output\");    
    OutputfName = strcat(outputDir, ...
        'grand_average_waveforms', ...
        '.mat');
    

% nElectrodes x nTimes x nConditions x nSubjects
    averages = nan(numel(allElectrodesInOrder), nTimes, numel(conditions), numel(subjects));
    
    reloadData = 1

    if reloadData == 1
        %% Loop through participants
        for subjectIdx = 1:numel(subjects)
            subject = subjects(subjectIdx);
            subjectName = num2str(subject, '%03.f');
            fprintf('Subject:\t%d\n',subject);
            
            % assign an integer value to each condition, so we can later put all the data in
            % one matrix and reference it with this number.
            
            for conditionIdx = 1:numel(conditions)
                condition = conditions(conditionIdx);
                currentFilename = strcat(subjectName, '_', condition, '_bc', '.mat');
                filePath = strcat(dataLocation, '\', currentFilename);
                t = load(filePath);
                EEGData = t.EEGData;
                EEGData.comments
                data = EEGData.data;
                avg =  squeeze(mean(data, 3));
                averages(:,:,conditionIdx, subjectIdx) = avg; %average across trials
                
                % times = EEGData.times;
                %                 for electrodeIdx = 14
                %                     for trialNum = 1:3
                %                         figure
                %                         plot(times, zeros(numel(times)))
                %                         hold on
                %
                %                         titleString = sprintf('waveform - subject %s - condition %s - electrode %s - trial %d', subjectName, condition, allElectrodesInOrder(electrodeIdx), trialNum);
                %                         sgtitle(titleString)
                %                         plot(EEGData.times, data(electrodeIdx, :,trialNum));
                %                         hold off
                %                     end
                %                 end
                
                %
                %
                %                 for electrodeIdx = 3:5 %:numel(allElectrodesInOrder)
                %                     figure
                %                 plot(times, zeros(numel(times)))
                %                 hold on
                
                %                     titleString = sprintf('Average waveforms - subject %s - condition %s - electrode %s', subjectName, condition, allElectrodesInOrder(electrodeIdx));
                %                     sgtitle(titleString)
                %                     plot(EEGData.times, averages(electrodeIdx,:,conditionIdx, subjectIdx))
                %                 hold off
                %                 end

            end
        end
        
        
        grandAverage = mean(averages, 4, 'omitnan');
        times = EEGData.times;
        
        rawData.grandAverage = grandAverage;
        rawData.times = times;
        save(OutputfName,'rawData','-v7.3');
        
    else
        load(OutputfName);
        grandAverage = rawData.grandAverage;
        times = rawData.times;
    end
    
    figure('units','normalized', 'WindowState', 'maximized')
    titleString = sprintf('Average waveforms');
    sgtitle(titleString)
    for r = 1:numel(regions)
        electrodeIdx = ismember(allElectrodesInOrder, regions{r});
        subplot(3,3,r);
        plot(times, zeros(numel(times)))
        hold on
        for c = 1:numel(conditions)
            timeIdx = 1:numel(times);
            dataToPlot = mean(grandAverage(electrodeIdx, timeIdx ,c), 1, 'omitnan');
            plot(times, dataToPlot)
        end
        hold off
        ylim([-15 10])
        xlim([times(1) times(end)]);
        ylabel('Potential (uV)')
        xlabel('Time (ms)')
        titleString = sprintf('%s', regionNames(r));
        title(titleString)
    end
    
    legend(conditions)
end
