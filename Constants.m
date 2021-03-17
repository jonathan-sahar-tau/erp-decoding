classdef Constants
    properties
        
        subjects;
        nSubjects;

        origConditions
        origLabels
        origLabelTranslation

        conditions
        labels
        labelTranslation

        nConditions
        relevantConditionsIdx
        bTranslateLabels
        uniqueLables
        nUniqueLables

        baseDir
        dataLocation
        outputDir
        resultsDir
        figuresDir
        eeglabPath

        % Suffixes
        % -----------
        initialDataSuffix = '_BC'
        data_suffix = '_filtered_30_Hz'
        result_suffix = '_allElectrodes'
        conditionDesc = "Congruency"


        allElectrodesInOrder = ["Vertical", "Horizontal", "Fp1", "AF7", "AF3", "F1", "F3", "F5", "F7", "FT7", "FC5", "FC3", "FC1", "C1", "C3", "C5", "T7", "TP7", "CP5", "CP3", "CP1", "P1", "P3", "P5", "P7", "P9", "PO7", "PO3", "O1", "Iz", "Oz", "POz", "Pz", "CPz", "Fpz", "Fp2", "AF8", "AF4", "AFz", "Fz", "F2", "F4", "F6", "F8", "FT8", "FC6", "FC4", "FC2", "FCz", "Cz", "C2", "C4", "C6", "T8", "TP8", "CP6", "CP4", "CP2", "P2", "P4", "P6", "P8", "P10", "PO8"]

        allElectrodes
        relevantElectrodes
        
        customElectrodeGroup
        customElectrodeGroupName
        
        leftFrontalElectrodes
        middleFrontalElectrodes
        rightFrontalElectrodes
        
        leftCentralElectrodes
        middleCentralElectrodes
        rightCentralElectrodes
        
        leftOccipitalElectrodes
        middleOcccipitalElectrodes
        rightOccipitalElectrodes
            
        frontalElectrodes
        centralElectrodes
        occipitalElectrodes

        

        
        nIter = 10; % # of iterations
        lowpassFrequencies = [0 30]; % low pass filter borders
        nCVBlocks = 3; % #of cross validation blocks
        origSamplingFreq = 512 ; % sampling rate of the preprocessed data

        svmTemplate = templateSVM('KernelFunction','linear');


    end
       methods
           function obj = Constants(subjects)

                % Data params
                % -----------

                obj.subjects = [11:12, 15:17, 19:27, 29, 30, 35:38, 41, 42];
                obj.origConditions = ["Cong_Int", "Inc_Int", "Cong_Scr", "Inc_Scr"];

                obj.origLabels =        [1, 2, 3, 4];
                obj.origLabelTranslation =  [1, 1, 2, 2];

                obj.nSubjects = numel(obj.subjects);

                % whether or note to aggregate different conditions, e.g. "ConInt" and "IncInt"
                % and treat them as "AllInc"
                obj.bTranslateLabels = 0;

                % choose which conditions to include in the analysis.
                % by default - use all conditions.
                obj.relevantConditionsIdx = 1:numel(obj.origConditions);

                % Paths
                % -----------
                obj.baseDir = "C:\Users\Jonathan\Google Drive\Msc neuroscience\lab\analysis_scripts\erp-decoding\";
                obj.dataLocation = strcat(obj.baseDir, "experiment_data\exp_itai_unconscious\");
                obj.outputDir = strcat(obj.baseDir, "output\");
                obj.resultsDir = strcat(obj.outputDir, 'mat-files\');
                obj.figuresDir = strcat(obj.outputDir, 'figures\');

                %                 obj.resultsDir = strcat(obj.outputDir, 'mat-files-temp\');

                obj.eeglabPath = strcat(obj.baseDir, "software\eeglab2020_0");

                obj = setLabelProperties(obj);

                % Electrodes params
                % -----------------

                % remove horizontal and vertical "electrodes"
                obj.allElectrodes = obj.allElectrodesInOrder(3:end)

                obj.leftFrontalElectrodes = ["Fp1", "AF7", "AF3", "F3", "F5", "F7"];
                obj.middleFrontalElectrodes = ["FpZ", "AFZ", "F1" "Fz", "F2"];
                obj.rightFrontalElectrodes = ["Fp2", "AF8", "AF4", "F4", "F6", "F8"];
                
                obj.leftCentralElectrodes = ["FT7", "FC5", "FC3", "C3", "C5", "T7", "TP7", "CP5", "CP3"];
                obj.middleCentralElectrodes = ["FC1", "FCz", "FC2", "C1", "C2", "Cz", "CP2", "CP1", "CPZ"];
                obj.rightCentralElectrodes = ["FC4", "FC6", "FT8", "C4", "C6", "T8", "CP4", "CP6", "TP8"];
                
                obj.leftOccipitalElectrodes = ["P3", "P5", "P7", "P9", "PO3", "PO7", "O1"];
                obj.middleOcccipitalElectrodes = ["Pz", "P1", "P2", "POz", "Oz", "Iz"];
                obj.rightOccipitalElectrodes = ["P4", "P6", "P8", "P10", "PO8"];
                
                obj.frontalElectrodes = union(union( obj.leftFrontalElectrodes,  obj.middleOcccipitalElectrodes),  obj.rightFrontalElectrodes);
                obj.centralElectrodes = union(union( obj.leftCentralElectrodes,  obj.middleOcccipitalElectrodes),  obj.rightCentralElectrodes);
                obj.occipitalElectrodes = union(union( obj.leftOccipitalElectrodes,  obj.middleOcccipitalElectrodes),  obj.rightOccipitalElectrodes);

                obj.relevantElectrodes = obj.allElectrodes;
            end %constructor


            function obj = setLabelProperties(obj)
                % narrow downs the condition & labels if required
                obj.conditions = obj.origConditions(obj.relevantConditionsIdx);
                obj.labels = obj.origLabels(obj.relevantConditionsIdx);
                obj.labelTranslation = obj.origLabelTranslation(obj.relevantConditionsIdx);
                obj.nConditions = numel(obj.conditions);

                if obj.bTranslateLabels == 1
                    obj.uniqueLables = unique(obj.labelTranslation);
                else
                    obj.uniqueLables = obj.labels; % # of different classes
                end
                    obj.nUniqueLables = numel(obj.uniqueLables); % # of different classes

            end %setLabelProperties

       end
end
