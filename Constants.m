classdef Constants
    properties %( Constant = true )
        
        subjects;
        nSubjects;

        origLabels
        translateLabels
        labelTranslation
        nUniqueLables
        conditionDescriptors
        conditions
        nConditions

        baseDir
        dataLocation
        outputDir
        resultsDir
        figuresDir
        eeglabPath
        
        allElectrodesInOrder = ["Vertical", "Horizontal", "Fp1", "AF7", "AF3", "F1", "F3", "F5", "F7", "FT7", "FC5", "FC3", "FC1", "C1", "C3", "C5", "T7", "TP7", "CP5", "CP3", "CP1", "P1", "P3", "P5", "P7", "P9", "PO7", "PO3", "O1", "Iz", "Oz", "POz", "Pz", "CPz", "Fpz", "Fp2", "AF8", "AF4", "AFz", "Fz", "F2", "F4", "F6", "F8", "FT8", "FC6", "FC4", "FC2", "FCz", "Cz", "C2", "C4", "C6", "T8", "TP8", "CP6", "CP4", "CP2", "P2", "P4", "P6", "P8", "P10", "PO8"]

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

        nIter = 1; % # of iterations
        frequencies = [0 30]; % low pass filter
        nCVBlocks = 3;
        window = 4; % 1 data point per 4 ms in the preprocessed data
        Fs = 512 ; % sampling rate of the preprocessed data for filtering
        svmTemplate = templateSVM('KernelFunction','linear');
        nAveragedTrials = 15;
        
        fieldsToDelete = {'filename', 'filepath', 'subject', 'group', ...
            'condition', 'setname', 'session', 'srate', ...
            'xmin', 'xmax', 'icaact', 'icawinv', ...
            'icasphere', 'icaweights', 'icachansind', 'chanlocs', ...
            'urchanlocs', 'chaninfo', 'ref', 'event', ...
            'urevent', 'eventdescription', 'epoch', 'epochdescription', ...
            'reject', 'stats', 'specdata', 'specicaact', ...
            'splinefile', 'icasplinefile', 'dipfit', 'history', ...
            'saved', 'etc', 'run'};


    end
       methods
           function obj = Constants(subjects)
%                    obj.subjects = [102]

                   obj.subjects = [102 104:106 108:112 114:116 118:120 122]
                                      
                % obj.baseDir("G:\My Drive\MudrikLab020818\Experiments_new\Jonathan\erp-decoding\")
                % obj.dataLocation = strcat(baseDir, "experiment_data\");
                % obj.dataLocation = strcat(baseDir, "experiment_data\unconscious_conditions\");

                obj.nSubjects = numel(obj.subjects);
                obj.baseDir = "C:\Users\Jonathan\Google Drive\Msc neuroscience\lab\analysis_scripts\erp-decoding\";
                obj.dataLocation = strcat(obj.baseDir, "experiment_data\conscious_conditions_long_exposure\");
                obj.outputDir = strcat(obj.baseDir, "output\");
                obj.resultsDir = strcat(obj.outputDir, 'mat-files\');
                obj.figuresDir = strcat(obj.outputDir, 'figures\');

                obj.eeglabPath = strcat(obj.baseDir, "software\eeglab2020_0");
                
                obj.leftFrontalElectrodes = ["Fp1", "AF7", "AF3", "F3", "F5", "F7"];
                obj.middleFrontalElectrodes = ["FpZ", "AFZ", "F1" "Fz", "F2"];
                obj.rightFrontalElectrodes = ["Fp2", "AF8", "AF4", "F4", "F6", "F8"];
                
                obj.leftCentralElectrodes = ["FT7", "FC5", "FC3", "C3", "C5", "T7", "TP7", "CP5", "CP3"];
                obj.middleCentralElectrodes = ["FC1", "FCz", "FC2", "C1", "C2", "Cz", "CP2", "CP1", "CPZ"];
                obj.rightCentralElectrodes = ["FC4", "FC6", "FT8", "C4", "C6", "T8", "CP4", "CP6", "TP8"];
                
                obj.leftOccipitalElectrodes = ["P3", "P5", "P7", "P9", "PO3", "PO7", "O1"];
                obj.middleOcccipitalElectrodes = ["Pz", "P1", "P2", "POz", "Oz", "Iz"];
                obj.rightOccipitalElectrodes = ["P4", "P6", "P8", "P10", "PO8"];
                % obj.rightOccipitalElectrodes = ["P4", "P6", "P8", "P10", "PO4", "PO8", "O2"];
                
                obj.frontalElectrodes = union(union( obj.leftFrontalElectrodes,  obj.middleOcccipitalElectrodes),  obj.rightFrontalElectrodes);
                obj.centralElectrodes = union(union( obj.leftCentralElectrodes,  obj.middleOcccipitalElectrodes),  obj.rightCentralElectrodes);
                obj.occipitalElectrodes = union(union( obj.leftOccipitalElectrodes,  obj.middleOcccipitalElectrodes),  obj.rightOccipitalElectrodes);

                obj.conditions = ["ConInt", "IncInt", "ConScr", "IncScr"]
                obj.origLabels =        [1, 2, 3, 4];
                obj.labelTranslation =  [1, 2, 3, 4];
                obj.nConditions = numel(obj.conditions);

                obj.translateLabels = 0; % whether to join together different conditions, e.g. "ConInt" and "IncInt" and tread them as "AllInc"

                if obj.translateLabels == 1
                    obj.nUniqueLables = numel(unique(obj.labelTranslation)); % # of different classes

                else
                    obj.nUniqueLables = numel(unique(obj.origLabels)); % # of different classes
                end

%                 obj.conditionDescriptors = nan(numel(obj.conditions));
                for i = 1:numel(obj.conditions)
                    conditionDescriptor.name = obj.conditions(i)
                    conditionDescriptor.label = obj.origLabels(i)
                    conditionDescriptor.labelTranslation = obj.labelTranslation(i)
                    obj.conditionDescriptors{i} = conditionDescriptor;
                end



           end
       end
end
