function import_data(subjects)
    % This function iterates over data that was exported from brainVision's Analyzer software,
    % and converts them into .m files using EEGLab's load function.
    C = Constants();
    subjects = C.subjects;
    nSubs = length(subjects);
    baseDir = C.baseDir;
    eeglabPath = C.eeglabPath;
    dataLocation = C.dataLocation;
    outputDir = C.outputDir;
    numElectrodes = 64;
    exportDataLocation = C.dataLocation;
    addpath(eeglabPath);

%% Loop through participants
    for subjectIdx = 1:numel(subjects)
        subject = subjects(subjectIdx);
        subjectName = num2str(subject, '%03.f');
        fprintf('\n\nSubject:\t%d\n\n',subject);
        
        for condition = C.conditions
            fprintf('\n\nCondition:\t%s\n\n',condition);
            currentFilename= strcat(subjectName, '_', condition, C.initialDataSuffix, '.vhdr');
            exportFileName = strcat(subjectName, '_', condition, '.mat');
            exportLocation = strcat(exportDataLocation, '\', exportFileName);
            
            fprintf('\n\nimporting file:\t%s\n\n',currentFilename);

            
            EEGData = (pop_loadbv(dataLocation, currentFilename, [], 1:numElectrodes));
            save(exportLocation, 'EEGData', '-v7.3')
        end
    end
end
