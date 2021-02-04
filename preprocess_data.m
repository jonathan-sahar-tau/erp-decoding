function preprocess_data(subjects)
   
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
        fprintf('Subject:\t%d\n',subject);
        
        for condition = C.conditions
            currentFilename= strcat('N300_', subjectName, '_', condition, '_filtered_6Hz', '.vhdr');
            exportFileName = strcat(subjectName, '_', condition, C.data_suffix, '.mat');
            exportLocation = strcat(exportDataLocation, '\', exportFileName);
            EEGData = (pop_loadbv(dataLocation, currentFilename, [], 1:numElectrodes));
            save(exportLocation, 'EEGData', '-v7.3')
        end
    end
end
