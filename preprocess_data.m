function preprocess_data(subjects)
   
    C = Constants();
    subjects = C.subjects;
    nSubs = length(subjects);
    baseDir = C.baseDir;
    eeglabPath = C.eeglabPath;
    dataLocation = C.dataLocation;
    outputDir = C.outputDir;
    numElectrodes = 64;
    exportDataLocation = strcat(baseDir, "analysis_scripts\erp-decoding\experiment_data\");
    addpath(eeglabPath);

%% Loop through participants
    for subjectIdx = 1:numel(subjects)
        subject = subjects(subjectIdx);
        subjectName = num2str(subject, '%03.f');
        fprintf('Subject:\t%d\n',subject);

        % assign an integer value to each condition, so we can later put all the data in
        % one matrix and reference it with this number.
        enumVal = 1;

        for condition = conditions
            currentFilename= strcat(subjectName, '_', condition, '_bc', '.vhdr');
            exportFileName = strcat(subjectName, '_', condition, '_bc', '.mat');
            exportLocation = strcat(exportDataLocation, '\', exportFileName);
            tempData = (pop_loadbv(dataLocation, currentFilename, [], 1:numElectrodes));
            save(exportLocation, 'tempData', '-v7.3')
        end
    end
end
