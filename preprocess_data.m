function preprocess_data(subjects)
   

    if nargin == 0
        subjects = [11 13 14 15 16 17 19 20 21];
    end
    
    subjects
    nSubs = length(subjects);

    % parameters to set
    conditions = ["cong_int", "cong_scr", "inc_int"];
    numElectrodes = 64;
    % baseDir("G:\My Drive\MudrikLab020818\Experiments_new\Jonathan\")
    baseDir = "C:\Users\Jonathan\Google Drive\Msc neuroscience\lab\";
    eeglabPath = strcat(baseDir, "software\eeglab2020_0")
    addpath(eeglabPath);

    dataLocation = strcat(baseDir, "data\experiment_data\");
    exportDataLocation = strcat(baseDir, "analysis_scripts\erp-decoding\experiment_data\");

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
