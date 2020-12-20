function preprocess_data(subjects)
   
    addpath("G:\My Drive\MudrikLab020818\Experiments_new\Jonathan\erp-decoding\software\eeglab2020_0")

    if nargin == 0
        subjects = [11 13 14 15 16 17 19 20 21];
    end
    
    subjects
    nSubs = length(subjects);

    % parameters to set
    conditions = ["cong_int", "cong_scr", "inc_int"];
    numElectrodes = 64;
    addpath("G:\My Drive\MudrikLab020818\Experiments_new\Jonathan\erp-decoding\software\eeglab2020_0")
    dataLocation = "G:\My Drive\MudrikLab020818\Experiments_new\Jonathan\exported data";
    exportDataLocation = "G:\My Drive\MudrikLab020818\Experiments_new\Jonathan\erp-decoding\experiment_data";

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
