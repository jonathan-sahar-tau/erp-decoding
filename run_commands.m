function run_commands()
bPlot = 0;
bUseCustomElectrodes = 0;

C = Constants();

eeglab;

% all parameter options - the whole parameter space
% ---------------------
filters = ["_filtered_6_Hz", "_filtered_30_Hz"];
congruencyElectrodes = ["allElectrodes", "centralElectrodes", "rightCentralElectrodes"];
intactnessElectrodes = ["allElectrodes", "occipitalElectrodes", "rightOccipitalElectrodes"];
descriptions = ["Intactness", "Congruency"];


% overriding with alternative parameter options
% ---------------------
filters = ["_filtered_30_Hz"];
congruencyElectrodes = ["centralElectrodes"];
intactnessElectrodes = ["occipitalElectrodes"];
descriptions = ["Congruency"];


if bUseCustomElectrodes == 1
    load('./output/mat-files/electrode-list-low-variance.mat'); %into lowVarianceElectrodes
    C.customElectrodeGroup = lowVarianceElectrodes;
    C.customElectrodeGroupName = "lowVarianceElectrodes";
end


% for plotting output of multiple runs on the same plot
% -----------------------------------------------------
if bPlot == 1
    figure('units','normalized', 'WindowState', 'maximized')
    titleString = sprintf("Mean prediction success and SE - %s ", descriptions(1));
    title(titleString,'Interpreter','none')
    hold on
end

for filter = filters
    C.data_suffix = filter; % set the data suffix to the filter we used
    for description = descriptions
        C.conditionDesc = description; % set the description
        if  description == "Intactness";
            % for intactness, we want to aggregate: <congIn, incInt> and <congScr, incScr>
            C.bTranslateLabels = 1;
            electrodeGroups = intactnessElectrodes;
        elseif  description == "Congruency";
            C.bTranslateLabels = 0;
            % for congruency, we want to use only intact conditions.
            C.relevantConditionsIdx = [1,2];
            electrodeGroups = congruencyElectrodes;
        end
     %-------------------------------   
        if bUseCustomElectrodes == 1; electrodeGroups = ["customElectrodeGroup"]; end
     %-------------------------------
        for electrodeGroup = electrodeGroups
            % set the result suffix to the name of the electrode subset
            C.result_suffix = strcat("_", electrodeGroup);
            % get the actual list of electrodes based on the name of the group
            C.relevantElectrodes = C.(electrodeGroup);
            %-------------------------------
            if bUseCustomElectrodes == 1;
            % set the result suffix to the name of the electrode subset
                C.result_suffix = strcat("_", C.customElectrodeGroupName);
            % get the actual list of electrodes
                C.relevantElectrodes = C.customElectrodeGroup; 
            end
            %-------------------------------
            C = C.setLabelProperties(); % apply the changes to the constans class (stupid matlab syntax)
            if bPlot == 1
                plot_mean_and_error_callable(C);
            else
                decode_eeg(C);
            end
        end
    end
    
end
end


