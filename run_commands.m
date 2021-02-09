function run_commands()
C = Constants();
% all parameter options
% ---------------------
filters = ["_filtered_6_Hz", "_filtered_30_Hz"];
congruencyElectrodes = ["allElectrodes", "centralElectrodes", "rightCentralElectrodes"];
intactnessElectrodes = ["allElectrodes", "occipitalElectrodes", "rightOccipitalElectrodes"];
descriptions = ["Intactness", "Congruency"];


% overriding with alternative parameter options
% ---------------------
% filters = ["_filtered_30_Hz"];
% congruencyElectrodes = ["allElectrodes", "rightCentralElectrodes"];
% intactnessElectrodes = ["allElectrodes", "rightOccipitalElectrodes"];
% descriptions = ["Intactness"];
descriptions = ["Congruency"];

bPlot = 1;

% choosing electrodes
% ---------------------
% electrodeGroups = congruencyElectrodes;
% electrodes = intactnessElectrodes;

% for plotting output of multiple runs on the same plot
% -----------------------------------------------------
if bPlot == 1
    figure('units','normalized', 'WindowState', 'maximized')
    titleString = sprintf("Mean prediction success and SE - %s ", descriptions(1));
    title(titleString,'Interpreter','none')
    hold on
end

for filter = filters
    C.data_suffix = filter;
    for description = descriptions
        C.conditionDesc = description;
        if  description == "Intactness";
            C.bTranslateLabels = 1;
            electrodeGroups = intactnessElectrodes;
        elseif  description == "Congruency";
            C.bTranslateLabels = 0;
            C.relevantConditionsIdx = [1,2];
            electrodeGroups = congruencyElectrodes;
        end

        for electrodeGroup = electrodeGroups
            C.result_suffix = strcat("_", electrodeGroup);
            C.relevantElectrodes = C.(electrodeGroup);
            C = C.setLabelProperties();
            if bPlot == 1
                plot_mean_and_error_callable(C);
            else
                decode_eeg(C);
            end
        end
    end
    
end
end


