function plot_mean_and_error_callable(C)
    % This is a version of plot_mean_and_error that does not spawn its own figure, but rather adds another plot onto the CALLER's latest active figure
    subjects = C.subjects;
    nSubjects = C.nSubjects;
    eeglabPath = C.eeglabPath;
    conditionDesc =  C.conditionDesc;

    resultsFile = strcat(C.resultsDir, ...
        C.conditionDesc, ...
        '_results', ...
        C.data_suffix, ...
        C.result_suffix, ...
        '.mat');

    load(resultsFile); % into decodingResults
    allResults = decodingResults.successRates;
    times = decodingResults.downsampledTimes;
    nClasses = decodingResults.nClasses;
    means = mean(allResults,1);
    SEs = std(allResults, 1)/sqrt(size(allResults, 1));

    color = rand(1,3);
    er = errorbar(times, means, SEs);
    plot(times, means, 'k')
    plot(times, repmat((1/nClasses * 100), 1, numel(times)), 'm--');
    plot(times, repmat((1/nClasses * 100 * 2), 1, numel(times)), 'm--');
    er.Color = color;

    ylim([20 90]);
    xlim([times(1) times(end)]);
    xlabel('Time')
    ylabel('Mean success rate %')
end
