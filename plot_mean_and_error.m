function plot_mean_and_error(C)
    if nargin ==  0
        C = Constants();
    end
    
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
    allResults = decodingResults.sucessRates;
    times = decodingResults.downsampledTimes;
    nClasses = decodingResults.nClasses;
    means = mean(allResults,1);
    % calculate standard error
    SEs = std(allResults, 1)/sqrt(size(allResults, 1));

    figure('units','normalized', 'WindowState', 'maximized')

    color = rand(1,3);
    er = errorbar(times, means, SEs, 'b')
    er.Color = color;
    hold on
    plot(times, means, 'k')
    plot(times, repmat((1/nClasses * 100), 1, numel(times)), 'm--');
    plot(times, repmat((1/nClasses * 100 * 2), 1, numel(times)), 'm--');
            hold off

    ylim([20 90]);
    xlim([times(1) times(end)]);
    xlabel('Time')
    ylabel('Mean success rate %')

    titleString = sprintf("Means and SEs\nCondition: %s%s%s", conditionDesc, C.data_suffix, C.result_suffix);
    title(titleString,'Interpreter','none')

    figureFileName = sprintf('mean-and-error-%d-%d-%s%s%s.jpg',subjects(1), subjects(end), conditionDesc, C.data_suffix, C.result_suffix);
    figureFileName = strcat(C.figuresDir, 'latest\', figureFileName);
    print(gcf, figureFileName, '-djpeg', '-r0');
end
