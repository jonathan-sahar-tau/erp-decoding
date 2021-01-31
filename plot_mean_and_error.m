function plot_mean_and_error(subjects)

    C = Constants();
    subjects = C.subjects;
    nSubjects = C.nSubjects;
    suffix = '';
    eeglabPath = C.eeglabPath;
    conditionDesc =  "Congruency"

    for s = 1:nSubjects
        sub = subjects(s);
        subjectName = num2str(sub, '%03.f')
        resultsFile = strcat(C.outputDir,  'mat-files\',...
            subjectName, '_', ...
            strjoin(C.conditions, '_'), ...
            '_results',  ...
            '.mat');
        tmp =  load(resultsFile);
        svmECOC = tmp.svmECOC;
        results{s} = svmECOC.successRatesTime;
        fprintf('subject %s, max success rate: %d, mean success rate: %d\n', subjectName,  max(svmECOC.successRatesTime), mean(svmECOC.successRatesTime));
    end

        allResults = 100 * cat(1, results{:});
        times = tmp.svmECOC.downsampledTimes;
        means = mean(allResults);
        SEs = std(allResults)/sqrt(size(allResults, 1));
        
        figure('units','normalized', 'WindowState', 'maximized')

        errorbar(times, means, SEs, 'b')
        hold on
        plot(times, means, 'k')
        plot(times, repmat((1/C.nUniqueLables * 100), 1, numel(times)), 'm--');
        plot(times, repmat((1/C.nUniqueLables * 100 * 2), 1, numel(times)), 'm--');
        hold off

        ylim([20 90]);
        xlim([times(1) times(end)]);
        xlabel('Time')
        ylabel('Mean success rate %')
        titleString = sprintf('Means and SEs\nCondition: %s', conditionDesc);
        title(titleString)

        figureFileName = sprintf('mean-and-error-%d-%d-%s.jpg',subjects(1), subjects(end), conditionDesc);
        figureFileName = strcat(C.figuresDir, 'latest\', figureFileName);
        print(gcf, figureFileName, '-djpeg', '-r0');
