function plot_results(subjects)
if nargin == 0
    %     subjects = [102 104] %:106 108 109]
    %     subjects = [201] %:206]
%     subjects = [11 13:17 19:22 24:27 29 30 35:42]
    subjects = [102 104:106 108:112 114:116 118:120 122]

end
conditionDesc =  ["Congruency"]
labels = [1,2]
    C = Constants();
    subjects = C.subjects;
    nSubjects = length(subjects);
    conditions = ["ConInt", "IncInt"] %, "ConScr", "IncScr"]    labels = [1, 2, 1, 2]; % intact = 1, scrambled = 2

    baseDir = C.baseDir;
    eeglabPath = C.eeglabPath;
    dataLocation = C.dataLocation;
    outputDir = C.outputDir;
    resultsDir = C.resultsDir;
    figuresDir = C.figuresDir;



for c = 1
%     conditions = conditionPairs{c}
    
    for s = 1:numel(subjects)
        sub = subjects(s);
        subjectName = num2str(sub, '%03.f')
        resultsFile = strcat(outputDir,  'mat-files\',...
            subjectName, '_', ...
            strjoin(conditions, '_'), ...
            '_results_all_electrodes',  ...
...         '_results_occipital',  ...
            '.mat');
        tmp =  load(resultsFile);
        svmECOC = tmp.svmECOC;
        results{s} = svmECOC.successRatesTime;
        fprintf('subject %s, max success rate: %d, mean success rate: %d\n', subjectName,  max(svmECOC.successRatesTime), mean(svmECOC.successRatesTime));
    end

        allResults = 100 * cat(1, results{:});
        times = tmp.svmECOC.downsampledTimes;
        means(c,:) = mean(allResults);
        RMSEs(c,:) = sqrt(sum((allResults - mean(allResults)).^2));
        SEs(c,:) = std(allResults)/sqrt(size(allResults, 1));
end
        
for c = 1
%         conditions = conditionPairs{c};
        figure('units','normalized', 'WindowState', 'maximized')

        errorbar(times, means(c,:), SEs(c, :), 'b')
        hold on
        plot(times, means(c,:), 'k')
        plot(times, repmat((1/numel(unique(labels)) * 100), 1, numel(times)), 'm--');
        plot(times, repmat((1/numel(unique(labels)) * 100 * 2), 1, numel(times)), 'm--');
        hold off

        ylim([20 90]);
        xlim([times(1) times(end)]);
        xlabel('Time')
        ylabel('Mean success rate %')
        titleString = sprintf('Means and SEs\nCondition: %s, all electrodes', conditionDesc(c));
        title(titleString)

        figureFileName = sprintf('mean-and-error-%d-%d-%s-all-electrodes.jpg',subjects(1), subjects(end), conditionDesc(c));
        figureFileName = strcat(outputDir, 'figures\latest\', figureFileName);
        print(gcf, figureFileName, '-djpeg', '-r0');
   end
end
