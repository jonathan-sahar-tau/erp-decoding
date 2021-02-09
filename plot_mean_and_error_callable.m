function plot_mean_and_error(C)

%     C = Constants();
    subjects = C.subjects;
    nSubjects = C.nSubjects;
    eeglabPath = C.eeglabPath;
    conditionDesc =  C.conditionDesc;

%     for s = 1:nSubjects
%         sub = subjects(s);
%         subjectName = num2str(sub, '%03.f');
%         resultsFile = strcat(C.resultsDir, ...
%                              subjectName, '_', ...
%                              C.conditionDesc, ...
%                              '_results', ...
%                              C.data_suffix, ...
%                              C.result_suffix, ...
%                              '.mat');
%         load(resultsFile); % into decoder
%         results{s} = decoder.successRatesTime;
% %         fprintf('subject %s, max success rate: %d, mean success rate: %d\n', subjectName,  max(decoder.successRatesTime), mean(decoder.successRatesTime));
%     end

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
        SEs = std(allResults, 1)/sqrt(size(allResults, 1));
        
%         figure('units','normalized', 'WindowState', 'maximized')
%         hold on
        color = rand(1,3);
        er = errorbar(times, means, SEs); %,'color', color);
        plot(times, means, 'k')
        plot(times, repmat((1/nClasses * 100), 1, numel(times)), 'm--');
        plot(times, repmat((1/nClasses * 100 * 2), 1, numel(times)), 'm--');
        er.Color = color;
%         hold off

        ylim([20 90]);
        xlim([times(1) times(end)]);
        xlabel('Time')
        ylabel('Mean success rate %')
%         
%         titleString = sprintf("Means and SEs\nCondition: %s%s%s", conditionDesc, C.data_suffix, C.result_suffix);
%         title(titleString,'Interpreter','none')

%         figureFileName = sprintf('mean-and-error-%d-%d-%s%s%s.jpg',subjects(1), subjects(end), conditionDesc, C.data_suffix, C.result_suffix);
%         figureFileName = strcat(C.figuresDir, 'latest\', figureFileName);
%         print(gcf, figureFileName, '-djpeg', '-r0');
end
