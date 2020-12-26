function plot_results(subjects)
if nargin == 0
    subjects = [19 20 21];
end

%     subjects = [11 13 14 15 16 17 19 20 21];
    conditions = ["cong_int", "cong_scr"]; %, "inc_int"];
    % baseDir("G:\My Drive\MudrikLab020818\Experiments_new\Jonathan\erp-decoding\")
    baseDir = "C:\Users\Jonathan\Google Drive\Msc neuroscience\lab\analysis_scripts\erp-decoding\";
    outputDir = strcat(baseDir, "output\");

    for s = 1:numel(subjects)
        sub = subjects(s);
        subjectName = num2str(sub, '%03.f')
        resultsFile = strcat(outputDir, ...
            ... 'prev_runs\', ...
            subjectName, '_', ...
            strjoin(conditions, '_'), ...
            '_results',  ...
            ... '_sliding_window', ...
            '.mat');
        tmp =  load(resultsFile);
        svmECOC = tmp.svmECOC;
        results{s} = svmECOC.successRatesTime;
        fprintf('subject %s, max success rate: %d, mean success rate: %d\n', subjectName,  max(svmECOC.successRatesTime), mean(svmECOC.successRatesTime));
    end

    allResults = cat(1, results{:});
    times = tmp.svmECOC.time;

   
        numPlotsPerFigure = 6;
    plotIdx = repmat(1:numPlotsPerFigure, 1, floor(numel(subjects)/numPlotsPerFigure + 1));
    for i = 1:numel(subjects)
        if  mod(i, numPlotsPerFigure) == 1
            figure('units','normalized')
        end
        subplot(2, numPlotsPerFigure/2, plotIdx(i));
        plot(times, allResults(i, :)*100);
        ylim([0 120]);
        xlim([times(1) times(end)]);
        xlabel('Time')
        ylabel('success rate %')
        titleString = sprintf('Sucess rate, subject %d', subjects(i));
        title(titleString)
    end
end
