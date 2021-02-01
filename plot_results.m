function plot_results(subjects)
    C = Constants();
    subjects = C.subjects;
    nSubjects = C.nSubjects;
    suffix = '';
    conditions = C.conditions;
%     conditions = ["cong_int", "cong_scr"]; %, "inc_int"];
    % baseDir("G:\My Drive\MudrikLab020818\Experiments_new\Jonathan\erp-decoding\")

    for s = 1:nSubjects
        sub = subjects(s);
        subjectName = num2str(sub, '%03.f')
        resultsFile = strcat(C.resultsDir, ...
            subjectName, '_', ...
            strjoin(conditions, '_'), ...
            '_results',  ...
            suffix,  ...
            '.mat');
        load(resultsFile); into "decoder"
        results{s} = decoder.successRatesTime;
    end

    allResults = cat(1, results{:});
    times = decoder.downsampledTimes;

   
    numPlotsPerFigure = 6;
    plotIdx = repmat(1:numPlotsPerFigure, 1, floor(nSubjects/numPlotsPerFigure + 1));
    for i = 1:nSubjects
        if  mod(i, numPlotsPerFigure) == 1
            figure('units','normalized')
            title("Decoding success rates")
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
