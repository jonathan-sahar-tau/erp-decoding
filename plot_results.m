function plot_results(subjects)
    if nargin ==  0
        C = Constants();
    end

    subjects = C.subjects;
    nSubjects = C.nSubjects;
    conditions = C.conditions;
    
    outputFileNameResults = strcat(C.resultsDir, ...
        C.conditionDesc, ...
        '_results', ...
        C.data_suffix, ...
        C.result_suffix, ...
        '.mat');
       
    subject = subjects(1);
    subjectName = num2str(subject, '%03.f');
    outputFileNameDecoder = strcat(C.decoderOutputDir, ...
        subjectName, '_', ...
        C.conditionDesc, ...
        '_decoder_params', ...
        C.data_suffix, ...
        C.result_suffix, ...
        '.mat');
    
    load(outputFileNameResults); %into decodingResults
    load(outputFileNameDecoder); %into decoder
    
    allResults = decodingResults.successRates;
    times = decoder.downsampledTimes;

   
    numPlotsPerFigure = 6;
    plotIdx = repmat(1:numPlotsPerFigure, 1, floor(nSubjects/numPlotsPerFigure + 1));
    for i = 1:nSubjects
        if  mod(i, numPlotsPerFigure) == 1
            figure('units','normalized')
            title("Decoding success rates")
        end
        subplot(2, numPlotsPerFigure/2, plotIdx(i));
        plot(times, allResults(i, :));
        ylim([0 120]);
        xlim([times(1) times(end)]);
        xlabel('Time')
        ylabel('success rate %')
        titleString = sprintf('Sucess rate, subject %d', subjects(i));
        title(titleString)
    end
end
