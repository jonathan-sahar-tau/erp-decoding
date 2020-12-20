% This function does ERP-based decoding for orientation reported in Experiment 1.

% Structure of analysis and names of some variables used in this script

% were borrowed from the analysis script used in Foster et al. (2016).

% Note: Randomization routine included in this analysis (e.g., random permutation)

% can produce results slightly different from results reported in the paper. 

% Gi-Yeul Bae 2017.10.3.

function svmECOC_our_data(subs)
% delete(gcp)
% parpool

subjects = [11 13]; % 14 15 16];
if nargin == 0
    subjects = [11 13 14 15 16 17 19 20 21];
end
subjects
nSubs = length(subs);

% parameters to set

svmECOC.conditions = ["cong_int", "cong_scr", "inc_int"];

svmECOC.nChans = numel(svmECOC.conditions); % # of channels

svmECOC.nBins = svmECOC.nChans; % # of stimulus bins

svmECOC.nIter = 10; % # of iterations

svmECOC.nBlocks = 3; % # of blocks for cross-validation

svmECOC.frequencies = [0 6]; % low pass filter

svmECOC.time = -200:20:800; % time points of interest in the analysis

svmECOC.window = 4; % 1 data point per 4 ms in the preprocessed data

svmECOC.Fs = 512 ; % samplring rate of the preprocessed data for filtering


svmECOC.nElectrodes = 64; % # of electrode included in the analysis

% for brevity in analysis

nChans = svmECOC.nChans;

nBins = svmECOC.nBins;

nIter = svmECOC.nIter;

nBlocks = svmECOC.nBlocks;

freqs = svmECOC.frequencies;

times = svmECOC.time;
 
nElectrodes = svmECOC.nElectrodes;

nSamps = length(svmECOC.time);

Fs = svmECOC.Fs;

t = templateSVM('Verbose', 1)

fieldsToDelete = {'filename', 'filepath', 'subject', 'group', ...
                'condition', 'setname', 'session', 'srate', ...
                'xmin', 'xmax', 'icaact', 'icawinv', ...
                'icasphere', 'icaweights', 'icachansind', 'chanlocs', ...
                'urchanlocs', 'chaninfo', 'ref', 'event', ...
                'urevent', 'eventdescription', 'epoch', 'epochdescription', ...
                'reject', 'stats', 'specdata', 'specicaact', ...
                'splinefile', 'icasplinefile', 'dipfit', 'history', ...
                'saved', 'etc', 'run'};

dataLocationLuck = strcat("/home/jonathan/google_drive/Msc neuroscience/lab", ...
                        "/data/data_luck_2018/Exp1_Data"); % set directory of data set
loadThis = strcat(dataLocationLuck,"/Decoding_DE_",num2str(currentSub),".mat");
load(loadThis)

dataLocation = strcat("/home/jonathan/google_drive/Msc neuroscience/lab", ...
                        "/data/experiment_data"); % set directory of data set
outputDir = strcat(dataLocation, '/'); %, 'output');

outputFile = strcat(outputDir, '/', datestr(now,'mm-dd-yyyy-HH-MM'), '_run-log');
fileId = fopen(outputFile, 'w');

nCVBlocks = 3;
averagedSuccessRates = nan(numel(subjects), 1);
allSuccessRates = nan(nCVBlocks, numel(subjects));
%% Loop through participants
    for subjectIdx = 1:numel(subjects)
        subject = subjects(subjectIdx);
        subjectName = num2str(subject, '%03.f');
        fprintf(fileId, 'Subject:\t%d\n',subject);
        % conditions = ["cong_int", "inc_int", "cong_scr", "inc_scr"]
        minTrialnum = intmax;

        % assign an integer value to each condition, so we can later put all the data in
        % one matrix and reference it with this number.
        enumVal = 1;

        for condition = conditions
            currentFilename= strcat(subjectName, '_', condition, '_bc', '.vhdr');
            tempData = (pop_loadbv(dataLocation, currentFilename, [], 1:64));

            % remove uneeded fields from the struct
            tempData = rmfield(tempData, fieldsToDelete);

            % swap the 1st and 2nd dimensions, effectively transposing the
            % "inner" (electrode) matrices to be nTrials x nSamples -
            % every row is the time course of one trial
            tempData.data = permute(tempData.data, [1,3,2]);
            EEG.(condition) = tempData;
            EEG.(condition).name = condition;
            EEG.(condition).enumVal = enumVal;
            enumVal = enumVal + 1;

            % get the minimal number of trials across conditions
            trialNums(subjectIdx) = EEG.(condition).trials;
            % set directory for decoding results.
        end

        % find the actual # of trials that can accomodate our desired # of averagedTrials,
        nTrialsPerAveragedTrial = floor(min(trialNums)/nAveragedTrials);
        nTrials = nTrialsPerAveragedTrial * nAveragedTrials;

        % normalize all conditions to have the same number of time points in each trial recording
        nSamples = EEG.(conditions(1)).pnts;

        allData = nan(nTrials * nChans, nElectrodes,nSamples)
        for condition = conditions
            % create a mapping for shuffling the trials: the value of index[i]
            % gives the index of trial that was shuffled into the i-th
            % position for this iteration
            EEG.(condition).averagedTrialIndex = repmat((1:nAveragedTrials)',nTrialsPerAveragedTrial, 1);
            enumVal =  EEG.(condition).enumVal;
            trialDataRange = (1 + (enumVal-1) *  nTrials : nTrials * enumVal);
            allData(trialDataRange, :, :) = EEG.(condition).data
        end

    data.eeg = data.eeg(data.channel == 1 | data.channel == 10);
    % where to save decoding output
    saveLocation = pwd; % set directory for decoding results.

    
    % set up locaiton bin of each trial

    channel = data.channel(data.channel == 1 | data.channel == 10); % tip location of sample teardrop

    svmECOC.posBin = channel'; % add to fm structure so it's saved

    posBin = svmECOC.posBin;
    
    % grab EEG data
    eegs = data.eeg(:,ReleventChan,:); 
    
    % set up time points
    tois = ismember(data.time.pre:4:data.time.post,svmECOC.time);
    nTimes = length(tois);
    
    % # of trials
    svmECOC.nTrials = length(posBin); nTrials = svmECOC.nTrials;

    % Preallocate Matrices
    svm_predict = nan(nIter,nSamps,nBlocks,nChans); % a matrix to save prediction from SVM
    tst_target = nan(nIter,nSamps,nBlocks,nChans);  % a matrix to save true target values

    svmECOC.blocks = nan(nTrials,nIter);  % create svmECOC.block to save block assignments

    % low-pass filtering
    filtData = nan(nTrials,nElectrodes,nTimes);


    for c = 1:nElectrodes
          filtData(:,c,:) = squeeze(eegs(:,c,:)); % low pass filter
    end

    size(filtData)
    % for c = 1:nElectrodes
    %       filtData(:,c,:) = eegfilt(squeeze(eegs(:,c,:)),Fs,freqs(1,1),freqs(1,2)); % low pass filter
    % end

    % Loop through each iteration
    tic % start timing iteration loop

    for iter = 1:nIter

        % preallocate arrays

        blocks = nan(size(posBin));

        shuffBlocks = nan(size(posBin));

        % count number of trials within each position bin

        clear binCnt

        for bin = 1:nBins

            binCnt(bin) = sum(posBin == bin); 

        end

        minCnt = min(binCnt); % # of trials for position bin with fewest trials

        nPerBin = floor(minCnt/nBlocks); % max # of trials such that the # of trials for each bin can be equated within each block

        % shuffle trials

        shuffInd = randperm(nTrials)'; % create shuffle index

        shuffBin = posBin(shuffInd); % shuffle trial order

        % take the 1st nPerBin x nBlocks trials for each position bin.

        for bin = 1:nBins;   

            idx = find(shuffBin == bin); % get index for trials belonging to the current bin

            idx = idx(1:nPerBin*nBlocks); % drop excess trials

            x = repmat((1:nBlocks)',nPerBin,1); shuffBlocks(idx) = x; % assign randomly order trials to blocks - actually assign block IDs to data points

        end


        % unshuffle block assignment

        blocks(shuffInd) = shuffBlocks;

        % save block assignment

        svmECOC.blocks(:,iter) = blocks; % block assignment

        svmECOC.nTrialsPerBlock = length(blocks(blocks == 1)); % # of trials per block

        % Average data for each position bin across blocks   

        posBins = 1:nBins;

        blockDat_filtData = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged & filtered EEG data with resampling at 50 Hz

        labels = nan(nBins*nBlocks,1);  % bin labels for averaged & filtered EEG data

        blockNum = nan(nBins*nBlocks,1); % block numbers for averaged & filtered EEG data

        bCnt = 1;

        for ii = 1:nBins

            for iii = 1:nBlocks

                blockDat_filtData(bCnt,:,:) = squeeze(mean(filtData(posBin==posBins(ii) & blocks==iii,:,tois),1));

                labels(bCnt) = ii;

                blockNum(bCnt) = iii;

                bCnt = bCnt+1;

            end

        end

        % Do SVM_ECOC at each time point
        for t = 1:nSamps

            % grab data for timepoint t

            toi = ismember(times,times(t)-svmECOC.window/2:times(t)+svmECOC.window/2);
            
            % average across time window of interest

            dataAtTimeT = squeeze(mean(blockDat_filtData(:,:,toi),3));  

            % Do SVM_ECOC for each block
            for i=1:nBlocks % loop through blocks, holding each out as the test set

                trnl = labels(blockNum~=i); % training labels

                tstl = labels(blockNum==i); % test labels

                trnD = dataAtTimeT(blockNum~=i,:);    % training data

                tstD = dataAtTimeT(blockNum==i,:);    % test data

                % SVM_ECOC
                fprintf("training SVM...")
                mdl = fitcecoc(trnD,trnl, 'Coding','onevsall','Learners','SVM');   %train support vector mahcine
                LabelPredicted = predict(mdl, tstD);       % predict classes for new data
                
                svm_predict(iter,t,i,:) = LabelPredicted;  % save predicted labels
                
                tst_target(iter,t,i,:) = tstl;             % save true target labels

            end % end of block

        end % end of time points

    end % end of iteration

    toc % stop timing the iteration loop


    OutputfName = strcat(saveLocation,'\Orientation_Results_ERPbased_',currentSub,'.mat');
    
    svmECOC.targets = tst_target;
    svmECOC.modelPredict = svm_predict; 

    svmECOC.nBlocks = nBlocks;

    save(OutputfName,'svmECOC','-v7.3');

end % end of subject loop
end % end of condition loop
