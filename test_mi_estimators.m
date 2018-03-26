%% Test the performance of different MI "style" estimators for how well they separate signals w/ ICA

clear;
clc;

if(ispc)
    folder = 'C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\ica\\synthetic\\full_sweep';
elseif(ismac)
    folder = '/Users/Kiran/ownCloud/PhD/sim_results/ica/synthetic/full_sweep';
else
    folder = '/home/kiran/ownCloud/PhD/sim_results/ica/synthetic/full_sweep';
end

rng(12345);

fHandles = {@cim; @KraskovMI_cc_mex; @KraskovMI_cc_mex; @KraskovMI_cc_mex; ...
            @apMI_interface; @vmeMI_interface};
fArgs = {{}; {1}; {6}; {20}; {}; {}};
estimatorNames = {'CIM', 'KNN-1', 'KNN-6', 'KNN-20', 'AP', 'vME'};

estimatorIdx = 1;

fHandle = fHandles{estimatorIdx};
fArg = fArgs{estimatorIdx};
estimatorName = estimatorNames{estimatorIdx};
outputFname = fullfile(folder,sprintf('ica_%s_results.mat',estimatorName));

numMCSims = 100;
nVec   = [50, 100, 200, 500, 1000];             % # samples
T   = [3, 4, 5];        % # periods for each signal
SNRVec = [-15, -10, -5, 0, 5, 10, 15];               % Signal SNR
d   = 3;                % # mixed observations
r   = 3;                % # independent/principal components

icaScoreVec = zeros(length(SNRVec), length(nVec), numMCSims);

% if the file exists, try to "resume" progress rather than redoing already
% computed computations
if(exist(outputFname,'file'))
    load(outputFname,'snrIdx','nIdx','mcSimNum','icaScoreVec');
    snrIdxStart = snrIdx;
    nIdxStart = nIdx;
    mcSimNumStart = mcSimNum + 1; % the value that is saved is what we last completed, so start at the next one
    if(mcSimNumStart>numMCSims)
        mcSimNumStart = 1;
        nIdxStart = nIdxStart + 1;
    end
    if(nIdxStart>length(nVec))
        nIdxStart = 1;
        snrIdxStart = snrIdxStart + 1;
    end
else
    snrIdxStart = 1;
    nIdxStart = 1;
    mcSimNumStart = 1;
end

dispstat('','init'); % One time only initialization
for snrIdx=snrIdxStart:length(SNRVec)
    dispstat(sprintf('Total Progress = %d/%d\n',snrIdx, length(SNRVec)),'keepthis', 'timestamp');
    SNR = SNRVec(snrIdx);
    for nIdx=nIdxStart:length(nVec)
        n = nVec(nIdx);
        % Generate ground truth
        t          = @(n,T) linspace(0,1,n) * 2 * pi * T;
        Ztrue      = zeros(3,n);
        Ztrue(1,:) = sin(t(n,T(1)));            % Sinusoid
        Ztrue(2,:) = sign(sin(t(n,T(2))));      % Square
        Ztrue(3,:) = sawtooth(t(n,T(3)));       % Sawtooth
        
        for mcSimNum=mcSimNumStart:numMCSims
            dispstat(sprintf('MonteCarlo[SNR=%d dB,n=%d] = %d/%d',SNR, n, mcSimNum, numMCSims),'timestamp');
            % Add noise
            sigma    = @(SNR,X) exp(-SNR / 20) * (norm(X(:)) / sqrt(numel(X)));
            Znoisy = Ztrue + sigma(SNR,Ztrue) * randn(size(Ztrue));

            % Generate mixed signals
            normRows = @(X) bsxfun(@rdivide,X,sum(X,2));
            A = normRows(rand(d,3));
            Zmixed = A * Znoisy;

            K = 3; 
            n_random_initializations = 3; 
            random_seed = 1;
            plot_figures = 0;
            
            [Rica, Wica, Rpca, Wpca] = min_mi_estimator_ICA(Zmixed, K, fHandle, fArg, ...
                n_random_initializations, random_seed, plot_figures);

            icaScore = compute_signal_similarities(Ztrue,Wica);
            icaScoreVec(snrIdx,nIdx,mcSimNum) = icaScore;
            
            % save the progress to ownCloud
            save(outputFname);
        end
        mcSimNumStart = 1;
    end
    nIdxStart = 1;
end

%% Analyse the results

clear;
clc;

if(ispc)
    folderNow = 'C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\ica\\synthetic';
elseif(ismac)
    folderNow = '/Users/Kiran/ownCloud/PhD/sim_results/ica/synthetic';
else
    folderNow = '/home/kiran/ownCloud/PhD/sim_results/ica/synthetic';
end

rng(12345);

fHandles = {@cim; @KraskovMI_cc_Mex; @KraskovMI_cc_mex; @KraskovMI_cc_mex; ...
            @apMI_interface; @vmeMI_interface};
fArgs = {{}; {1}; {6}; {20}; {}; {}};
estimatorNames = {'CIM', 'KNN-1', 'KNN-6', 'KNN-20', 'AP', 'vME'};

estimatorsToAnalyze = [1, 2, 3, 4];

% Plot Style 1 -- Plot different SNR's against a given N
snrsToPlot = [0,5,10,15];
nToPlot = 200;
figure;
for ii=1:length(snrsToPlot)    
    snrToPlot = snrsToPlot(ii);

    estimatorNamesOutput = cell(1,length(estimatorsToAnalyze));
    for estimatorIdx=1:length(estimatorsToAnalyze)
        estimator = estimatorsToAnalyze(estimatorIdx);
        estimatorName = estimatorNames{estimatorIdx};
        outputFname = fullfile(folderNow,sprintf('ica_%s_results.mat',estimatorName));
        load(outputFname);

        % find the index associated w/ the n and SNR to plot
        nIdx = find(nVec==nToPlot);
        snrIdx = find(SNRVec==snrToPlot);

        plotData(:,estimatorIdx) = icaScoreVec(snrIdx,nIdx,:);   
        estimatorNamesOutput{estimatorIdx} = estimatorName;
    end

    subplot(2,2,ii);
    boxplot(plotData,estimatorNamesOutput,'Notch','on');
    xlabel('MI Estimator');
    ylabel('ICA Signal Separation Score');
    title(sprintf('n=%d SNR=%d dB',nToPlot,snrToPlot));
end

% Plot Style 2 -- Plot different N's given an SNR
snrToPlot = 15;
nVecToPlot = [100, 200, 500];
figure;
for ii=1:length(nVecToPlot)    
    nToPlot = nVecToPlot(ii);

    estimatorNamesOutput = cell(1,length(estimatorsToAnalyze));
    for estimatorIdx=1:length(estimatorsToAnalyze)
        estimator = estimatorsToAnalyze(estimatorIdx);
        estimatorName = estimatorNames{estimatorIdx};
        outputFname = fullfile(folderNow,sprintf('ica_%s_results.mat',estimatorName));
        load(outputFname);

        % find the index associated w/ the n and SNR to plot
        nIdx = find(nVec==nToPlot);
        snrIdx = find(SNRVec==snrToPlot);

        plotData(:,estimatorIdx) = icaScoreVec(snrIdx,nIdx,:);   
        estimatorNamesOutput{estimatorIdx} = estimatorName;
    end

    subplot(1,3,ii);
    boxplot(plotData,estimatorNamesOutput,'Notch','on');
    if(ii==2)
        xlabel('MI Estimator');
    end
    if(ii==1)
        ylabel('ICA Signal Separation Score');
    end
    title(sprintf('n=%d SNR=%d dB',nToPlot,snrToPlot));
end
