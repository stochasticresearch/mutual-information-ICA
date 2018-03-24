clear;
clc;

% Knobs
n   = 1000;             % # samples
T   = [3, 4, 5];        % # periods for each signal
SNR = 10;               % Signal SNR
d   = 3;                % # mixed observations
r   = 3;                % # independent/principal components

% Generate ground truth
t        = @(n,T) linspace(0,1,n) * 2 * pi * T;
Ztrue(1,:) = sin(t(n,T(1)));            % Sinusoid
Ztrue(2,:) = sign(sin(t(n,T(2))));      % Square
Ztrue(3,:) = sawtooth(t(n,T(3)));       % Sawtooth

% Add noise
sigma    = @(SNR,X) exp(-SNR / 20) * (norm(X(:)) / sqrt(numel(X)));
Znoisy = Ztrue + sigma(SNR,Ztrue) * randn(size(Ztrue));

% Generate mixed signals
normRows = @(X) bsxfun(@rdivide,X,sum(X,2));
A = normRows(rand(d,3));
Zmixed = A * Znoisy;

K = 3; 
n_random_initializations = 1; 
random_seed = 1;
plot_figures = 0;
% fHandle = @cim;
% fArgs = {};
fHandle = @KraskovMI_cc_mex;
fArgs = {1};
[Rica, Wica, Rpca, Wpca] = min_mi_estimator_ICA(Zmixed, K, fHandle, fArgs, ...
    n_random_initializations, random_seed, plot_figures);

icaScore = compute_signal_similarities(Ztrue,Wica);
pcaScore = compute_signal_similarities(Ztrue,Wpca);
fprintf('ICA=%0.02f  ||  PCA=%0.02f\n', icaScore, pcaScore);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cm = hsv(max([3, r, d]));
%cm = linspecer(max([3, r, d]));

figure();

% Truth
subplot(2,2,1);
for i = 1:3
%     plot(Znoisy(i,:),'-','Color',cm(i,:)); hold on;
    plot(Znoisy(i,:),'-'); hold on;
end
title('True signals');
axis tight;

% Observations
subplot(2,2,2);
for i = 1:d
%     plot(Zmixed(i,:),'-','Color',cm(i,:)); hold on;
    plot(Zmixed(i,:),'-'); hold on;
end
title('Observed mixed signals');
axis tight;

% MI ICA
subplot(2,2,3);
for i = 1:r
%     plot(Wica(i,:),'-','Color',cm(i,:)); hold on;
    plot(Wica(i,:),'-'); hold on;
end
title(sprintf('[MI-ICA] - Score=%0.02f',icaScore));
axis tight;

% PCA
subplot(2,2,4);
for i = 1:r
%     plot(Wpca(i,:),'-','Color',cm(i,:)); hold on;
    plot(Wpca(i,:),'-'); hold on;
end
title(sprintf('[PCA] - Score=%0.02f',pcaScore));
axis tight;