function [score, signal_mapping, signal_scores] = compute_signal_similarities(X_ref, X_unmixed)
% we expect X and Y to be matrices of the size [D x N],
% where D is the # of signals, and N is the # of samples
% per signal

num_signals = size(X_ref,1);
scores = zeros(1,num_signals);
signal_mapping = containers.Map('KeyType','int32','ValueType','int32');
signal_scores = containers.Map('KeyType','int32','ValueType','double');
for ii=1:num_signals
    x_ref = X_ref(ii,:)';
    bestCorr = -inf;
    bestIdx = 0;
    for jj=1:num_signals
        if(~any(ismember(cell2mat(signal_mapping.values),jj)))
            % means we haven't mapped this signal yet, so we can consider
            % it
            x_test = X_unmixed(jj,:)';
            score_metric = abs(corr(x_ref,x_test,'type','spearman'));
            if(score_metric>bestCorr)
                bestCorr = score_metric;
                bestIdx = jj;
            end
        end
    end
    signal_mapping(ii) = bestIdx;
    signal_scores(ii) = bestCorr;
end

score = mean(cell2mat(signal_scores.values));

end