
function [winnerIdx, T, chanSumm] = pick_winner_channel(chanSumm)
% Robust within-group scoring to select winner

n = numel(chanSumm);
chanLabel = strings(n,1);
subID = strings(n,1);
sessNum = nan(n,1);

prevPooled = nan(n,1);
promPooled = nan(n,1);
madFreq    = nan(n,1);
burstProm  = nan(n,1);
burstCov   = nan(n,1);

for i=1:n
    chanLabel(i) = chanSumm(i).chanLabelSafe;
    subID(i)     = chanSumm(i).subIDSafe;
    sessNum(i)   = chanSumm(i).sessNum;

    prevPooled(i) = mean(chanSumm(i).prev_by_epoch, 'omitnan');
    promPooled(i) = chanSumm(i).prom_pooled;
    madFreq(i)    = chanSumm(i).madFreq_pooled;
    burstProm(i)  = chanSumm(i).burstProm_median;
    burstCov(i)   = chanSumm(i).burstCoverage;
end

% Robust z (median/MAD)
zPrev  = robust_z(prevPooled);
zProm  = robust_z(promPooled);
zMadF  = robust_z(madFreq);
zBProm = robust_z(burstProm);
zBCov  = robust_z(burstCov);

score = zPrev + zProm;

% Package table sorted
T = table(chanLabel, subID, sessNum, prevPooled, promPooled, madFreq, burstProm, burstCov, ...
    zPrev, zProm, zMadF, zBProm, zBCov, score);

[~, ord] = sort(score, 'descend', 'MissingPlacement','last');
T = T(ord,:);

% Winner = first row
winnerChanLabel = T.chanLabel(1);

% map back to index in chanSumm
winnerIdx = find(strcmp(chanLabel, winnerChanLabel), 1, 'first');

% Save score into chanSumm (optional)
for i=1:n
    chanSumm(i).score = score(i);
end

end

