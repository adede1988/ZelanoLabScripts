function thresh_analysis()
% Threshold task analysed as bias-calibrated intensity sensitivity:
% (med - air) and (high - air), each expressed as a proportion of the
% ~800 px rating track so the units are interpretable.
P      = ohrv_config();
outDir = P.work;
TRACK = 800;
O = readtable(fullfile(outDir,'olfactory_scores.csv'),'TextType','string');
G = readtable(fullfile(outDir,'report_sessions.csv'),'TextType','string');

O.subj = extractBetween(O.sessID,"NMH_","_"+digitsPattern);
O.sess = double(extractAfter(O.sessID,"NMH_"+lettersPattern+"_"));
O.subj(O.subj=="TPB") = "TB";
O.med  = O.thresh_low_cal  / TRACK;     % med  - air, proportion of track
O.high = O.thresh_high_cal / TRACK;     % high - air
O.sens = mean([O.med O.high],2,'omitnan');   % intensity sensitivity index
O.grad = O.high - O.med;                     % concentration gradient

fprintf('=== per-session bias-calibrated intensity sensitivity ===\n');
fprintf('%-6s %5s %9s %9s %9s %9s\n','subj','sess','med-air','high-air','mean','gradient');
X = sortrows(O(~isnan(O.sess),:),{'subj','sess'});
for i=1:height(X)
    fprintf('%-6s %5d %9.3f %9.3f %9.3f %9.3f\n', X.subj(i), X.sess(i), ...
        X.med(i), X.high(i), X.sens(i), X.grad(i));
end

%% change scores, session 1 -> 2, joined to the coupling measures
M = innerjoin(G, O(:,{'subj','sess','med','high','sens','grad'}), 'Keys',{'subj','sess'});
M = sortrows(M,{'subj','sess'});
D = table();
us = unique(M.subj);
for k=1:numel(us)
    m = M(M.subj==us(k),:);
    for j=1:height(m)-1
        D=[D; table(us(k),m.sess(j),m.sess(j+1), ...
            m.med(j+1)-m.med(j), m.high(j+1)-m.high(j), ...
            m.sens(j+1)-m.sens(j), m.grad(j+1)-m.grad(j), ...
            m.b_vol(j+1)-m.b_vol(j), m.thetaVol(j+1)-m.thetaVol(j), ...
            m.b_len(j+1)-m.b_len(j), m.medRSA(j+1)-m.medRSA(j), ...
            'VariableNames',{'subj','from','to','dMed','dHigh','dSens','dGrad', ...
                             'd_b_vol','d_thetaVol','d_b_len','d_medRSA'})]; %#ok<AGROW>
    end
end
Q = D(D.from==1 & D.to==2,:);
fprintf('\n=== change scores, session 1 -> 2  (n=%d) ===\n', height(Q));
disp(Q);

fprintf('\n=== coupling change vs THRESHOLD sensitivity change ===\n');
fprintf('%-12s %22s %10s %10s\n','coupling','olfactory measure','Spearman','Pearson');
cm = {'d_b_vol','depth slope (raw)'; 'd_thetaVol','depth slope (corr)'; 'd_b_len','duration slope'};
om = {'dMed','(med - air)'; 'dHigh','(high - air)'; 'dSens','sensitivity index'; 'dGrad','med->high gradient'};
for c=1:size(cm,1)
    for o=1:size(om,1)
        g = isfinite(Q.(cm{c,1})) & isfinite(Q.(om{o,1}));
        if sum(g)<4, continue; end
        fprintf('%-12s %22s %10.3f %10.3f\n', cm{c,2}, om{o,2}, ...
            corr(Q.(cm{c,1})(g),Q.(om{o,1})(g),'Type','Spearman'), ...
            corr(Q.(cm{c,1})(g),Q.(om{o,1})(g)));
    end
end
writetable(Q, fullfile(outDir,'thresh_changes.csv'));
fprintf('\nTHRESH DONE\n');
end
