function rsa_stab()
% Rank-stability under measurement error. The jackknife removes points; it
% cannot see that two adjacent points may be separated by far less than their
% standard errors. Perturb each change score by its own SE and see how often
% the rank correlation survives.

P      = ohrv_config();

sp     = P.code;
outDir = P.work;
Q = readtable(fullfile(outDir, 'gate2_data.csv'), 'TextType', 'string');
load(fullfile(outDir, 'panel_state.mat'), 'keepM');
rng(11); B = 5000; n = height(Q);

fprintf('rank stability under SE perturbation (B=%d, n=%d)\n\n', B, n);
fprintf('%-12s %8s %8s %8s %10s %10s\n', 'metric','rho_obs','median','p05','P(rho=1)','P(rho>=.8)');
for k = 1:numel(keepM)
    d  = Q.(['d_' keepM{k}]); se = Q.(['se_' keepM{k}]);
    r0 = corr(d, Q.dOlf, 'Type', 'Spearman');
    r  = nan(B,1);
    for b = 1:B
        r(b) = corr(d + se .* randn(n,1), Q.dOlf, 'Type', 'Spearman');
    end
    fprintf('%-12s %+8.3f %+8.3f %+8.3f %10.3f %10.3f\n', ...
        keepM{k}, r0, median(r), prctile(r, 5), mean(r >= 0.999), mean(r >= 0.8));
end

% how separated are adjacent points on the primary metric?
fprintf('\nadjacent-pair separation on logRMSSD (sorted by dOlf):\n');
[~, o] = sort(Q.dOlf);
d = Q.d_logRMSSD(o); se = Q.se_logRMSSD(o); s = Q.subj(o);
for i = 1:n-1
    gap = d(i+1) - d(i); pooled = sqrt(se(i)^2 + se(i+1)^2);
    fprintf('  %-3s -> %-3s   gap=%+7.3f   pooled SE=%6.3f   gap/SE=%5.2f\n', ...
        s(i), s(i+1), gap, pooled, gap/pooled);
end
end
