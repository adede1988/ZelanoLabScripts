function fig_duration()
% Duration-slope change vs olfactory change, split by interval.
P      = ohrv_config();
outDir = P.work;
figDir = P.figs;
D = readtable(fullfile(outDir,'report_changes.csv'),'TextType','string');
TEAL=[.09 .43 .40]; MAD=[.64 .18 .28]; GREY=[.45 .50 .48];

sets = { D(D.from==1 & D.to==2,:), 'session 1 \rightarrow 2   (~31 days)', TEAL;
         D(D.from==2 & D.to==3,:), 'session 2 \rightarrow 3   (~140 days)', MAD };

f = figure('Position',[60 60 1050 430],'Color','w');
for s = 1:2
    W = sets{s,1}; col = sets{s,3};
    subplot(1,2,s); hold on
    errorbar(W.dOlf, W.d_theta, W.se_theta, 'o', 'MarkerFaceColor',col,'Color',col, ...
        'MarkerSize',10,'LineWidth',1.4,'CapSize',8);
    for i=1:height(W)
        text(W.dOlf(i)+0.07, W.d_theta(i), W.subj(i),'FontSize',11,'FontWeight','bold');
    end
    r  = corr(W.dOlf, W.d_theta, 'Type','Spearman');
    rp = corr(W.dOlf, W.d_theta);
    p  = polyfit(W.dOlf, W.d_theta, 1);
    xx = linspace(min(W.dOlf)-.25, max(W.dOlf)+.35, 10);
    plot(xx, polyval(p,xx), '--','Color',GREY,'LineWidth',1.2);
    yline(0,'k:'); xline(0,'k:');
    xlabel('\Delta olfactory composite','FontSize',11);
    ylabel('\Delta duration slope  \theta_{duration}','FontSize',11);
    title(sprintf('%s\nn = %d,  \\rho = %+.2f   (Pearson %+.2f)', sets{s,2}, height(W), r, rp), ...
        'FontSize',11);
    xlim([min(W.dOlf)-.3 max(W.dOlf)+.45]); ylim([-0.75 1.0]);
    grid on; box on
    fprintf('%s: n=%d  Spearman=%+.3f  Pearson=%+.3f  meanSE=%.3f  var=%.4f  kappa=%.2f\n', ...
        sets{s,2}, height(W), r, rp, mean(W.se_theta), var(W.d_theta), ...
        mean(W.se_theta.^2)/var(W.d_theta));
end
exportgraphics(f, fullfile(figDir,'fig5_duration.png'), 'Resolution', 150);

% same thing for the depth slope, for side-by-side comparison
fprintf('\nfor comparison, DEPTH slope:\n');
for s = 1:2
    W = sets{s,1};
    fprintf('%s: n=%d  Spearman=%+.3f  kappa=%.2f\n', sets{s,2}, height(W), ...
        corr(W.dOlf, W.d_thetaVol,'Type','Spearman'), ...
        mean(W.se_thetaVol.^2)/var(W.d_thetaVol));
end
fprintf('FIG5 DONE\n');
end
