function fig_raw()
% Raw, uncorrected slopes - no surrogate subtraction - vs olfactory change.
P      = ohrv_config();
outDir = P.work;
figDir = P.figs;
D = readtable(fullfile(outDir,'report_changes.csv'),'TextType','string');
TEAL=[.09 .43 .40]; MAD=[.64 .18 .28]; GREY=[.45 .50 .48];

ivs  = { D(D.from==1 & D.to==2,:), '1 \rightarrow 2  (~31 d)';
         D(D.from==2 & D.to==3,:), '2 \rightarrow 3  (~140 d)' };
mets = { 'd_b_len', 'raw duration slope  b_{length}', TEAL;
         'd_b_vol', 'raw depth slope  b_{depth}',     MAD };

f = figure('Position',[50 50 1040 780],'Color','w');
k = 0;
for m = 1:2
  for s = 1:2
    k = k + 1; W = ivs{s,1}; col = mets{m,3};
    d = W.(mets{m,1});
    subplot(2,2,k); hold on
    plot(W.dOlf, d, 'o','MarkerFaceColor',col,'Color',col,'MarkerSize',10,'LineWidth',1.3);
    for i=1:height(W)
        text(W.dOlf(i)+0.07, d(i), W.subj(i),'FontSize',10,'FontWeight','bold');
    end
    r  = corr(W.dOlf, d, 'Type','Spearman');
    rp = corr(W.dOlf, d);
    p  = polyfit(W.dOlf, d, 1);
    xx = linspace(min(W.dOlf)-.25, max(W.dOlf)+.4, 10);
    plot(xx, polyval(p,xx),'--','Color',GREY,'LineWidth',1.2);
    yline(0,'k:'); xline(0,'k:');
    xlabel('\Delta olfactory composite','FontSize',10);
    ylabel(['\Delta ' mets{m,2}],'FontSize',10);
    title(sprintf('%s   ·   n = %d\n\\rho = %+.2f    Pearson = %+.2f', ...
        ivs{s,2}, height(W), r, rp),'FontSize',11);
    xlim([min(W.dOlf)-.3 max(W.dOlf)+.5]); grid on; box on
    fprintf('%-10s  %-22s n=%d  Spearman=%+.3f  Pearson=%+.3f\n', ...
        mets{m,1}, ivs{s,2}, height(W), r, rp);
  end
end
exportgraphics(f, fullfile(figDir,'fig6_raw.png'), 'Resolution', 150);
fprintf('FIG6 DONE\n');
end
