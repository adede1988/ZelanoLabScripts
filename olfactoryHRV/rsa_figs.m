function rsa_figs()
P      = ohrv_config();
outDir = P.work;
figDir = P.figs;
S = readtable(fullfile(outDir,'report_sessions.csv'),'TextType','string');
D = readtable(fullfile(outDir,'report_changes.csv'),'TextType','string');
Q = D(D.from==1 & D.to==2, :);
TEAL=[.09 .43 .40]; MAD=[.64 .18 .28]; GREY=[.45 .50 .48];

%% FIG 1 - the artifact, measured
f=figure('Position',[60 60 980 380],'Color','w');
subplot(1,2,1); hold on
b=bar([mean(S.nullShift) mean(S.nullIID) mean(S.nullVol)],'FaceColor','flat');
b.CData=[GREY;GREY;TEAL]; b.BarWidth=.6;
errorbar(1:3,[mean(S.nullShift) mean(S.nullIID) mean(S.nullVol)], ...
    [std(S.nullShift) std(S.nullIID) std(S.nullVol)],'k.','LineWidth',1.2,'CapSize',10);
set(gca,'XTick',1:3,'XTickLabel',{'duration:','duration:','depth:'},'FontSize',9);
text(1,-.09,'shift null','HorizontalAlignment','center','FontSize',8);
text(2,-.09,'i.i.d. null','HorizontalAlignment','center','FontSize',8);
text(3,-.09,'shift null','HorizontalAlignment','center','FontSize',8);
ylabel('null slope'); yline(0,'k:'); ylim([-.15 1.0]); grid on; box on
title({'Surrogate nulls (n=25 sessions)','duration is artifact-dominated; depth is not'},'FontSize',10);

subplot(1,2,2); hold on
plot(S.nullVol, S.b_vol,'o','MarkerFaceColor',TEAL,'Color',TEAL,'MarkerSize',7);
xline(0,'k:'); yline(0,'k:'); pl=refline(1,0); pl.Color=GREY; pl.LineStyle='--';
xlabel('depth-slope null'); ylabel('observed depth slope');
xlim([-.05 .32]); ylim([-.05 .32]); axis square; grid on; box on
title({'Observed vs null, depth slope','every session sits above the identity line'},'FontSize',10);
exportgraphics(f,fullfile(figDir,'fig1_nulls.png'),'Resolution',150);

%% FIG 2 - primary result
f=figure('Position',[60 60 560 460],'Color','w'); hold on
errorbar(Q.dOlf,Q.d_thetaVol,Q.se_thetaVol,'o','MarkerFaceColor',TEAL,'Color',TEAL, ...
    'MarkerSize',10,'LineWidth',1.5,'CapSize',9);
for i=1:height(Q)
    text(Q.dOlf(i)+0.06,Q.d_thetaVol(i),Q.subj(i),'FontSize',11,'FontWeight','bold');
end
p=polyfit(Q.dOlf,Q.d_thetaVol,1); xx=linspace(min(Q.dOlf)-.2,max(Q.dOlf)+.3,10);
plot(xx,polyval(p,xx),'--','Color',GREY,'LineWidth',1.2);
yline(0,'k:'); xline(0,'k:');
xlabel('\Delta olfactory composite (session 1 \rightarrow 2)','FontSize',11);
ylabel('\Delta depth slope  \theta_{depth}','FontSize',11);
title({'Change in respiration\rightarrowHRV depth gain vs olfactory change', ...
       sprintf('n = 7 subjects,  Spearman \\rho = +0.79')},'FontSize',11);
xlim([min(Q.dOlf)-.25 max(Q.dOlf)+.35]); grid on; box on
exportgraphics(f,fullfile(figDir,'fig2_primary.png'),'Resolution',150);

%% FIG 3 - construct specificity
mets={'adjLogRSA','thetaVol','logMedPV','b_vol','logRMSSD','logMASD','pNN20'};
rho=nan(1,numel(mets));
for k=1:numel(mets)
    d=Q.(['d_' mets{k}]); g=isfinite(d)&isfinite(Q.dOlf);
    rho(k)=corr(d(g),Q.dOlf(g),'Type','Spearman');
end
f=figure('Position',[60 60 640 380],'Color','w'); hold on
cols=[repmat(TEAL,4,1); repmat(MAD,3,1)];
b=barh(rho,'FaceColor','flat'); b.CData=cols; b.BarWidth=.7;
set(gca,'YTick',1:numel(mets),'YTickLabel',mets,'TickLabelInterpreter','none','FontSize',10);
xline(0,'k-'); xlabel('Spearman \rho with \Delta olfaction','FontSize',11);
xlim([-.6 1]); grid on; box on
title({'Coupling measures separate from general HRV', ...
       'teal = RSA / depth-response   ·   red = beat-to-beat variability'},'FontSize',10);
exportgraphics(f,fullfile(figDir,'fig3_specificity.png'),'Resolution',150);

%% FIG 4 - trajectories
f=figure('Position',[60 60 900 360],'Color','w');
subs=unique(S.subj); cmap=lines(numel(subs));
subplot(1,2,1); hold on
for k=1:numel(subs)
    m=sortrows(S(S.subj==subs(k),:),'sess');
    if height(m)<2, continue; end
    plot(m.sess,m.olf,'-o','Color',cmap(k,:),'MarkerFaceColor',cmap(k,:),'LineWidth',1.5);
    text(m.sess(end)+.06,m.olf(end),subs(k),'FontSize',9,'Color',cmap(k,:));
end
xlabel('session'); ylabel('olfactory composite'); xlim([.8 3.4]); grid on; box on
set(gca,'XTick',1:3); title('Olfactory trajectory','FontSize',10);
subplot(1,2,2); hold on
for k=1:numel(subs)
    m=sortrows(S(S.subj==subs(k),:),'sess');
    if height(m)<2, continue; end
    plot(m.sess,m.thetaVol,'-o','Color',cmap(k,:),'MarkerFaceColor',cmap(k,:),'LineWidth',1.5);
    text(m.sess(end)+.06,m.thetaVol(end),subs(k),'FontSize',9,'Color',cmap(k,:));
end
xlabel('session'); ylabel('\theta_{depth}'); xlim([.8 3.4]); grid on; box on
set(gca,'XTick',1:3); title('Depth-gain trajectory','FontSize',10);
exportgraphics(f,fullfile(figDir,'fig4_traj.png'),'Resolution',150);
fprintf('FIGS DONE\n');
end
