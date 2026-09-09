function fig_grant()
% Preliminary-data figures: descriptive, minimal inference.
P      = ohrv_config();
outDir = P.work;
figDir = P.figs;
S = readtable(fullfile(outDir,'report_sessions.csv'),'TextType','string');
D = readtable(fullfile(outDir,'report_changes.csv'),'TextType','string');
Q = D(D.from==1 & D.to==2,:);
TEAL=[.09 .43 .40]; MAD=[.64 .18 .28]; GREY=[.42 .48 .46];
set(0,'DefaultAxesFontSize',11);

%% FIG A - olfactory recovery
f=figure('Position',[40 40 1160 340],'Color','w');
subs=unique(S.subj); cm=lines(numel(subs));
vars={'cue_d','discrimination  (d'')'; 'O15_acc','identification  (proportion)'; 'olf','composite  (z)'};
for v=1:3
    subplot(1,3,v); hold on
    for k=1:numel(subs)
        m=sortrows(S(S.subj==subs(k),:),'sess');
        y=m.(vars{v,1});
        if sum(isfinite(y))<2, continue; end
        plot(m.sess,y,'-o','Color',cm(k,:),'MarkerFaceColor',cm(k,:),'LineWidth',1.8,'MarkerSize',7);
        j=find(isfinite(y),1,'last');
        text(m.sess(j)+.07,y(j),subs(k),'FontSize',9.5,'Color',cm(k,:),'FontWeight','bold');
    end
    xlabel('session'); ylabel(vars{v,2}); set(gca,'XTick',1:3); xlim([.8 3.45]);
    grid on; box on
    if v==2, title('Olfactory recovery is large, and highly variable between patients','FontSize',12); end
end
exportgraphics(f,fullfile(figDir,'grantA_olf.png'),'Resolution',160);

%% FIG B - the two central observations
f=figure('Position',[40 40 1160 470],'Color','w');
panels = { 'd_adjLogRSA', ...
           {'change in RSA at a standard breath','(\Delta log RSA, breath size held fixed)'}, ...
           {'How much cardiac modulation a breath of', 'FIXED size produces'};
           'd_b_vol', ...
           {'change in cardiac response to breath depth','(\Delta slope of log RSA on log inhaled volume)'}, ...
           {'How much EXTRA modulation a deeper', 'breath produces'} };
for s = 1:2
    subplot(1,2,s); hold on
    y = Q.(panels{s,1}); g = isfinite(y) & isfinite(Q.dOlf);
    x = Q.dOlf(g); yy = y(g);
    p = polyfit(x,yy,1); xx = linspace(-.85,1.75,10);
    sp = range(yy)*0.12;
    fill([xx fliplr(xx)],[polyval(p,xx)-sp fliplr(polyval(p,xx)+sp)],TEAL, ...
        'FaceAlpha',.10,'EdgeColor','none');
    plot(xx,polyval(p,xx),'-','Color',TEAL,'LineWidth',1.9);
    plot(x,yy,'o','MarkerFaceColor',TEAL,'Color','w','MarkerSize',15,'LineWidth',1.5);
    lab = Q.subj(g);
    for i=1:numel(x)
        text(x(i),yy(i)-range(yy)*0.085,lab(i),'FontSize',10.5,'FontWeight','bold', ...
            'HorizontalAlignment','center','Color',[.2 .25 .24]);
    end
    yline(0,':','Color',GREY); xline(0,':','Color',GREY);
    xlabel('change in olfactory function  (composite z, session 1\rightarrow2)','FontSize',11.5);
    ylabel(panels{s,2},'FontSize',11.5);
    title(panels{s,3},'FontSize',12.5);
    xlim([-.85 1.8]); ylim([min(yy)-range(yy)*.22 max(yy)+range(yy)*.28]); grid on; box on
    text(-.78, max(yy)+range(yy)*.19, ...
        sprintf('n = %d   \\rho = %.2f   r = %.2f', numel(x), ...
        corr(x,yy,'Type','Spearman'), corr(x,yy)), ...
        'FontSize',11.5,'FontWeight','bold','Color',TEAL);
end
exportgraphics(f,fullfile(figDir,'grantB_main.png'),'Resolution',160);

%% FIG C - convergence + dissociation
f=figure('Position',[40 40 1080 380],'Color','w');
subplot(1,2,1); hold on
lab={'RSA at a standard breath','peak-valley RSA','depth slope (corrected)','depth slope (raw)'};
val=[corr(Q.dOlf,Q.d_adjLogRSA,'Type','Spearman','Rows','complete'), ...
     corr(Q.dOlf,Q.d_logMedPV, 'Type','Spearman','Rows','complete'), ...
     corr(Q.dOlf,Q.d_thetaVol, 'Type','Spearman'), ...
     corr(Q.dOlf,Q.d_b_vol,    'Type','Spearman')];
b=barh(val,'FaceColor',TEAL); b.BarWidth=.6;
set(gca,'YTick',1:4,'YTickLabel',lab); xlim([0 1]); xline(0,'k-');
xlabel('rank correlation with olfactory change'); grid on; box on
title('Four ways of measuring it agree','FontSize',12);

subplot(1,2,2); hold on
plot(Q.dOlf,Q.d_b_vol,'o','MarkerFaceColor',TEAL,'Color',TEAL,'MarkerSize',11);
plot(Q.dOlf,Q.d_b_len,'s','MarkerFaceColor',MAD,'Color',MAD,'MarkerSize',11);
p1=polyfit(Q.dOlf,Q.d_b_vol,1); p2=polyfit(Q.dOlf,Q.d_b_len,1);
xx=linspace(-.7,1.6,10);
plot(xx,polyval(p1,xx),'-','Color',TEAL,'LineWidth',1.8);
plot(xx,polyval(p2,xx),'-','Color',MAD,'LineWidth',1.8);
yline(0,':','Color',GREY);
legend({'breath DEPTH','breath DURATION'},'Location','northeast','FontSize',10.5);
xlabel('change in olfactory function'); ylabel('change in slope');
title({'Depth and duration move oppositely:','the response shifts from time to volume'},'FontSize',12);
grid on; box on; xlim([-.75 1.7]);
exportgraphics(f,fullfile(figDir,'grantC_conv.png'),'Resolution',160);

fprintf('r(depth,olf) Pearson=%.2f Spearman=%.2f\n', corr(Q.dOlf,Q.d_b_vol), corr(Q.dOlf,Q.d_b_vol,'Type','Spearman'));
fprintf('improved: %d of %d\n', sum(Q.dOlf>0), height(Q));
fprintf('GRANT FIGS DONE\n');
end
