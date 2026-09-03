function fig_grant_thresh()
% Threshold task as bias-calibrated intensity sensitivity: (high - air).
P      = ohrv_config();
outDir = P.work;
figDir = P.figs;
TRACK = 800;
O = readtable(fullfile(outDir,'olfactory_scores.csv'),'TextType','string');
Q = readtable(fullfile(outDir,'thresh_changes.csv'),'TextType','string');
O.subj = extractBetween(O.sessID,"NMH_","_"+digitsPattern);
O.sess = double(extractAfter(O.sessID,"NMH_"+lettersPattern+"_"));
O.subj(O.subj=="TPB")="TB";
O.high = O.thresh_high_cal / TRACK;
TEAL=[.09 .43 .40]; GREY=[.42 .48 .46];
set(0,'DefaultAxesFontSize',11);

f=figure('Position',[40 40 1160 460],'Color','w');

% -- left: sensitivity trajectories
subplot(1,2,1); hold on
subs=unique(O.subj); cm=lines(numel(subs));
used=nan(1,numel(subs)); lastS=nan(1,numel(subs));
for k=1:numel(subs)
    m=sortrows(O(O.subj==subs(k) & ~isnan(O.sess),:),'sess');
    y=m.high; if sum(isfinite(y))<2, continue; end
    plot(m.sess,y,'-o','Color',cm(k,:),'MarkerFaceColor',cm(k,:),'LineWidth',1.8,'MarkerSize',7);
    j=find(isfinite(y),1,'last');
    used(k) = y(j); %#ok<AGROW>
    off = 0;
    for q = 1:k-1                      % stagger labels that would collide
        if abs(used(q) - y(j)) < 0.022 && m.sess(j) == lastS(q), off = off + 0.028; end
    end
    lastS(k) = m.sess(j); %#ok<AGROW>
    text(m.sess(j)+.07, y(j)+off, subs(k),'FontSize',10,'Color',cm(k,:),'FontWeight','bold');
end
yline(0,':','Color',GREY);
xlabel('session'); set(gca,'XTick',1:3); xlim([.8 3.45]);
ylabel({'intensity sensitivity  (high - air)','proportion of rating scale'},'FontSize',11.5);
title({'Odor intensity sensitivity, bias-calibrated','against the blank-air trial'},'FontSize',12.5);
grid on; box on

% -- right: coupling change vs sensitivity change
subplot(1,2,2); hold on
x=Q.dHigh; y=Q.d_b_vol; g=isfinite(x)&isfinite(y);
p=polyfit(x(g),y(g),1); xx=linspace(-.18,.36,10);
fill([xx fliplr(xx)],[polyval(p,xx)-.035 fliplr(polyval(p,xx)+.035)],TEAL,'FaceAlpha',.10,'EdgeColor','none');
plot(xx,polyval(p,xx),'-','Color',TEAL,'LineWidth',1.9);
plot(x(g),y(g),'o','MarkerFaceColor',TEAL,'Color','w','MarkerSize',15,'LineWidth',1.5);
lab=Q.subj(g); xs=x(g); ys=y(g);
for i=1:numel(xs)
    text(xs(i),ys(i)-.030,lab(i),'FontSize',10.5,'FontWeight','bold', ...
        'HorizontalAlignment','center','Color',[.2 .25 .24]);
end
yline(0,':','Color',GREY); xline(0,':','Color',GREY);
xlabel({'change in intensity sensitivity','\Delta (high - air), session 1\rightarrow2'},'FontSize',11.5);
ylabel({'change in cardiac response to breath depth'},'FontSize',11.5);
title({'Coupling tracks the most peripheral','olfactory measure available'},'FontSize',12.5);
xlim([-.18 .38]); ylim([-.30 .29]); grid on; box on
text(-.155,.25,sprintf('n = %d   \\rho = %.2f', sum(g), corr(xs,ys,'Type','Spearman')), ...
    'FontSize',12,'FontWeight','bold','Color',TEAL);
exportgraphics(f,fullfile(figDir,'grantD_thresh.png'),'Resolution',160);
fprintf('rho=%.3f  n=%d\nGRANT THRESH FIG DONE\n', corr(xs,ys,'Type','Spearman'), sum(g));
end
