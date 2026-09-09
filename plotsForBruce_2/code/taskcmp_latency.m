function taskcmp_latency()
% TASKCMP_LATENCY  Timing of theta (4-8 Hz) and gamma (25-58 Hz) power in the first
% 1000 ms after inhale onset, paired audiobook vs focus, CONTROL, CP EXCLUDED (n=5;
% CP breathed atypically and had the fewest focus breaths). Four timing metrics:
%   peakLatMs  : latency of the max of a 50-ms-smoothed trace in [50,950] ms
%   riseSlope  : OLS slope (z/s) from 0 ms to the peak
%   centroidMs : center-of-mass latency = sum(t*max(y,0))/sum(max(y,0)) over 0-1000 ms
%   xcorrLagMs : lag at the peak of xcorr(audiobook, focus) within +-500 ms
%                (+ = audiobook LATER than focus). One value per subject (paired->one-sample vs 0).
% Theta time courses from scratch_theta_pac.mat; gamma whole-band mean-z from spectro2.
proj='C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2';
tdir=fullfile(proj,'out','tables'); fdir=fullfile(proj,'out','figs','taskcmp');
sdir=fullfile(proj,'out','spectro2'); fs=500; conds={'audiobook','focus'}; EXCL={'CP'};

% ---- theta time courses per subject ----
Z=load(fullfile(tdir,'scratch_theta_pac.mat')); Sres=Z.Sres; NTS=Z.NTS;
subs=setdiff(unique({Sres.subj},'stable'), EXCL, 'stable');
TH=struct(); for c=1:2, TH.(conds{c})=nan(numel(subs),NTS); end
for si=1:numel(subs), for c=1:2
    m=strcmp({Sres.subj},subs{si})&strcmp({Sres.cond},conds{c}); if ~any(m), continue; end
    TH.(conds{c})(si,:)=mean(cell2mat({Sres(m).thTS}'),1,'omitnan');
end, end
tTh=(0:NTS-1)/fs*1000;

% ---- gamma whole-band mean-z time courses per subject (spectro2) ----
fa=dir(fullfile(sdir,'*OBE*__audiobook.mat')); GG=struct('subj',{{}},'cond',{{}},'tc',{{}}); tG=[];
for i=1:numel(fa)
    base=strrep(fa(i).name,'__audiobook.mat',''); p=strsplit(base,'_'); subj=p{4};
    if any(strcmp(subj,EXCL)), continue; end
    for c=1:2
        f=fullfile(sdir,[base '__' ternary(c==1,'audiobook','focusedBreathing') '.mat']);
        if ~exist(f,'file'), continue; end; o=load(f); o=o.out; if o.nBreaths==0, continue; end
        if isempty(tG), tG=o.tMs; end
        nb=o.nBreaths*o.nBaseSamp; mp=o.mid.sumPow/o.nBreaths; bmu=o.mid.baseSum/nb;
        bv=o.mid.baseSumSq/nb-bmu.^2; bv(bv<=0)=eps; z=(mp-bmu)./sqrt(bv/o.nBreaths); z=max(min(z,10),-10);
        fm=o.mid.freqs>=25&o.mid.freqs<=58; GG.subj{end+1}=subj; GG.cond{end+1}=conds{c}; GG.tc{end+1}=mean(z(fm,:),1);
    end
end
gm=tG>=0&tG<=1000; tGa=tG(gm);
GA=struct(); for c=1:2, GA.(conds{c})=nan(numel(subs),sum(gm)); end
for si=1:numel(subs), for c=1:2
    m=strcmp(GG.subj,subs{si})&strcmp(GG.cond,conds{c}); if ~any(m), continue; end
    GA.(conds{c})(si,:)=mean(cell2mat(cellfun(@(x)x(gm),GG.tc(m),'uni',0)'),1,'omitnan');
end, end

% ---- per-subject-per-condition metrics ----
R=addband(table(),'theta',TH,tTh,subs,fs);
R=addband(R,'gamma',GA,tGa,subs,fs);
writetable(R, fullfile(tdir,'taskcmp_latency.csv'));

% ---- xcorr lag per subject (audiobook vs focus) ----
LX=table();
LX=addlag(LX,'theta',TH,subs,fs); LX=addlag(LX,'gamma',GA,subs,fs);
writetable(LX, fullfile(tdir,'taskcmp_xcorrlag.csv'));

% ---- stats ----
St=table();
for bnd={'theta','gamma'}, for mm={'peakLatMs','riseSlope','centroidMs'}
    a=R.(mm{1})(strcmp(R.band,bnd{1})&strcmp(R.cond,'audiobook'));
    f=R.(mm{1})(strcmp(R.band,bnd{1})&strcmp(R.cond,'focus'));
    ok=isfinite(a)&isfinite(f); a=a(ok); f=f(ok); d=a-f; n=numel(d);
    dz=mean(d)/std(d); [~,pt]=ttest(a,f); pw=signrank(a,f);
    St=[St; table(bnd(1),mm(1),n,mean(a),mean(f),mean(d),dz,pt,pw, ...
        'VariableNames',{'band','metric','n','mean_audiobook','mean_focus','diff','dz','t_p','wilcox_p'})]; %#ok<AGROW>
end, end
% xcorr lag: one-sample vs 0
for bnd={'theta','gamma'}
    v=LX.xcorrLagMs(strcmp(LX.band,bnd{1})); v=v(isfinite(v)); n=numel(v);
    dz=mean(v)/std(v); [~,pt]=ttest(v); pw=signrank(v);
    St=[St; table(bnd(1),{'xcorrLagMs'},n,NaN,NaN,mean(v),dz,pt,pw, ...
        'VariableNames',{'band','metric','n','mean_audiobook','mean_focus','diff','dz','t_p','wilcox_p'})]; %#ok<AGROW>
end
writetable(St, fullfile(tdir,'taskcmp_latency_stats.csv')); disp(St);

% ---- plots ----
fig=figure('Position',[20 20 1400 900],'Color','w','Visible','off');
paired_panel(1, R, St, 'theta','peakLatMs','theta peak latency (ms)');
paired_panel(2, R, St, 'theta','centroidMs','theta centroid (ms)');
paired_panel(3, R, St, 'gamma','peakLatMs','gamma peak latency (ms)');
paired_panel(4, R, St, 'gamma','centroidMs','gamma centroid (ms)');
subplot(3,2,5); lagbar(LX,'theta',St); subplot(3,2,6); lagbar(LX,'gamma',St);
sgtitle('Timing of theta & gamma power, first 1000 ms (control, CP excluded, n=5)','FontWeight','bold');
exportgraphics(fig, fullfile(fdir,'latency_slope.png'),'Resolution',140); close(fig);
fprintf('wrote taskcmp_latency.csv, taskcmp_xcorrlag.csv, taskcmp_latency_stats.csv, latency_slope.png\n');
end

function R=addband(R,band,TC,t,subs,fs)
sm=round(0.05*fs); win=t>=50 & t<=950; wpos=t>=0 & t<=1000;
for c={'audiobook','focus'}, M=TC.(c{1});
  for si=1:numel(subs)
    y=M(si,:); if all(isnan(y)), continue; end
    ys=movmean(y,sm,'omitnan');
    yy=ys; yy(~win)=-inf; [pv,pk]=max(yy); pl=t(pk);
    i0=find(t>=0,1,'first'); idx=i0:pk;
    slope=NaN; if numel(idx)>=3, cf=polyfit(t(idx)/1000, ys(idx),1); slope=cf(1); end
    yp=max(ys(wpos),0); tt=t(wpos); cen=NaN; if sum(yp)>0, cen=sum(tt.*yp)/sum(yp); end
    R=[R; table(string(subs{si}),string(band),string(c{1}),pl,pv,slope,cen, ...
        'VariableNames',{'subject','band','cond','peakLatMs','peakVal','riseSlope','centroidMs'})]; %#ok<AGROW>
  end
end
end
function LX=addlag(LX,band,TC,subs,fs)
maxlag=round(0.5*fs);
for si=1:numel(subs)
    A=TC.audiobook(si,:); F=TC.focus(si,:); if all(isnan(A))||all(isnan(F)), continue; end
    A=A-mean(A,'omitnan'); F=F-mean(F,'omitnan'); A(isnan(A))=0; F(isnan(F))=0;
    [c,lg]=xcorr(A,F,maxlag,'coeff'); [~,ix]=max(c); lag=lg(ix)/fs*1000;   % + = audiobook later
    LX=[LX; table(string(subs{si}),string(band),lag,'VariableNames',{'subject','band','xcorrLagMs'})]; %#ok<AGROW>
end
end
function paired_panel(pos, R, St, bnd, mm, ylab)
subplot(3,2,pos); hold on;
a=R.(mm)(strcmp(R.band,bnd)&strcmp(R.cond,'audiobook')); f=R.(mm)(strcmp(R.band,bnd)&strcmp(R.cond,'focus'));
for i=1:numel(a), plot([1 2],[a(i) f(i)],'-o','Color',[.6 .6 .6],'MarkerFaceColor',[.6 .6 .6],'MarkerSize',4); end
plot([1 2],[mean(a,'omitnan') mean(f,'omitnan')],'-o','Color','k','LineWidth',2.5,'MarkerFaceColor','k','MarkerSize',7);
sr=St(strcmp(St.band,bnd)&strcmp(St.metric,mm),:); xlim([.7 2.3]); set(gca,'XTick',[1 2],'XTickLabel',{'audiobook','focus'});
ylabel(ylab); title(sprintf('%s  (dz=%.2f, p=%.2f)', regexprep(ylab,' \(.*',''), sr.dz, sr.t_p),'FontWeight','bold'); set(gca,'FontSize',10);
end
function lagbar(LX, bnd, St)
v=LX.xcorrLagMs(strcmp(LX.band,bnd)); hold on; bar(1:numel(v), v,'FaceColor',[.4 .5 .8]); yline(0,'k-'); yline(mean(v),'r--','LineWidth',1.5);
sr=St(strcmp(St.band,bnd)&strcmp(St.metric,'xcorrLagMs'),:);
ylabel('xcorr lag (ms, + = audiobook later)'); set(gca,'XTick',1:numel(v));
title(sprintf('%s xcorr lag  mean=%.0f ms (dz=%.2f, p=%.2f)', bnd, mean(v), sr.dz, sr.t_p),'FontWeight','bold'); set(gca,'FontSize',10);
end
function y=ternary(c,a,b), if c, y=a; else, y=b; end, end
