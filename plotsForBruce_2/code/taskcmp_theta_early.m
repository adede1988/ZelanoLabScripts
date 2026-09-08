function taskcmp_theta_early()
% Theta (4-8 Hz) baseline (audiobook) vs ATB (focused breathing), 0-500 ms after
% inhale onset, CONTROL only (CP excluded), paired across sessions.
%   POWER: superlet band-z from the spectro2 LOW band (mean over 4-8 Hz, mean over
%          0-500 ms) -- the canonical superlet -> myChanZscore pipeline (no reload).
%   ITPC : 4-8 Hz Butterworth+Hilbert phase, per breath locked to inhale onset;
%          inter-trial phase coherence across breaths, 0-500 ms. n-MATCHED between
%          conditions (subsample to min breath count, 30 reps) to remove the ITPC
%          trial-count bias, since baseline/ATB have unequal breath counts.
proj='C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2';
addpath(fullfile(proj,'code')); addpath('C:/Users/Adam/Documents/GitHub/ZelanoLabScripts');
addpath('C:/Users/Adam/Documents/GitHub/Superlets/matlab-pure');
tdir=fullfile(proj,'out','tables'); fdir=fullfile(proj,'out','figs','taskcmp'); sdir=fullfile(proj,'out','spectro2');
THETA=[4 8]; fs=500; EXCL={'CP'}; WIN=[0 500]; NREP=30; rng(7);
GREY=[0.62 0.62 0.62]; GREEN=[0.133 0.545 0.133];
best=readtable(fullfile(tdir,'macbp_best.csv'),'TextType','string');
best=best(best.cohort=="OBE" & best.task=="breathingTask",:);
base='R:/Neurology/Zelano_Lab/Lab_Common/OBEControl/';
isAudio=@(s) strcmpi(s,'audio')|strcmpi(s,'audiobook');
isFocus=@(s) strcmpi(s,'focus')|strcmpi(s,'naturalFocus')|strcmpi(s,'slowFocus')|strcmpi(s,'focusedBreathing');
e0=round(-0.2*fs); e1=round(0.8*fs); tMs=(e0:e1)/fs*1000; wM=tMs>=WIN(1)&tMs<=WIN(2);
POW=table(); ITP=table(); IcurveB=[]; IcurveF=[];
for r=1:height(best)
  id=char(best.sessID(r)); part=char(best.participant(r)); bl=char(best.bestLabel(r));
  if any(strcmpi(part,EXCL)), continue; end
  % ---- theta POWER from spectro2 low band (no reload) ----
  for cc={'baseline','ATB'}
    fn=fullfile(sdir,[id '__' ternIf(strcmp(cc{1},'baseline'),'audiobook','focusedBreathing') '.mat']);
    if ~exist(fn,'file'), continue; end; o=load(fn); o=o.out; if o.nBreaths==0, continue; end
    nb=o.nBreaths*o.nBaseSamp; mp=o.low.sumPow/o.nBreaths; bmu=o.low.baseSum/nb;
    bv=o.low.baseSumSq/nb-bmu.^2; bv(bv<=0)=eps; z=(mp-bmu)./sqrt(bv/o.nBreaths);
    z=min(max(z,-10),10);   % winsorize per spectro2 z-map convention
    fm=o.low.freqs>=THETA(1)&o.low.freqs<=THETA(2); bz=mean(z(fm,:),1);
    w2=o.tMs>=WIN(1)&o.tMs<=WIN(2);
    POW=[POW; table(string(part),string(cc{1}),mean(bz(w2)),max(bz(w2)),'VariableNames',{'sess','cond','thetaZ_mean','thetaZ_max'})]; %#ok<AGROW>
  end
  % ---- theta ITPC from reload (bandpass+Hilbert phase) ----
  dd=dir(fullfile(base,id,'preProc','*breathing*.mat')); if isempty(dd), continue; end
  od=[]; for a=1:6, try, S=load(fullfile(dd(1).folder,dd(1).name)); f=fieldnames(S); od=S.(f{1}); clear S; break; catch e, fprintf('[retry %d] ',a); pause(6*a); end, end
  if isempty(od), fprintf('%s LOAD FAIL\n',id); continue; end
  labs=cellfun(@(x)char(string(x)),od.labels,'uni',0); ci=find(strcmpi(labs,bl),1); if isempty(ci), continue; end
  x=double(od.data(ci,:)); x=fillmissing(x,'linear'); x=fillmissing(x,'nearest'); N=numel(x);
  [bB,aB]=butter(4,THETA/(fs/2),'bandpass'); ph=angle(hilbert(filtfilt(bB,aB,x)));
  bd=od.behDat; io=round(double(bd.finalOnset)); tc=strtrim(string(getcol(bd,'task')));
  PA=[]; PF=[];
  for k=1:numel(io)
    s0=io(k)+e0; s1=io(k)+e1; if isnan(io(k))||s0<1||s1>N, continue; end
    seg=ph(s0:s1);
    if isAudio(tc(k)), PA=[PA; seg]; elseif isFocus(tc(k)), PF=[PF; seg]; end %#ok<AGROW>
  end
  if size(PA,1)<10||size(PF,1)<10, fprintf('%s thin (%d/%d)\n',id,size(PA,1),size(PF,1)); continue; end
  nmin=min(size(PA,1),size(PF,1));
  ia=nmatch_itpc(PA,nmin,NREP); ifc=nmatch_itpc(PF,nmin,NREP);   % n-matched ITPC curves
  IcurveB=[IcurveB; ia]; IcurveF=[IcurveF; ifc]; %#ok<AGROW>
  ITP=[ITP; table(string(part),string('baseline'),mean(ia(wM)),nmin,'VariableNames',{'sess','cond','itpc','nMatched'}); ...
            table(string(part),string('ATB'),mean(ifc(wM)),nmin,'VariableNames',{'sess','cond','itpc','nMatched'})]; %#ok<AGROW>
  fprintf('%s (%s): ITPC base %.3f ATB %.3f (n=%d matched)\n', id, part, mean(ia(wM)), mean(ifc(wM)), nmin);
end
writetable(POW, fullfile(tdir,'taskcmp_theta_early_power.csv'));
writetable(ITP, fullfile(tdir,'taskcmp_theta_early_itpc.csv'));

% ---- paired tests ----
St=table();
St=[St; pairedrow('theta power (band-z mean 0-500)', POW,'thetaZ_mean')];
St=[St; pairedrow('theta power (band-z max 0-500)',  POW,'thetaZ_max')];
St=[St; pairedrow('theta ITPC (0-500, n-matched)',   ITP,'itpc')];
disp(St); writetable(St, fullfile(tdir,'taskcmp_theta_early_stats.csv'));

% ---- figure: power mean, power max, ITPC (paired boxplots) + ITPC curve ----
fig=figure('Position',[20 20 1500 560],'Color','w','Visible','off');
subplot(1,4,1); pbox(POW,'thetaZ_mean','theta power band-z (mean 0-500)',St,'theta power (band-z mean 0-500)',GREY,GREEN);
subplot(1,4,2); pbox(POW,'thetaZ_max','theta power band-z (max 0-500)',St,'theta power (band-z max 0-500)',GREY,GREEN);
subplot(1,4,3); pbox(ITP,'itpc','theta ITPC (0-500, n-matched)',St,'theta ITPC (0-500, n-matched)',GREY,GREEN);
subplot(1,4,4); hold on; tt=tMs;
mb=mean(IcurveB,1); mf=mean(IcurveF,1); sb=std(IcurveB,0,1)/sqrt(size(IcurveB,1)); sf=std(IcurveF,0,1)/sqrt(size(IcurveF,1));
fill([tt fliplr(tt)],[mb+sb fliplr(mb-sb)],GREY,'FaceAlpha',.2,'EdgeColor','none');
fill([tt fliplr(tt)],[mf+sf fliplr(mf-sf)],GREEN,'FaceAlpha',.2,'EdgeColor','none');
hB=plot(tt,mb,'-','Color',GREY*0.6,'LineWidth',2.5); hA=plot(tt,mf,'-','Color',GREEN,'LineWidth',2.5);
xline(0,'k-'); xline(500,'k:'); xlabel('ms from onset','FontSize',12); ylabel('4-8 Hz ITPC','FontSize',12);
legend([hB hA],{'baseline','ATB'},'FontSize',11,'Location','best'); title('ITPC time course','FontWeight','bold'); set(gca,'FontSize',11);
sgtitle('Control: theta (4-8 Hz) 0-500 ms, baseline vs ATB (paired)','FontSize',15,'FontWeight','bold');
exportgraphics(fig, fullfile(fdir,'theta_early.png'),'Resolution',150); close(fig);
fprintf('wrote taskcmp_theta_early_{power,itpc,stats}.csv, theta_early.png\n');
end

function I=nmatch_itpc(Ph,n,R)
% mean ITPC curve over R subsamples of n breaths (n-matched, reduces trial-count bias)
T=size(Ph,2); acc=zeros(1,T);
for r=1:R, idx=randperm(size(Ph,1),n); acc=acc+abs(mean(exp(1i*Ph(idx,:)),1)); end
I=acc/R;
end
function row=pairedrow(name,Tb,mc)
subs=unique(Tb.sess,'stable'); a=[]; f=[];
for i=1:numel(subs), va=Tb.(mc)(Tb.sess==subs(i)&Tb.cond=="baseline"); vf=Tb.(mc)(Tb.sess==subs(i)&Tb.cond=="ATB");
  if ~isempty(va)&&~isempty(vf), a(end+1)=va(1); f(end+1)=vf(1); end, end %#ok<AGROW>
d=a-f; n=numel(d); if n>=3, dz=mean(d)/std(d); [~,pt]=ttest(a,f); pw=signrank(a,f); else, dz=NaN;pt=NaN;pw=NaN; end
row=table(string(name),n,mean(a),mean(f),mean(d),dz,pt,pw,'VariableNames',{'metric','n','mean_baseline','mean_ATB','diff','dz','t_p','wilcox_p'});
end
function pbox(Tb,mc,ylab,St,nm,GREY,GREEN)
subs=unique(Tb.sess,'stable'); hold on; a=[]; f=[];
for i=1:numel(subs), va=Tb.(mc)(Tb.sess==subs(i)&Tb.cond=="baseline"); vf=Tb.(mc)(Tb.sess==subs(i)&Tb.cond=="ATB");
  if isempty(va)||isempty(vf), continue; end
  a(end+1)=va(1); f(end+1)=vf(1); plot([1 2],[va(1) vf(1)],'-o','Color',[.6 .6 .6],'MarkerFaceColor',[.6 .6 .6],'MarkerSize',4); end %#ok<AGROW>
boxchart(ones(numel(a),1),a(:),'BoxFaceColor',GREY); boxchart(2*ones(numel(f),1),f(:),'BoxFaceColor',GREEN);
sr=St(St.metric==nm,:); set(gca,'XTick',[1 2],'XTickLabel',{'baseline','ATB'},'FontSize',11,'LineWidth',1.1);
ylabel(ylab,'FontSize',11,'FontWeight','bold'); title(sprintf('dz=%.2f p=%.2f',sr.dz,sr.t_p),'FontSize',11);
end
function y=ternIf(c,a,b), if c, y=a; else, y=b; end, end
function s=getcol(bd,nm), n=height(bd); s=repmat({''},n,1); if ~any(strcmp(bd.Properties.VariableNames,nm)), return; end
  v=bd.(nm); for i=1:n, try, if iscell(v), e=v{i}; else, e=v(i); end, t=string(e); if ~isscalar(t), if isempty(t), t=""; else, t=t(1); end, end, s{i}=char(t); catch, end, end, end
