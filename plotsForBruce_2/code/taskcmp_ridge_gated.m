function taskcmp_ridge_gated()
% Ridge-frequency sweep GATED on supra-threshold gamma (ridge z-power zp>3),
% canonical pipeline (superlet c1=3 ord[3 30] F=22:62; per-freq whole-window z Zw;
% [fr,~,zp]=ridge_track on 25-58 Hz). Control only (CP excluded), baseline vs ATB.
% Time-resolved gated sweep = at each t, mean ridge freq over breaths with zp(t)>3.
% Also mean gated freq per breath (over the breath's z>3 points). Caches fr+zp per breath.
proj='C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2';
addpath(fullfile(proj,'code')); addpath('C:/Users/Adam/Documents/GitHub/ZelanoLabScripts');
addpath('C:/Users/Adam/Documents/GitHub/Superlets/matlab-pure');
tdir=fullfile(proj,'out','tables'); fdir=fullfile(proj,'out','figs','taskcmp');
F=22:1:62; C1=3; ORD=[3 30]; MULT=1; RB=[25 58]; PEN=1.0; RBW=2; fs=500; EXCL={'CP'}; ZT=3;
bIdx=find(F>=RB(1)&F<=RB(2)); Fb=F(bIdx)';
e0=round(-0.5*fs); e1=round(4.0*fs); tMs=(e0:e1)/fs*1000; nT=numel(tMs);
best=readtable(fullfile(tdir,'macbp_best.csv'),'TextType','string');
best=best(best.cohort=="OBE" & best.task=="breathingTask",:);
base='R:/Neurology/Zelano_Lab/Lab_Common/OBEControl/';
isAudio=@(s) strcmpi(s,'audio')|strcmpi(s,'audiobook');
isFocus=@(s) strcmpi(s,'focus')|strcmpi(s,'naturalFocus')|strcmpi(s,'slowFocus')|strcmpi(s,'focusedBreathing');
GREY=[0.62 0.62 0.62]; GREEN=[0.133 0.545 0.133];
cache=struct('sess',{{}},'cond',{{}},'fr',{{}},'zp',{{}},'breathMs',[]);
% per session x condition: gated-sweep sum/count over breaths, and scalar mean gated freq
sw=containers.Map('KeyType','char','ValueType','any');
scalarGF=table('Size',[0 4],'VariableTypes',{'string','string','double','double'},'VariableNames',{'sess','cond','gatedFreq','nBreath'});
for r=1:height(best)
  id=char(best.sessID(r)); part=char(best.participant(r)); bl=char(best.bestLabel(r));
  if any(strcmpi(part,EXCL)), continue; end
  dd=dir(fullfile(base,id,'preProc','*breathing*.mat')); if isempty(dd), continue; end
  od=[]; for a=1:6, try, S=load(fullfile(dd(1).folder,dd(1).name)); f=fieldnames(S); od=S.(f{1}); clear S; break; catch e, fprintf('[retry %d] ',a); pause(6*a); end, end
  if isempty(od), fprintf('%s LOAD FAIL (all retries)\n',id); continue; end
  labs=cellfun(@(x)char(string(x)),od.labels,'uni',0); ci=find(strcmpi(labs,bl),1); if isempty(ci), continue; end
  x=double(od.data(ci,:)); x=fillmissing(x,'linear'); x=fillmissing(x,'nearest'); N=numel(x);
  P=slt_power_cont(x, fs, F, C1, ORD, MULT);
  bd=od.behDat; io=round(double(bd.finalOnset)); en=round(double(bd.endTim)); tc=strtrim(string(getcol(bd,'task')));
  for cc={'baseline','ATB'}
    if strcmp(cc{1},'baseline'), sel=isAudio(tc); else, sel=isFocus(tc); end
    ks=find(sel(:).'); if isempty(ks), continue; end
    swSum=zeros(1,nT); swCnt=zeros(1,nT); gfAll=[];
    for k=ks
      s0=io(k)+e0; s1=io(k)+e1; if isnan(io(k))||isnan(en(k))||s0<1||s1>N, continue; end
      breathMs=(en(k)-io(k))/fs*1000; if breathMs<500, continue; end
      Pb=P(bIdx, s0:s1); Zw=(Pb-mean(Pb,2))./max(std(Pb,0,2),eps);
      [fr,~,zp]=ridge_track(max(Zw,0), Fb, PEN, 1, RBW); fr=fr(:).'; zp=zp(:).';
      valid = tMs>=0 & tMs<=min(breathMs,4000) & zp>ZT;
      swSum(valid)=swSum(valid)+fr(valid); swCnt(valid)=swCnt(valid)+1;
      % scalar: mean gated freq over this breath's supra-threshold, post-onset, pre-end points
      wnd = tMs>=200 & tMs<=min(breathMs,4000); g = wnd & zp>ZT;
      if nnz(g)>=3, gfAll(end+1)=mean(fr(g)); end %#ok<AGROW>
      cache.sess{end+1}=part; cache.cond{end+1}=cc{1}; cache.fr{end+1}=fr; cache.zp{end+1}=zp; cache.breathMs(end+1)=breathMs; %#ok<AGROW>
    end
    key=[part '_' cc{1}]; sweep=swSum./max(swCnt,1); sweep(swCnt<3)=NaN;
    sw(key)=struct('sweep',sweep,'cnt',swCnt);
    if ~isempty(gfAll), scalarGF=[scalarGF; table(string(part),string(cc{1}),median(gfAll),numel(gfAll),'VariableNames',{'sess','cond','gatedFreq','nBreath'})]; end %#ok<AGROW>
  end
  fprintf('%s (%s) done\n', id, part);
end
save(fullfile(tdir,'taskcmp_ridge_full.mat'),'cache','tMs','-v7.3');

% ---- assemble session gated sweeps ----
subs=unique(cellfun(@(k) k(1:end-9), keys(sw),'uni',0),'stable'); % strip '_baseline'/'_ATB' (both 8/3 chars)
subs=unique(scalarGF.sess,'stable');
SB=nan(numel(subs),nT); SF=nan(numel(subs),nT);
for i=1:numel(subs)
  kb=[char(subs(i)) '_baseline']; kf=[char(subs(i)) '_ATB'];
  if isKey(sw,kb), SB(i,:)=sw(kb).sweep; end
  if isKey(sw,kf), SF(i,:)=sw(kf).sweep; end
end
show=tMs>=0 & tMs<=1500; tt=tMs(show);

fig=figure('Position',[20 20 1400 620],'Color','w','Visible','off');
subplot(1,2,1); hold on;
for i=1:numel(subs), plot(tt,SB(i,show),'-','Color',[GREY 0.4]); plot(tt,SF(i,show),'-','Color',[GREEN 0.4]); end
hB=plot(tt,mean(SB(:,show),1,'omitnan'),'-','Color',GREY*0.5,'LineWidth',3.5);
hA=plot(tt,mean(SF(:,show),1,'omitnan'),'-','Color',GREEN,'LineWidth',3.5);
xline(0,'k-'); xlabel('time from inhale onset (ms)','FontSize',14,'FontWeight','bold');
ylabel('gated ridge frequency (Hz), z>3','FontSize',14,'FontWeight','bold');
legend([hB hA],{'baseline','ATB'},'FontSize',13,'Location','best'); title('Gated ridge-frequency sweep by session','FontSize',14,'FontWeight','bold'); set(gca,'FontSize',12,'LineWidth',1.1);
subplot(1,2,2); hold on; D=SF-SB; nS=sum(all(isfinite(SB(:,show)),2)&all(isfinite(SF(:,show)),2));
mD=mean(D(:,show),1,'omitnan'); sD=std(D(:,show),0,1,'omitnan')./sqrt(sum(isfinite(D(:,show)),1));
fill([tt fliplr(tt)],[mD+sD fliplr(mD-sD)],[0.2 0.5 0.2],'FaceAlpha',0.2,'EdgeColor','none');
plot(tt,mD,'-','Color',[0.13 0.45 0.13],'LineWidth',3.5); yline(0,'k:'); xline(0,'k-');
xlabel('time from inhale onset (ms)','FontSize',14,'FontWeight','bold'); ylabel('\Delta gated ridge freq (ATB-baseline, Hz)','FontSize',14,'FontWeight','bold');
title('Session-wise gated-sweep difference (mean\pmSEM)','FontSize',14,'FontWeight','bold'); set(gca,'FontSize',12,'LineWidth',1.1);
sgtitle('Control: ridge-frequency sweep GATED on ridge z>3, baseline vs ATB','FontSize',15,'FontWeight','bold');
exportgraphics(fig, fullfile(fdir,'ridge_gated_sweep.png'),'Resolution',150); close(fig);

% ---- scalar mean gated freq paired test ----
a=[]; f=[];
for i=1:numel(subs)
  va=scalarGF.gatedFreq(scalarGF.sess==subs(i)&scalarGF.cond=="baseline"); vf=scalarGF.gatedFreq(scalarGF.sess==subs(i)&scalarGF.cond=="ATB");
  if ~isempty(va)&&~isempty(vf), a(end+1)=va(1); f(end+1)=vf(1); end %#ok<AGROW>
end
d=a-f; n=numel(d); dz=mean(d)/std(d); [~,pt]=ttest(a,f); pw=signrank(a,f);
St=table(n,mean(a),mean(f),mean(d),dz,pt,pw,'VariableNames',{'n','mean_baseline_Hz','mean_ATB_Hz','diff','dz','t_p','wilcox_p'});
disp(St); writetable(St, fullfile(tdir,'taskcmp_ridge_gated_stats.csv')); writetable(scalarGF, fullfile(tdir,'taskcmp_ridge_gated_session.csv'));
fprintf('wrote ridge_gated_sweep.png, ridge_gated_stats.csv, ridge_gated_session.csv, ridge_full.mat (scalar gated freq: base %.1f vs ATB %.1f Hz)\n', mean(a), mean(f));
end
function s=getcol(bd,nm), n=height(bd); s=repmat({''},n,1); if ~any(strcmp(bd.Properties.VariableNames,nm)), return; end
  v=bd.(nm); for i=1:n, try, if iscell(v), e=v{i}; else, e=v(i); end, t=string(e); if ~isscalar(t), if isempty(t), t=""; else, t=t(1); end, end, s{i}=char(t); catch, end, end, end
