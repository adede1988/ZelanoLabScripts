function taskcmp_ampratio()
% gammaHRV's ampRatio_inhaleExhale metric, computed for CONTROL breathing sessions,
% baseline (audiobook) vs ATB (focused breathing). Faithful to
% gammaHRV/matlab/runGammaExtraction.m: best macBP channel (FOOOF 25-58 peak),
% bp = butter(4,[25 58]) filtfilt, env = abs(hilbert(bp)); per breath
%   inhale window = finalOnset -> bm_exhaleOnsets ; exhale window = bm_exhaleOnsets -> endTim
%   ratio = mean(env inhale) / mean(env exhale).
% Aggregates per session x condition (median ratio + geometric-mean via mean log2),
% paired test across sessions. CP excluded (breathing quality).
proj='C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2';
tdir=fullfile(proj,'out','tables'); fdir=fullfile(proj,'out','figs','taskcmp');
GAMMA=[25 58]; EXCL={'CP'};
best=readtable(fullfile(tdir,'macbp_best.csv'),'TextType','string');
best=best(best.cohort=="OBE" & best.task=="breathingTask",:);
base='R:/Neurology/Zelano_Lab/Lab_Common/OBEControl/';
isAudio=@(s) strcmpi(s,'audio')|strcmpi(s,'audiobook');
isFocus=@(s) strcmpi(s,'focus')|strcmpi(s,'naturalFocus')|strcmpi(s,'slowFocus')|strcmpi(s,'focusedBreathing');

Rows=table(); PB=table();
for r=1:height(best)
  id=char(best.sessID(r)); part=char(best.participant(r)); bl=char(best.bestLabel(r));
  if any(strcmpi(part,EXCL)), continue; end
  dd=dir(fullfile(base,id,'preProc','*breathing*.mat')); if isempty(dd), fprintf('%s NO FINAL\n',id); continue; end
  fp=fullfile(dd(1).folder, dd(1).name);
  od=[]; for a=1:3, try, S=load(fp); f=fieldnames(S); od=S.(f{1}); clear S; break; catch e, pause(3); end, end
  if isempty(od), fprintf('%s LOAD FAIL\n',id); continue; end
  fs=500; if isfield(od,'fs'), fs=double(od.fs); end
  labs=cellfun(@(x)char(string(x)),od.labels,'uni',0); ci=find(strcmpi(labs,bl),1);
  if isempty(ci), fprintf('%s no %s\n',id,bl); continue; end
  x=double(od.data(ci,:)); x=fillmissing(x,'linear'); x=fillmissing(x,'nearest'); nSamp=numel(x);
  [bB,aB]=butter(4, GAMMA/(fs/2), 'bandpass'); env=abs(hilbert(filtfilt(bB,aB,x)));
  bd=od.behDat; vn=bd.Properties.VariableNames; g=@(c) double(bd.(c));
  io=round(g('finalOnset'));
  if ismember('bm_exhaleOnsets',vn), eo=round(g('bm_exhaleOnsets')); else, eo=round(g('exMinTim')); end
  ee=round(g('endTim')); tc=strtrim(string(getcol(bd,'task')));
  nb=height(bd); rat=nan(nb,1);
  for k=1:nb
    iS=io(k); iE=eo(k); eS=eo(k); eE=ee(k);
    if any(isnan([iS iE eE])) || iE<=iS+1 || eE<=eS+1 || iS<1 || eE>nSamp, continue; end
    mIn=mean(env(iS:iE)); mEx=mean(env(eS:eE)); if mEx>0, rat(k)=mIn/mEx; end
  end
  for cc={'baseline','ATB'}
    if strcmp(cc{1},'baseline'), sel=isAudio(tc); else, sel=isFocus(tc); end
    v=rat(sel(:) & isfinite(rat)); if numel(v)<3, continue; end
    Rows=[Rows; table(string(id),string(part),cc(1),numel(v),median(v),mean(log2(v)),mean(v), ...
      'VariableNames',{'sessID','participant','cond','nBreath','medRatio','meanLog2','meanRatio'})]; %#ok<AGROW>
    PB=[PB; table(repmat(string(id),numel(v),1),repmat(cc(1),numel(v),1),v,'VariableNames',{'sessID','cond','ratio'})]; %#ok<AGROW>
  end
  fprintf('%s (%s) bl=%s: base med=%.2f ATB med=%.2f\n', id, part, bl, ...
    med_c(Rows,id,'baseline'), med_c(Rows,id,'ATB'));
end
writetable(Rows, fullfile(tdir,'taskcmp_ampratio_session.csv'));
writetable(PB, fullfile(tdir,'taskcmp_ampratio_perbreath.csv'));

% ---- paired tests across sessions (baseline vs ATB) ----
sids=unique(Rows.sessID); St=table();
for mm={'medRatio','meanLog2'}
  a=[]; f=[];
  for i=1:numel(sids)
    ra=Rows.(mm{1})(Rows.sessID==sids(i) & Rows.cond=="baseline");
    rf=Rows.(mm{1})(Rows.sessID==sids(i) & Rows.cond=="ATB");
    if ~isempty(ra)&&~isempty(rf), a(end+1)=ra(1); f(end+1)=rf(1); end %#ok<AGROW>
  end
  d=a-f; n=numel(d); dz=mean(d)/std(d); [~,pt]=ttest(a,f); pw=signrank(a,f);
  St=[St; table(mm(1),n,mean(a),mean(f),mean(d),dz,pt,pw, ...
    'VariableNames',{'metric','n','mean_baseline','mean_ATB','diff','dz','t_p','wilcox_p'})]; %#ok<AGROW>
end
disp(St); writetable(St, fullfile(tdir,'taskcmp_ampratio_stats.csv'));

% ---- paired boxplot (median ratio per session) ----
fig=figure('Position',[40 40 520 640],'Color','w','Visible','off'); hold on;
a=[]; f=[];
for i=1:numel(sids)
  ra=Rows.medRatio(Rows.sessID==sids(i)&Rows.cond=="baseline"); rf=Rows.medRatio(Rows.sessID==sids(i)&Rows.cond=="ATB");
  if ~isempty(ra)&&~isempty(rf), a(end+1)=ra(1); f(end+1)=rf(1); plot([1 2],[ra(1) rf(1)],'-o','Color',[.6 .6 .6],'MarkerFaceColor',[.6 .6 .6],'MarkerSize',5); end %#ok<AGROW>
end
boxchart(ones(numel(a),1), a(:),'BoxFaceColor',[.62 .62 .62]); boxchart(2*ones(numel(f),1), f(:),'BoxFaceColor',[.13 .55 .13]);
yline(1,'k:','LineWidth',1); set(gca,'XTick',[1 2],'XTickLabel',{'baseline','ATB'},'FontSize',14);
ylabel('gamma inhale/exhale amp ratio (median/session)','FontSize',13,'FontWeight','bold');
sr=St(strcmp(St.metric,'medRatio'),:);
title(sprintf('Control gamma inhale/exhale ratio\n(dz=%.2f, p=%.2f, n=%d)', sr.dz, sr.t_p, sr.n),'FontSize',13);
set(gca,'LineWidth',1.1); exportgraphics(fig, fullfile(fdir,'ampratio_inhaleExhale.png'),'Resolution',150); close(fig);
fprintf('wrote taskcmp_ampratio_{session,perbreath,stats}.csv + ampratio_inhaleExhale.png\n');
end
function m=med_c(R,id,c), v=R.medRatio(R.sessID==string(id)&R.cond==string(c)); if isempty(v), m=NaN; else, m=v(1); end, end
function s=getcol(bd,nm), n=height(bd); s=repmat({''},n,1); if ~any(strcmp(bd.Properties.VariableNames,nm)), return; end
  v=bd.(nm); for i=1:n, try, if iscell(v), e=v{i}; else, e=v(i); end, t=string(e); if ~isscalar(t), if isempty(t), t=""; else, t=t(1); end, end, s{i}=char(t); catch, end, end, end
