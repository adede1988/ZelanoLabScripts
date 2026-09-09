function taskcmp_ridge_extrema()
% Per-breath ridge-frequency extrema in the 200 ms -> end-of-breath window,
% canonical ridge pipeline (superlet c1=3 ord[3 30] F=22:62; per-freq whole-window
% z Zw; ridge_track on 25-58 Hz). Control only (CP excluded), baseline vs ATB.
% Per breath (window = [200 ms, endTim]): ridge freq range (max-min), and the
% latency of the min and of the max as % of the breath. Aggregate per session x
% condition (median), paired across 8 sessions. Caches per-breath ridge traces.
proj='C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2';
addpath(fullfile(proj,'code')); addpath('C:/Users/Adam/Documents/GitHub/ZelanoLabScripts');
addpath('C:/Users/Adam/Documents/GitHub/Superlets/matlab-pure');
tdir=fullfile(proj,'out','tables'); fdir=fullfile(proj,'out','figs','taskcmp');
F=22:1:62; C1=3; ORD=[3 30]; MULT=1; RB=[25 58]; PEN=1.0; RBW=2; fs=500; EXCL={'CP'};
bIdx=find(F>=RB(1)&F<=RB(2)); Fb=F(bIdx)';
e0=round(-0.5*fs); e1=round(4.0*fs); tMs=(e0:e1)/fs*1000;   % epoch -500..4000 ms
best=readtable(fullfile(tdir,'macbp_best.csv'),'TextType','string');
best=best(best.cohort=="OBE" & best.task=="breathingTask",:);
base='R:/Neurology/Zelano_Lab/Lab_Common/OBEControl/';
isAudio=@(s) strcmpi(s,'audio')|strcmpi(s,'audiobook');
isFocus=@(s) strcmpi(s,'focus')|strcmpi(s,'naturalFocus')|strcmpi(s,'slowFocus')|strcmpi(s,'focusedBreathing');
PBt=table(); traces=struct('sess',{{}},'cond',{{}},'fr',{{}}); traces.breathMs=[];
for r=1:height(best)
  id=char(best.sessID(r)); part=char(best.participant(r)); bl=char(best.bestLabel(r));
  if any(strcmpi(part,EXCL)), continue; end
  dd=dir(fullfile(base,id,'preProc','*breathing*.mat')); if isempty(dd), continue; end
  od=[]; for a=1:3, try, S=load(fullfile(dd(1).folder,dd(1).name)); f=fieldnames(S); od=S.(f{1}); clear S; break; catch e, pause(3); end, end
  if isempty(od), fprintf('%s LOAD FAIL\n',id); continue; end
  labs=cellfun(@(x)char(string(x)),od.labels,'uni',0); ci=find(strcmpi(labs,bl),1); if isempty(ci), continue; end
  x=double(od.data(ci,:)); x=fillmissing(x,'linear'); x=fillmissing(x,'nearest'); N=numel(x);
  P=slt_power_cont(x, fs, F, C1, ORD, MULT);
  bd=od.behDat; io=round(double(bd.finalOnset)); en=round(double(bd.endTim)); tc=strtrim(string(getcol(bd,'task')));
  nA=0; nF=0;
  for k=1:numel(io)
    s0=io(k)+e0; s1=io(k)+e1; if isnan(io(k))||isnan(en(k))||s0<1||s1>N, continue; end
    breathMs=(en(k)-io(k))/fs*1000; if breathMs<500, continue; end
    Pb=P(bIdx, s0:s1); Zw=(Pb-mean(Pb,2))./max(std(Pb,0,2),eps);
    fr=ridge_track(max(Zw,0), Fb, PEN, 1, RBW); fr=fr(:).';
    w=tMs>=200 & tMs<=min(breathMs,4000); if nnz(w)<10, continue; end
    frw=fr(w); tw=tMs(w);
    [mx,ix]=max(frw); [mn,in]=min(frw);
    cond=''; if isAudio(tc(k)), cond='baseline'; nA=nA+1; elseif isFocus(tc(k)), cond='ATB'; nF=nF+1; else, continue; end
    PBt=[PBt; table(string(part),string(cond),mx-mn,tw(in),tw(ix),100*tw(in)/breathMs,100*tw(ix)/breathMs,breathMs, ...
      'VariableNames',{'sess','cond','freqRange','minLat_ms','maxLat_ms','minPct','maxPct','breathMs'})]; %#ok<AGROW>
    traces.sess{end+1}=part; traces.cond{end+1}=cond; traces.fr{end+1}=fr; traces.breathMs(end+1)=breathMs; %#ok<AGROW>
  end
  fprintf('%s (%s): %d baseline, %d ATB\n', id, part, nA, nF);
end
writetable(PBt, fullfile(tdir,'taskcmp_ridge_extrema_perbreath.csv'));
save(fullfile(tdir,'taskcmp_ridge_extrema.mat'),'traces','tMs','-v7.3');

% ---- session-level (median per session x condition), paired ----
mets={'freqRange','minPct','maxPct','minLat_ms','maxLat_ms'};
subs=unique(PBt.sess,'stable'); Sess=table(); St=table();
for mc=mets
  a=[]; f=[];
  for i=1:numel(subs)
    va=PBt.(mc{1})(PBt.sess==subs(i) & PBt.cond=="baseline"); vf=PBt.(mc{1})(PBt.sess==subs(i) & PBt.cond=="ATB");
    if isempty(va)||isempty(vf), continue; end
    a(end+1)=median(va,'omitnan'); f(end+1)=median(vf,'omitnan'); %#ok<AGROW>
    if strcmp(mc{1},'freqRange'), Sess=[Sess; table(subs(i),median(va,'omitnan'),median(vf,'omitnan'),'VariableNames',{'sess','base_range','ATB_range'})]; end %#ok<AGROW>
  end
  d=a-f; n=numel(d); dz=mean(d)/std(d); [~,pt]=ttest(a,f); pw=signrank(a,f);
  St=[St; table(mc(1),n,mean(a),mean(f),mean(d),dz,pt,pw,'VariableNames',{'metric','n','mean_baseline','mean_ATB','diff','dz','t_p','wilcox_p'})]; %#ok<AGROW>
end
disp(St); writetable(St, fullfile(tdir,'taskcmp_ridge_extrema_stats.csv'));

% ---- boxplots ----
GREY=[0.62 0.62 0.62]; GREEN=[0.133 0.545 0.133];
panels={'freqRange','ridge freq range (Hz), 200ms->end'; 'minPct','min latency (% of breath)'; 'maxPct','max latency (% of breath)'};
fig=figure('Position',[20 20 1350 560],'Color','w','Visible','off');
for p=1:3, subplot(1,3,p); hold on; mc=panels{p,1};
  a=[]; f=[];
  for i=1:numel(subs), va=PBt.(mc)(PBt.sess==subs(i)&PBt.cond=="baseline"); vf=PBt.(mc)(PBt.sess==subs(i)&PBt.cond=="ATB");
    if isempty(va)||isempty(vf), continue; end
    aa=median(va,'omitnan'); ff=median(vf,'omitnan'); a(end+1)=aa; f(end+1)=ff; %#ok<AGROW>
    plot([1 2],[aa ff],'-o','Color',[.6 .6 .6],'MarkerFaceColor',[.6 .6 .6],'MarkerSize',5); end
  boxchart(ones(numel(a),1),a(:),'BoxFaceColor',GREY); boxchart(2*ones(numel(f),1),f(:),'BoxFaceColor',GREEN);
  sr=St(strcmp(St.metric,mc),:); set(gca,'XTick',[1 2],'XTickLabel',{'baseline','ATB'},'FontSize',12,'LineWidth',1.1);
  ylabel(panels{p,2},'FontSize',12,'FontWeight','bold');
  title(sprintf('dz=%.2f, p=%.2f', sr.dz, sr.t_p),'FontSize',12);
end
sgtitle('Control: ridge-frequency extrema in 200ms->breath-end, baseline vs ATB (paired, n=8)','FontSize',14,'FontWeight','bold');
exportgraphics(fig, fullfile(fdir,'ridge_extrema.png'),'Resolution',150); close(fig);
fprintf('wrote ridge_extrema_perbreath.csv, ridge_extrema_stats.csv, ridge_extrema.png, ridge_extrema.mat\n');
end
function s=getcol(bd,nm), n=height(bd); s=repmat({''},n,1); if ~any(strcmp(bd.Properties.VariableNames,nm)), return; end
  v=bd.(nm); for i=1:n, try, if iscell(v), e=v{i}; else, e=v(i); end, t=string(e); if ~isscalar(t), if isempty(t), t=""; else, t=t(1); end, end, s{i}=char(t); catch, end, end, end
