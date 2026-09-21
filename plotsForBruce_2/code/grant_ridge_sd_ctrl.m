function grant_ridge_sd_ctrl()
% Control (OBE) breathing: per-breath SD of the tracked 25-58 Hz ridge frequency over the
% FIRST 1000 ms post-inhale (ungated), audiobook vs focusedBreathing. Tests whether ridge-
% frequency variability (frequency modulation) is higher in focused breathing (paired,
% session level). Replicates the pipeline ridge (extract_gamma_session): superlet 22:62 Hz,
% per-frequency whole-epoch z, ridge_track on max(z,0). Also reports a z>2-gated variant.
codeDir=fileparts(mfilename('fullpath')); proj=fileparts(codeDir); repo=fileparts(proj);
addpath(repo); addpath(codeDir);
T=fullfile(proj,'out','tables');
idx=readtable(fullfile(T,'session_index.csv'),'TextType','string');
mbB=readtable(fullfile(T,'macbp_best.csv'),'TextType','string'); mbB=mbB(mbB.task=="breathingTask",:);
fs=500; Fv=(22:1:62).'; c1=3; ord=[3 30]; mult=1; band=[25 58]; penalty=1.0; ridgeBW=2; floorZ=2;
e0=round(-1000/1000*fs); e1=round(4000/1000*fs); tMs=(e0:e1)/fs*1000;
win = tMs>=0 & tMs<=1000;                       % first 1000 ms post-inhale
bIdx=find(Fv>=band(1)&Fv<=band(2)); Fb=Fv(bIdx);
tasks={'audiobook','focusedBreathing'};

pbf=dir(fullfile(proj,'out','gamma','perbreath','*OBE*audiobook.csv'));
sids=unique(string(erase({pbf.name},'__audiobook.csv')));
S=table('Size',[0 3],'VariableTypes',{'string','string','double'},'VariableNames',{'sessID','ab_sd','fb_sd'});
gated=table('Size',[0 3],'VariableTypes',{'string','string','double'},'VariableNames',{'sessID','ab_sd','fb_sd'});
for i=1:numel(sids)
  s=sids(i);
  fp=idx.finalPath(idx.sessID==s & idx.task=="breathingTask"); if isempty(fp), continue; end
  bl=mbB.bestLabel(mbB.sessID==s); if isempty(bl), continue; end
  D=load(char(fp(1))); vv=fieldnames(D); od=D.(vv{1}); D=[];
  labs=od.labels; ci=find(strcmpi(labs,char(bl(1))),1); if isempty(ci), clear od; continue; end
  sig=od.data(ci,:); N=numel(sig); od=[];
  P=slt_power_cont(sig, fs, Fv, c1, ord, mult);   % nF x N
  msd=nan(1,2); msdg=nan(1,2);
  for t=1:2
    cp=fullfile(proj,'out','gamma','perbreath',char(s)+"__"+tasks{t}+".csv");
    if ~isfile(cp), continue; end
    pb=readtable(cp,'TextType','string'); pb=pb(pb.goodBreath==1 & isfinite(pb.onsetSample),:);
    sds=[]; sdsg=[];
    for k=1:height(pb)
      on=pb.onsetSample(k); a=on+e0; b=on+e1; if a<1||b>N, continue; end
      Pb=P(bIdx, a:b); mu=mean(Pb,2); sdv=std(Pb,0,2); sdv(sdv<=0)=eps; Zw=(Pb-mu)./sdv;
      [frW,~,zpW]=ridge_track(max(Zw,0), Fb, penalty, 1, ridgeBW); frW=frW(:).'; zpW=zpW(:).';
      fw=frW(win); if numel(fw)>=3, sds(end+1)=std(fw); end %#ok<AGROW>
      gg=win & (zpW>floorZ); if nnz(gg)>=3, sdsg(end+1)=std(frW(gg)); end %#ok<AGROW>
    end
    if ~isempty(sds), msd(t)=mean(sds); end
    if ~isempty(sdsg), msdg(t)=mean(sdsg); end
  end
  S=[S; {s, msd(1), msd(2)}]; gated=[gated; {s, msdg(1), msdg(2)}]; %#ok<AGROW>
  fprintf('  %-22s  ungated ab=%.3f fb=%.3f | gated ab=%.3f fb=%.3f\n', s, msd(1),msd(2),msdg(1),msdg(2));
  clear sig P;
end
S=renamevars(S,{'ab_sd','fb_sd'},{'audiobook','focusedBreathing'});
gated=renamevars(gated,{'ab_sd','fb_sd'},{'audiobook','focusedBreathing'});
writetable(S, fullfile(T,'grant_ridgeSD_ctrl_0_1000.csv'));

rep=@(lbl,Tb) report(lbl,Tb);
fprintf('\n==== ridge-frequency SD, first 1000 ms post-inhale, control (paired) ====\n');
rep('UNGATED (all ridge points 0-1000 ms)', S);
rep('GATED (ridge z>2 only)', gated);
end

function report(lbl, Tb)
b=Tb(~isnan(Tb.audiobook)&~isnan(Tb.focusedBreathing),:); d=b.focusedBreathing-b.audiobook;
[~,p,~,st]=ttest(b.focusedBreathing,b.audiobook);
dz=mean(d)/std(d);
fprintf('%s: n=%d | audiobook %.3f±%.3f, focus %.3f±%.3f | focus-audiobook=%.3f, dz=%.2f, t(%d)=%.2f, p=%.4f\n', ...
  lbl, height(b), mean(b.audiobook),std(b.audiobook)/sqrt(height(b)), mean(b.focusedBreathing),std(b.focusedBreathing)/sqrt(height(b)), ...
  mean(d), dz, st.df, st.tstat, p);
[~,pw]=deal(0); try, pw=signrank(b.focusedBreathing,b.audiobook); fprintf('   Wilcoxon signed-rank p=%.4f\n',pw); catch, end
end
