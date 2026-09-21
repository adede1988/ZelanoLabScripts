function grant_ridge_overlay_centered()
% Window-centered overlay of mean gamma-frequency trajectories, audiobook vs focused
% breathing, for Control (AD_2 excluded) and Dupi session-1, from the saved maps
% (grant_breathing_R.mat). EACH participant-condition curve is shifted by ITS OWN mean over
% the -100..+1000 ms window, so every individual line has zero window-mean and must cross 0
% (no line can sit entirely above/below 0). This removes both between-participant baseline
% offsets AND the per-participant condition-level offset, leaving PURE trajectory shape /
% frequency modulation. Group mean = mean of the individual window-centered curves. Two panels.
codeDir=fileparts(mfilename('fullpath')); proj=fileparts(codeDir); repo=fileparts(proj);
addpath(repo); addpath(codeDir); T=fullfile(proj,'out','tables'); Fg=fullfile(proj,'out');
L=load(fullfile(T,'grant_breathing_R.mat')); R=L.R; tMsF=L.tMsF; F=L.F;
band=[25 58]; bIdx=find(F>=band(1)&F<=band(2)); Fb=F(bIdx);
win=tMsF>=-100 & tMsF<=1000; tW=tMsF(win);
tasks={'audiobook','focusedBreathing'};

% --- accumulate per-participant, per-condition mean-frequency curves ---
% pk('group|participant') -> struct with fields ab (cell of session curves), fb (cell)
pk=containers.Map('KeyType','char','ValueType','any');
for i=1:numel(R)
  r=R{i}; if isempty(r), continue; end
  if strcmp(r.group,'control')
    if strcmp(r.sessID,'260326_OBE_NWU_AD_2'), continue; end   % 72 h fasted, excluded
    g='Control';
  elseif strcmp(r.group,'dupiS1'), g='Dupi S1';
  else, continue; end
  k=sprintf('%s|%s', g, r.participant);
  if ~isKey(pk,k), pk(k)=struct('ab',{{}},'fb',{{}}); end
  S=pk(k);
  for t=1:2
    d=r.(tasks{t}); if d.n==0, continue; end
    mZ=max(d.sumZ(bIdx,win)/d.n, 0);              % band mean z (>=0)
    w=mZ./(sum(mZ,1)+eps); fr=sum(Fb.*w,1);       % power-weighted mean freq per time
    if t==1, S.ab{end+1}=fr; else, S.fb{end+1}=fr; end
  end
  pk(k)=S;
end

% --- per participant: session-mean each condition, center by own grand mean (both conds) ---
groups={'Control','Dupi S1'};
CEN=containers.Map('KeyType','char','ValueType','any');   % group|t -> matrix (participants x time)
for gi=1:2, for t=1:2, CEN(sprintf('%s|%d',groups{gi},t))=[]; end, end
ks=keys(pk); allv=[];
for j=1:numel(ks)
  p=split(ks{j},'|'); g=p{1}; S=pk(ks{j});
  ab=[]; fb=[];
  if ~isempty(S.ab), ab=mean(cat(1,S.ab{:}),1); end
  if ~isempty(S.fb), fb=mean(cat(1,S.fb{:}),1); end
  % center EACH condition curve by its OWN window mean -> every line is zero-mean over the
  % window (must cross 0, cannot sit entirely above/below), isolating trajectory shape.
  if ~isempty(ab), abC=ab-mean(ab); CEN(sprintf('%s|1',g))=[CEN(sprintf('%s|1',g)); abC]; allv=[allv abC]; end %#ok<AGROW>
  if ~isempty(fb), fbC=fb-mean(fb); CEN(sprintf('%s|2',g))=[CEN(sprintf('%s|2',g)); fbC]; allv=[allv fbC]; end %#ok<AGROW>
end
m=max(abs(prctile(allv,[2 98]))); m=max(ceil(m),4); yl=[-m m];   % shared, symmetric

% --- plot ---
col=[0 0 0; 0.85 0.1 0.1];   % audiobook=black, focused=red
f=figure('Position',[40 40 1500 640],'Color','w','Visible','off');
tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for gi=1:2
  ax=nexttile; hold(ax,'on'); H=gobjects(1,2); np=0;
  for t=1:2
    M=CEN(sprintf('%s|%d',groups{gi},t));
    for rIdx=1:size(M,1), plot(ax,tW,M(rIdx,:),'-','Color',[col(t,:) 0.30],'LineWidth',1.5); end
    if ~isempty(M), H(t)=plot(ax,tW,mean(M,1),'-','Color',col(t,:),'LineWidth',6); np=max(np,size(M,1)); end
  end
  yline(ax,0,'k:','LineWidth',1.5); xline(ax,0,'k--','LineWidth',2);
  xlim(ax,[tW(1) tW(end)]); ylim(ax,yl); set(ax,'FontSize',22,'LineWidth',2.5);
  title(ax,sprintf('%s (n=%d)',groups{gi},np),'FontSize',24,'FontWeight','bold');
  xlabel(ax,'Time from +200 ms (ms)','FontWeight','bold','FontSize',22);
  if gi==1, ylabel(ax,'\Delta gamma frequency (Hz)','FontWeight','bold','FontSize',23);
    lg=legend(H,{'audiobook','focused'},'Location','northwest','FontSize',20); lg.Box='off'; end
end
title(tl,'Window-centered mean gamma-frequency ridges: audiobook vs focused breathing','FontSize',24,'FontWeight','bold');
exportgraphics(f, fullfile(Fg,'grant_ridge_overlay_centered.png'),'Resolution',180); close(f);
fprintf('wrote grant_ridge_overlay_centered.png\n');
end
