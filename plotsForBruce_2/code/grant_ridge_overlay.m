function grant_ridge_overlay()
% Overlay of participant mean gamma-frequency trajectories (25-58 Hz power-weighted mean
% frequency of each participant's MEAN z-map), audiobook vs focused breathing, for Control
% (AD_2 excluded) and Dupi session-1, from the saved maps (grant_breathing_R.mat). Two panels.
codeDir=fileparts(mfilename('fullpath')); proj=fileparts(codeDir); repo=fileparts(proj);
addpath(repo); addpath(codeDir); T=fullfile(proj,'out','tables'); Fg=fullfile(proj,'out');
L=load(fullfile(T,'grant_breathing_R.mat')); R=L.R; tMsF=L.tMsF; F=L.F;
band=[25 58]; bIdx=find(F>=band(1)&F<=band(2)); Fb=F(bIdx);
win=tMsF>=-100 & tMsF<=1000; tW=tMsF(win);
tasks={'audiobook','focusedBreathing'};

key=containers.Map('KeyType','char','ValueType','any');
for i=1:numel(R)
  r=R{i}; if isempty(r), continue; end
  if strcmp(r.group,'control'), if strcmp(r.sessID,'260326_OBE_NWU_AD_2'), continue; end, g='Control';
  elseif strcmp(r.group,'dupiS1'), g='Dupi S1'; else, continue; end
  for t=1:2
    d=r.(tasks{t}); if d.n==0, continue; end
    mZ=max(d.sumZ(bIdx,win)/d.n, 0);              % band mean z (>=0)
    w=mZ./(sum(mZ,1)+eps); fr=sum(Fb.*w,1);       % power-weighted mean freq per time
    k=sprintf('%s|%d|%s', g, t, r.participant);
    if ~isKey(key,k), key(k)=[]; end; key(k)=[key(k); fr];
  end
end

col=[0 0 0; 0.85 0.1 0.1];   % audiobook=black, focused=red
groups={'Control','Dupi S1'};
f=figure('Position',[40 40 1500 640],'Color','w','Visible','off');
tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for gi=1:2
  ax=nexttile; hold(ax,'on'); H=gobjects(1,2);
  for t=1:2
    ks=keys(key); M=[];
    for j=1:numel(ks)
      p=split(ks{j},'|');
      if strcmp(p{1},groups{gi}) && str2double(p{2})==t
        fr=mean(key(ks{j}),1); M=[M; fr]; %#ok<AGROW>
        plot(ax,tW,fr,'-','Color',[col(t,:) 0.30],'LineWidth',1.5);
      end
    end
    if ~isempty(M), H(t)=plot(ax,tW,mean(M,1),'-','Color',col(t,:),'LineWidth',6); end
  end
  xline(ax,0,'k--','LineWidth',2); xlim(ax,[tW(1) tW(end)]); ylim(ax,[36 48]);
  set(ax,'FontSize',22,'LineWidth',2.5);
  title(ax,sprintf('%s (n=%d)',groups{gi},groupN(key,groups{gi})),'FontSize',24,'FontWeight','bold');
  xlabel(ax,'Time from +200 ms (ms)','FontWeight','bold','FontSize',22);
  if gi==1, ylabel(ax,'Gamma frequency (Hz)','FontWeight','bold','FontSize',24);
    lg=legend(H,{'audiobook','focused'},'Location','northwest','FontSize',20); lg.Box='off'; end
end
title(tl,'Participant mean gamma-frequency ridges: audiobook vs focused breathing','FontSize',24,'FontWeight','bold');
exportgraphics(f, fullfile(Fg,'grant_ridge_overlay.png'),'Resolution',180); close(f);
fprintf('wrote grant_ridge_overlay.png\n');
end

function n=groupN(key,g)
ks=keys(key); ps={};
for j=1:numel(ks), p=split(ks{j},'|'); if strcmp(p{1},g), ps{end+1}=p{3}; end, end %#ok<AGROW>
n=numel(unique(ps));
end
