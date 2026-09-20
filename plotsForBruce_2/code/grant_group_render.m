function grant_group_render(matPath, taskRow, winMs, climPct, fLim, zeroAboveHz, climFixed, climLo, scope, xlab)
% GRANT_GROUP_RENDER  Render a group OB spectrogram from a saved meanZ .mat (produced by
% grant_rodent_group). Same style as the single-trial figure; per-task color scale so a
% weaker condition's gamma is highlighted on its own scale. taskRow passed explicitly so
% older maps (without a stored taskRow field) render fine.
%   matPath  path to grant_group_meanZ[_<task>].mat  (fields: meanZ,meanRsp,tMsF,F,nP,nTot)
%   taskRow  'cueTask' | 'focusedBreathing' | 'audiobook'
%   winMs    [preMs postMs] (default [100 1000]);  climPct (default 99.5);  fLim (default [25 60])
if nargin<3||isempty(winMs), winMs=[100 1000]; end
if nargin<4||isempty(climPct), climPct=99.5; end
if nargin<5||isempty(fLim), fLim=[25 60]; end
outdir=fullfile(fileparts(fileparts(mfilename('fullpath'))),'out');  % machine-agnostic (.../plotsForBruce_2/out)
S=load(matPath); meanZ=S.meanZ; meanRsp=S.meanRsp; tMsF=S.tMsF; F=S.F; nP=S.nP; nTot=S.nTot;
if nargin>=6 && ~isempty(zeroAboveHz), meanZ(F>zeroAboveHz,:)=0; end   % demo: null line-noise band
lbl=taskRow; if strcmp(taskRow,'cueTask'),lbl='cued sniff'; elseif strcmp(taskRow,'focusedBreathing'),lbl='focused breathing'; end
xl=[-winMs(1) winMs(2)]; dt=tMsF>=xl(1)&tMsF<=xl(2); fm=F>=fLim(1)&F<=fLim(2);
sub=meanZ(fm,dt); climv=prctile(sub(:),climPct); if ~isfinite(climv)||climv<=0, climv=max([sub(:);1]); end
if nargin>=7 && ~isempty(climFixed), climv=climFixed; end   % override with a fixed clim (e.g. match another task)
if nargin<8||isempty(climLo), climLo=0; end                 % lower clim bound
if nargin<9||isempty(scope), scope='control'; end           % 'control' (OBE) | 'all' (OBE+Dupi) | 'dupiS1' | 'dupiS23'
if nargin<10||isempty(xlab), xlab='Time from inhale (ms)'; end
pos=[0.15 0.26 0.68 0.66];
f=figure('Position',[40 40 1260 520],'Color','w','Visible','off');
ax=axes(f,'Position',pos);
imagesc(ax,tMsF,F,meanZ,[climLo climv]); axis(ax,'xy'); colormap(ax,jet); xlim(ax,xl); ylim(ax,fLim); hold(ax,'on');
r=meanRsp(dt); ry=fLim(1)+(r-min(r))/(max(r)-min(r)+eps)*(fLim(2)-fLim(1));
plot(ax,tMsF(dt),ry,'w-','LineWidth',5);
set(ax,'FontSize',30,'LineWidth',3.5); ylabel(ax,'Hz','FontWeight','bold','FontSize',42);
xlabel(ax,xlab,'FontWeight','bold','FontSize',38); xline(ax,0,'w--','LineWidth',3);
cb=colorbar(ax); ax.Position=pos; cb.Position=[0.845 pos(2) 0.024 pos(4)];
cb.LineWidth=3; cb.FontSize=26; cb.Label.String='power (z)'; cb.Label.FontSize=30;
fprintf('  (%s: n=%d, %d breaths)\n', lbl, nP, nTot);
exportgraphics(f, fullfile(outdir,sprintf('grant_rodentOB_GROUP_%s_%s.png',scope,taskRow)),'Resolution',200); close(f);
fprintf('wrote out/grant_rodentOB_GROUP_%s_%s.png (clim=%.2f, n=%d, breaths=%d)\n', scope, taskRow, climv, nP, nTot);
end
