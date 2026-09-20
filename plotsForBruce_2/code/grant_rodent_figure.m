function grant_rodent_figure(sessID, onsetSample, winMs, climPct, mode, fLim, climRange)
% GRANT_RODENT_FIGURE  Human control cueTask analog of the rodent OB gamma figure
% (Kay/Freeman-style): respiration, OB (macBP) raw LFP, and OB superlet spectrogram
% for a single sniff trial, styled for a small grant panel (bold axes, large fonts).
% Best macBP channel per session (macbp_best.csv).
%
%   sessID       e.g. '250313_OBE_NMH_CS'
%   onsetSample  sniff onset sample (@500 Hz) from the per-breath CSV
%   winMs        [preMs postMs] DISPLAY window around onset (default [100 1000])
%   climPct      spectrogram color-scale upper percentile of the DISPLAYED patch (default 99)
%   mode         'preview' -> stacked PNG to scratchpad; 'final' -> panels 1-4 to out/
%   fLim         [fLo fHi] DISPLAY frequency limits, Hz (default [25 60])
%
% The superlet + within-frequency z are computed over a fixed generous window so the
% normalization is stable; winMs / fLim only crop the view. 'final' writes four PNGs:
% _1resp, _2bulbRaw, _3bulbSpec, and _4stack (the three x-aligned).
if nargin<3||isempty(winMs), winMs=[100 1000]; end
if nargin<4||isempty(climPct), climPct=99; end
if nargin<5||isempty(mode), mode='preview'; end
if nargin<6||isempty(fLim), fLim=[25 60]; end
proj='C:\Users\Adam\Documents\GitHub\ZelanoLabScripts\plotsForBruce_2';
repo='C:\Users\Adam\Documents\GitHub\ZelanoLabScripts';
addpath(repo); addpath(fullfile(proj,'code'));
scratch='C:\Users\Adam\AppData\Local\Temp\claude\C--Users-Adam-Documents-GitHub-ZelanoLabScripts-plotsForBruce-2\2dd116fc-2bdb-4993-a330-8afe33295cdb\scratchpad';
T=fullfile(proj,'out','tables');
idx=readtable(fullfile(T,'session_index.csv'),'TextType','string');
mb =readtable(fullfile(T,'macbp_best.csv'),'TextType','string');

fs=500; F=(15:1:150).'; c1=3; ord=[3 30]; mult=1;
compPre=1000; compPost=3000; padMs=600;
cpre=round(compPre/1000*fs); cpost=round(compPost/1000*fs); pad=round(padMs/1000*fs);
tMsF=(-cpre:cpost)/fs*1000;

fp=idx.finalPath(idx.sessID==string(sessID) & idx.task=="cueTask"); fp=char(fp(1));
D=load(fp); vv=fieldnames(D); od=D.(vv{1}); clear D;
labs=od.labels; bl=char(mb.bestLabel(mb.sessID==string(sessID) & mb.task=="cueTask"));
ci=find(strcmpi(labs,bl),1);
isRsp=cellfun(@(x) ~isempty(regexpi(x,'rsp','once')), labs);
ra=od.data(isRsp,:); ri=1; if isfield(od,'rspIDX'),ri=double(od.rspIDX);end
rf=1; if isfield(od,'rspFlip'),rf=double(od.rspFlip);end
sig=od.data(ci,:); rsp=ra(ri,:)*rf; clear od;

a=onsetSample-cpre-pad; b=onsetSample+cpost+pad;
sw=sig(a:b); rw=rsp(a:b);
P=slt_power_cont(sw,fs,F,c1,ord,mult);
Z=myChanZscore(P.').'; crop=(pad+1):(size(Z,2)-pad); Z=Z(:,crop); Z(~isfinite(Z))=0;
rspF=rw(crop);
[bb,aa]=butter(4,1/(fs/2),'high'); rawHPfull=filtfilt(bb,aa,double(sw)); rawHP=rawHPfull(crop);  % 1 Hz high-pass (broadband)

xl=[-winMs(1) winMs(2)];
dt = tMsF>=xl(1) & tMsF<=xl(2);            % crop plotted traces to the view (no spill past axes)
tD=tMsF(dt); rspD=rspF(dt); rawD=rawHP(dt);
fm = F>=fLim(1) & F<=fLim(2);
sub = Z(fm, dt); climv=prctile(sub(:),climPct); if ~isfinite(climv)||climv<=0, climv=max([sub(:);1]); end
specClim=[0 climv]; if nargin>=7 && ~isempty(climRange), specClim=climRange; end

if strcmpi(mode,'preview')
    fig=figure('Position',[40 40 820 1040],'Color','w','Visible','off');
    draw_stack(fig, tD,rspD,rawD, tMsF,F,Z, xl,fLim,specClim, false);
    annotation('textbox',[0 0.965 1 0.03],'String',sprintf('%s onset=%d win=[-%d %d] f=[%d %d]',sessID,onsetSample,winMs(1),winMs(2),fLim(1),fLim(2)),'EdgeColor','none','HorizontalAlignment','center','FontSize',11,'Interpreter','none');
    exportgraphics(fig,fullfile(scratch,sprintf('preview_%s_%d.png',sessID,onsetSample)),'Resolution',130); close(fig);
    fprintf('wrote preview_%s_%d.png (clim=%.1f)\n',sessID,onsetSample,climv);
else
    outdir=fullfile(proj,'out'); pfx=sprintf('grant_rodentOB_%s_%d',sessID,onsetSample);
    % Panel 1: respiration
    f1=figure('Position',[40 40 1100 320],'Color','w','Visible','off');
    plot(tD,rspD,'k','LineWidth',4); xlim(xl); box off; set(gca,'FontSize',30,'LineWidth',3.5,'XTickLabel',[],'YTick',[]);
    ylabel('Resp','FontWeight','bold','FontSize',42); xline(0,'k--','LineWidth',3);
    exportgraphics(f1,fullfile(outdir,[pfx '_1resp.png']),'Resolution',200); close(f1);
    % Panel 2: OB ephys raw (>20 Hz)
    f2=figure('Position',[40 40 1250 320],'Color','w','Visible','off');
    plot(tD,rawD,'k','LineWidth',2); xlim(xl); box off; set(gca,'FontSize',30,'LineWidth',3.5,'XTickLabel',[],'YTick',[],'Position',[0.14 0.17 0.72 0.75]);
    mxr=max(abs(rawD)); ylim([-1.6*mxr 1.6*mxr]);   % wider Y range with headroom
    ylabel('OB ephys (\muV)','FontWeight','bold','FontSize',42); xline(0,'k--','LineWidth',3);
    add_scalebar(gca, xl, 50);
    exportgraphics(f2,fullfile(outdir,[pfx '_2bulbRaw.png']),'Resolution',200); close(f2);
    % Panel 3: OB spectrogram
    f3=figure('Position',[40 40 1180 360],'Color','w','Visible','off');
    imagesc(tMsF,F,Z,specClim); axis xy; colormap(jet); xlim(xl); ylim(fLim);
    set(gca,'FontSize',30,'LineWidth',3.5); ylabel('Hz','FontWeight','bold','FontSize',42);
    xlabel('Time from inhale (ms)','FontWeight','bold','FontSize',38); xline(0,'w--','LineWidth',3);
    cb=colorbar; cb.LineWidth=3; cb.FontSize=26; cb.Label.String='power (z)'; cb.Label.FontSize=30;
    exportgraphics(f3,fullfile(outdir,[pfx '_3bulbSpec.png']),'Resolution',200); close(f3);
    % Panel 4: stacked, x-aligned
    f4=figure('Position',[40 40 1150 1150],'Color','w','Visible','off');
    draw_stack(f4, tD,rspD,rawD, tMsF,F,Z, xl,fLim,specClim, true);
    exportgraphics(f4,fullfile(outdir,[pfx '_4stack.png']),'Resolution',200); close(f4);
    fprintf('wrote %s_{1resp,2bulbRaw,3bulbSpec,4stack}.png to out/ (clim=%.1f)\n',pfx,climv);
end
end

% ---- stacked 3-panel with manually aligned x-axes (shared left/width; colorbar manual) ----
function draw_stack(fig, tD,rspD,rawD, tMsF,F,Z, xl,fLim,clim2, big)
if big, fsA=30; fsL=40; lw=3.5; lwt=3; else, fsA=15; fsL=18; lw=2; lwt=2; end
left=0.175; w=0.65; h=0.23; yb=0.125; gap=0.06;
pS=[left yb w h]; pR=[left yb+h+gap w h]; pP=[left yb+2*(h+gap) w h];
% resp (top)
axR=axes(fig,'Position',pP); plot(axR,tD,rspD,'k','LineWidth',lwt); xlim(axR,xl); box(axR,'off');
set(axR,'FontSize',fsA,'LineWidth',lw,'XTickLabel',[],'YTick',[]); ylabel(axR,'Resp','FontWeight','bold','FontSize',fsL);
xline(axR,0,'k--','LineWidth',lw);
% raw (middle)
axW=axes(fig,'Position',pR); plot(axW,tD,rawD,'k','LineWidth',max(1.5,lwt-1)); xlim(axW,xl); box(axW,'off');
mxr=max(abs(rawD)); ylim(axW,[-1.6*mxr 1.6*mxr]);   % wider Y range with headroom
set(axW,'FontSize',fsA,'LineWidth',lw,'XTickLabel',[],'YTick',[]); ylabel(axW,'OB ephys (\muV)','FontWeight','bold','FontSize',fsL);
xline(axW,0,'k--','LineWidth',lw); add_scalebar(axW, xl, 50);
% spectrogram (bottom)
axS=axes(fig,'Position',pS); imagesc(axS,tMsF,F,Z,clim2); axis(axS,'xy'); colormap(axS,jet);
xlim(axS,xl); ylim(axS,fLim); set(axS,'FontSize',fsA,'LineWidth',lw); ylabel(axS,'Hz','FontWeight','bold','FontSize',fsL);
xlabel(axS,'Time from inhale (ms)','FontWeight','bold','FontSize',fsL); xline(axS,0,'w--','LineWidth',lw);
cb=colorbar(axS); axS.Position=pS; cb.Position=[left+w+0.015 yb 0.022 h];
cb.LineWidth=lw; cb.FontSize=max(14,fsA-2); cb.Label.String='power (z)';
end

% ---- amplitude scale bar just outside the right edge (per-object clipping off; axes clipping stays on) ----
function add_scalebar(ax, xl, sb)
x0=xl(2)+0.02*diff(xl); y0=-sb/2;
l=line(ax,[x0 x0],[y0 y0+sb],'Color','k','LineWidth',6); set(l,'Clipping','off');
t=text(ax,x0+0.014*diff(xl),y0+sb/2,sprintf('%d\\muV',sb),'FontSize',26,'FontWeight','bold','HorizontalAlignment','left');
set(t,'Clipping','off');
end
