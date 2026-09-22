function taskcmp_olf_maps()
% Olfactory vs non-olfactory gamma, CONTROL (OBE) only, CP dropped.
% ONE decomposition for everything (canonical superlet c1=3, order [3 30], F=22:62);
% z-score = myChanZscore (per-freq, baseline -500..-100 ms, denom baseStd/sqrt(nBreaths),
% winsorize +/-10) -> across-breath z-MAP per task block.
%
% Blocks: olfactory = cueTask + threshTask + O15 (all sniffs, finalOnset-locked);
%         non-olfactory = audiobook + focusedBreathing (breath blocks, inhale-onset-locked).
% Outputs:
%   out/tables/taskcmp_olf_blocks.mat   (per-block z-maps + meta)
%   out/tables/taskcmp_olf_roi.csv      (long: block x gammaType median-z, for the LMM)
%   out/figs/taskcmp/olf_vs_nonolf_spectrograms.png  (olfactory / non-olf / difference)
proj='C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2';
addpath(fullfile(proj,'code')); addpath('C:/Users/Adam/Documents/GitHub/ZelanoLabScripts');
addpath('C:/Users/Adam/Documents/GitHub/Superlets/matlab-pure');
tdir=fullfile(proj,'out','tables'); fdir=fullfile(proj,'out','figs','taskcmp');
F=22:1:62; C1=3; ORD=[3 30]; MULT=1; fs=500; EXCL={'CP'}; CLIPZ=10;
e0=round(-1.0*fs); e1=round(4.0*fs); tMs=(e0:e1)/fs*1000; nT=numel(tMs);
bMask=tMs>=-500 & tMs<=-100; nBaseSamp=nnz(bMask);
% ROIs
fmHi=F>=40 & F<=58;  tmEarly=tMs>=0 & tMs<=1000;      % early high gamma
fmLo=F>=25 & F<=40;  tmLate =tMs>=1000 & tMs<=2000;   % late low gamma
% storage/display window
showMask=tMs>=-500 & tMs<=3000; tShow=tMs(showMask);

best=readtable(fullfile(tdir,'macbp_best.csv'),'TextType','string');
best=best(best.cohort=="OBE",:);
base='R:/Neurology/Zelano_Lab/Lab_Common/OBEControl/';
isAudio=@(s) strcmpi(s,'audio')|strcmpi(s,'audiobook');
isFocus=@(s) strcmpi(s,'focus')|strcmpi(s,'naturalFocus')|strcmpi(s,'slowFocus')|strcmpi(s,'focusedBreathing');
glob=containers.Map({'cueTask','threshTask','O15','breathingTask'}, ...
                    {'*cueTask*.mat','*threshold*.mat','*O15*.mat','*breathing*.mat'});

B=struct('sessID',{},'participant',{},'task',{},'category',{},'block',{}, ...
         'nBreaths',{},'z',{},'earlyHigh',{},'lateLow',{});
% medZ = myChanZscore ROI median (SEM denom, ~scales with sqrt(nBreaths));
% medZ_ninv = medZ/sqrt(nBreaths) = per-breath-scale z, trial-count-invariant (robustness).
ROI=table('Size',[0 8],'VariableTypes',{'string','string','string','string','string','double','double','double'}, ...
          'VariableNames',{'sessID','participant','task','category','gammaType','medZ','medZ_ninv','nBreaths'});

for r=1:height(best)
  id=char(best.sessID(r)); part=char(best.participant(r)); bl=char(best.bestLabel(r)); tk=char(best.task(r));
  if any(strcmpi(part,EXCL)), continue; end
  if ~isKey(glob,tk), continue; end
  dd=dir(fullfile(base,id,'preProc',glob(tk)));
  if isempty(dd), fprintf('%s %s: no final\n',id,tk); continue; end
  fp=pick_final(dd);
  od=[]; for a=1:6, try, S=load(fp); f=fieldnames(S); od=S.(f{1}); clear S; break; catch, fprintf('[retry %d] ',a); pause(6*a); end, end
  if isempty(od), fprintf('%s %s LOAD FAIL\n',id,tk); continue; end
  labs=cellfun(@(x)char(string(x)),od.labels,'uni',0); ci=find(strcmpi(labs,bl),1);
  if isempty(ci), fprintf('%s %s: label %s missing\n',id,tk,bl); continue; end
  x=double(od.data(ci,:)); x=fillmissing(x,'linear'); x=fillmissing(x,'nearest'); N=numel(x);
  P=slt_power_cont(x, fs, F, C1, ORD, MULT);   % [nF x N] raw power, canonical [3 30]
  bd=od.behDat; fo=round(coerce(bd.finalOnset));

  % define blocks for this final
  if strcmp(tk,'breathingTask')
    tc=strtrim(string(getcol(bd,'task')));
    blocks={'audiobook', fo(isAudio(tc)); 'focusedBreathing', fo(isFocus(tc))};
    catg='nonolf';
  else
    blocks={tk, fo};   % all sniffs
    catg='olf';
  end

  for bi=1:size(blocks,1)
    bname=blocks{bi,1}; io=blocks{bi,2}; io=io(~isnan(io)&io>0);
    [z,nB]=zmap(P,io,e0,e1,bMask,nBaseSamp);   % UN-winsorized myChanZscore z-map
    if nB<5, fprintf('%s %s: thin (%d breaths)\n',id,bname,nB); continue; end
    eh=median(z(fmHi,tmEarly),'all'); ll=median(z(fmLo,tmLate),'all'); rn=sqrt(nB);
    B(end+1)=struct('sessID',id,'participant',part,'task',tk,'category',catg,'block',bname, ...
                    'nBreaths',nB,'z',z(:,showMask),'earlyHigh',eh,'lateLow',ll); %#ok<AGROW>
    ROI=[ROI; table(string(id),string(part),string(bname),string(catg),"earlyHigh",eh,eh/rn,nB, ...
                    'VariableNames',ROI.Properties.VariableNames); ...
             table(string(id),string(part),string(bname),string(catg),"lateLow", ll,ll/rn,nB, ...
                    'VariableNames',ROI.Properties.VariableNames)]; %#ok<AGROW>
    fprintf('%s | %-16s | %-6s | n=%3d | earlyHigh=%+.2f lateLow=%+.2f\n', id, bname, catg, nB, eh, ll);
  end
end

save(fullfile(tdir,'taskcmp_olf_blocks.mat'),'B','tShow','F','-v7.3');
writetable(ROI, fullfile(tdir,'taskcmp_olf_roi.csv'));
fprintf('\n%d blocks (%d olf, %d nonolf); %d distinct sessions\n', numel(B), ...
    sum(strcmp({B.category},'olf')), sum(strcmp({B.category},'nonolf')), numel(unique({B.sessID})));

% ---- averaged spectrograms (mean z-map over blocks within category) ----
olf = cat_mean(B,'olf'); non = cat_mean(B,'nonolf'); dif = olf-non;
fShow=F; tt=tShow;
fig=figure('Position',[20 20 1650 480],'Color','w','Visible','off');
climAbs=max([abs(prctile(olf(:),[2 98])) abs(prctile(non(:),[2 98]))]);
titles={sprintf('Olfactory (cue+thresh+O15), %d blocks',sum(strcmp({B.category},'olf'))), ...
        sprintf('Non-olfactory (audiobook+focus), %d blocks',sum(strcmp({B.category},'nonolf'))), ...
        'Difference (olfactory - non-olfactory)'};
maps={olf,non,dif}; clims={[-climAbs climAbs],[-climAbs climAbs],[]};
for p=1:3
  subplot(1,3,p);
  imagesc(tt, fShow, maps{p}); axis xy; hold on;
  yl=[25 58]; ylim(yl); xlim([-500 2500]);
  if p<3, clim(clims{p}); else, dc=max(abs(prctile(dif(:),[2 98]))); clim([-dc dc]); end
  colormap(gca, redblue()); colorbar;
  % ROI boxes
  rectangle('Position',[0 40 1000 18],'EdgeColor','k','LineWidth',1.6,'LineStyle','-');       % early high
  rectangle('Position',[1000 25 1000 15],'EdgeColor','k','LineWidth',1.6,'LineStyle','--');    % late low
  xline(0,'k-','LineWidth',1); yline(40,'k:');
  xlabel('time from onset (ms)','FontSize',12,'FontWeight','bold'); ylabel('Hz','FontSize',12,'FontWeight','bold');
  title(titles{p},'FontSize',12,'FontWeight','bold'); set(gca,'FontSize',11);
end
sgtitle('Control: olfactory vs non-olfactory gamma z (myChanZscore, superlet [3 30]); solid box=early high, dashed=late low','FontSize',13,'FontWeight','bold');
exportgraphics(fig, fullfile(fdir,'olf_vs_nonolf_spectrograms.png'),'Resolution',150); close(fig);
fprintf('wrote taskcmp_olf_blocks.mat, taskcmp_olf_roi.csv, olf_vs_nonolf_spectrograms.png\n');
end

% ---- helpers ----
function [z,nB]=zmap(P,io,e0,e1,bMask,nBaseSamp)
nF=size(P,1); N=size(P,2); nT=e1-e0+1;
sumPow=zeros(nF,nT); baseSum=zeros(nF,1); baseSumSq=zeros(nF,1); nB=0;
for k=1:numel(io)
  a=io(k)+e0; b=io(k)+e1; if a<1||b>N, continue; end
  Pe=P(:,a:b); sumPow=sumPow+Pe; Pb=Pe(:,bMask);
  baseSum=baseSum+sum(Pb,2); baseSumSq=baseSumSq+sum(Pb.^2,2); nB=nB+1;
end
if nB<1, z=zeros(nF,nT); return; end
nb=nB*nBaseSamp; baseMu=baseSum/nb; baseVar=baseSumSq/nb-baseMu.^2; baseVar(baseVar<=0)=eps;
z=(sumPow/nB - baseMu)./sqrt(baseVar/nB);          % myChanZscore across-breath z-map (UN-winsorized)
end
function M=cat_mean(B,catg)
% winsorize each block's z at +/-10 (display robustness, per spectro2) THEN average
idx=find(strcmp({B.category},catg)); Z=cat(3,B(idx).z); Z=min(max(Z,-10),10); M=mean(Z,3,'omitnan');
end
function fp=pick_final(dd)
% prefer a file whose name contains 'preproc' (the final), else the first
names=lower({dd.name}); k=find(contains(names,'preproc'),1); if isempty(k),k=1; end
fp=fullfile(dd(k).folder,dd(k).name);
end
function v=coerce(c)
if isnumeric(c), v=double(c(:)); return; end
if iscell(c), v=nan(numel(c),1); for i=1:numel(c), try, e=c{i}; if isnumeric(e)&&~isempty(e), v(i)=double(e(1)); else, v(i)=double(string(e)); end, catch, v(i)=NaN; end, end, return; end
try, v=double(string(c)); v=v(:); catch, v=nan(numel(c),1); end
end
function s=getcol(bd,nm), n=height(bd); s=repmat({''},n,1); if ~any(strcmp(bd.Properties.VariableNames,nm)), return; end
  v=bd.(nm); for i=1:n, try, if iscell(v), e=v{i}; else, e=v(i); end, t=string(e); if ~isscalar(t), if isempty(t), t=""; else, t=t(1); end, end, s{i}=char(t); catch, end, end, end
function c=redblue()
n=256; c=interp1([0;0.5;1], [0.23 0.30 0.75; 1 1 1; 0.75 0.20 0.20], linspace(0,1,n)');
end
