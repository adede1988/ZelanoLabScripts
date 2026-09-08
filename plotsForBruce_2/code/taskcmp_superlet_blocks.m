function taskcmp_superlet_blocks()
% CANONICAL pipeline (matches the spectro2 / peak-latency analysis), control only,
% baseline (audiobook) vs ATB (focused breathing):
%   superlet power (slt_power_cont, c1=3, ord[2 4], F=25:1:58) on best-macBP channel
%   -> within-frequency z to the pooled per-CONDITION baseline (-500..-100 ms):
%        z[f,t] = (pow[f,t] - baseMu[f]) / baseStd[f], baseMu/baseStd pooled over that
%        condition's breaths x baseline samples (myChanZscore convention, baseStd form)
%   -> collapse band: bz[t] = mean over the 34 freqs.
% Per breath = bz (not averaged); block trajectory = mean of bz over that condition's breaths.
% Writes two figures + caches all traces to taskcmp_superlet_blocks.mat. CP excluded.
proj='C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2';
addpath(fullfile(proj,'code')); addpath('C:/Users/Adam/Documents/GitHub/ZelanoLabScripts');
addpath('C:/Users/Adam/Documents/GitHub/Superlets/matlab-pure');
tdir=fullfile(proj,'out','tables'); fdir=fullfile(proj,'out','figs','taskcmp');
F=25:1:58; C1=3; ORD=[2 4]; MULT=1; fs=500; EXCL={'CP'}; CLIP=15;
GREY=[0.62 0.62 0.62]; GREEN=[0.133 0.545 0.133];
e0=round(-0.5*fs); e1=round(2.0*fs); tMs=(e0:e1)/fs*1000; nT=numel(tMs); preM=tMs>=-500&tMs<=-100;
best=readtable(fullfile(tdir,'macbp_best.csv'),'TextType','string');
best=best(best.cohort=="OBE" & best.task=="breathingTask",:);
base='R:/Neurology/Zelano_Lab/Lab_Common/OBEControl/';
isAudio=@(s) strcmpi(s,'audio')|strcmpi(s,'audiobook');
isFocus=@(s) strcmpi(s,'focus')|strcmpi(s,'naturalFocus')|strcmpi(s,'slowFocus')|strcmpi(s,'focusedBreathing');
conds={'baseline','ATB'};
PB=struct('cond',{{}},'sess',{{}},'bz',{{}});   % per-breath band-z
BLK=struct('cond',{{}},'sess',{{}},'bz',{{}});   % per-session-condition block mean
for r=1:height(best)
  id=char(best.sessID(r)); part=char(best.participant(r)); bl=char(best.bestLabel(r));
  if any(strcmpi(part,EXCL)), continue; end
  dd=dir(fullfile(base,id,'preProc','*breathing*.mat')); if isempty(dd), continue; end
  od=[]; for a=1:3, try, S=load(fullfile(dd(1).folder,dd(1).name)); f=fieldnames(S); od=S.(f{1}); clear S; break; catch e, pause(3); end, end
  if isempty(od), fprintf('%s LOAD FAIL\n',id); continue; end
  labs=cellfun(@(x)char(string(x)),od.labels,'uni',0); ci=find(strcmpi(labs,bl),1); if isempty(ci), continue; end
  x=double(od.data(ci,:)); x=fillmissing(x,'linear'); x=fillmissing(x,'nearest'); N=numel(x);
  P=slt_power_cont(x, fs, F, C1, ORD, MULT);   % [nF x N] superlet power
  bd=od.behDat; io=round(double(bd.finalOnset)); tc=strtrim(string(getcol(bd,'task')));
  for cix=1:2
    if cix==1, sel=isAudio(tc); else, sel=isFocus(tc); end
    ks=find(sel(:).'); if isempty(ks), continue; end
    % gather this condition's epochs [nF x nT x nBreath]
    epo=[]; keep=[];
    for k=ks
      s0=io(k)+e0; s1=io(k)+e1; if isnan(io(k))||s0<1||s1>N, continue; end
      epo=cat(3, epo, P(:, s0:s1)); keep(end+1)=k; %#ok<AGROW>
    end
    if isempty(epo), continue; end
    nb=size(epo,3);
    baseAll=reshape(epo(:,preM,:), numel(F), []);   % [nF x (nBase*nBreath)]
    baseMu=mean(baseAll,2); baseStd=std(baseAll,0,2); baseStd(baseStd<=0)=eps;
    bzMat=nan(nb,nT);
    for b=1:nb
      zf=(epo(:,:,b)-baseMu)./baseStd;              % within-freq z to pooled baseline
      bz=mean(zf,1); bz=max(min(bz,CLIP),-CLIP);     % collapse 34 freqs
      bzMat(b,:)=bz;
      PB.cond{end+1}=conds{cix}; PB.sess{end+1}=part; PB.bz{end+1}=bz;
    end
    BLK.cond{end+1}=conds{cix}; BLK.sess{end+1}=part; BLK.bz{end+1}=mean(bzMat,1,'omitnan');
    fprintf('%s (%s) %s: %d breaths\n', id, part, conds{cix}, nb);
  end
end
save(fullfile(tdir,'taskcmp_superlet_blocks.mat'),'PB','BLK','tMs','F','conds');

show=tMs>=-200 & tMs<=1500; tt=tMs(show);
% ---- Figure 1: block (session-mean) trajectories ----
fig=figure('Position',[30 30 1200 700],'Color','w','Visible','off'); hold on;
gm=@(c) sess_mean(BLK,c,show);
for cix=1:2, col=ternIf(cix==1,GREY,GREEN); idx=find(strcmp(BLK.cond,conds{cix}));
  for i=idx, plot(tt, BLK.bz{i}(show),'-','Color',[col 0.5],'LineWidth',1.3); end, end
[mB,~]=gm('baseline'); [mA,~]=gm('ATB');
hB=plot(tt,mB,'-','Color',GREY*0.5,'LineWidth',3.5); hA=plot(tt,mA,'-','Color',GREEN,'LineWidth',3.5);
xline(0,'k-','LineWidth',1); yline(0,'k:');
xlabel('time from inhale onset (ms)','FontSize',15,'FontWeight','bold');
ylabel('gamma band-z (superlet, 25-58 Hz, myChanZscore)','FontSize',15,'FontWeight','bold');
nB=sum(strcmp(BLK.cond,'baseline')); nA=sum(strcmp(BLK.cond,'ATB'));
legend([hB hA],{sprintf('baseline (%d sessions)',nB),sprintf('ATB (%d sessions)',nA)},'FontSize',14,'Location','northeast');
title('Control: per-session BLOCK gamma trajectories (canonical superlet pipeline)','FontSize',15,'FontWeight','bold');
set(gca,'FontSize',13,'LineWidth',1.2); xlim([-200 1500]);
exportgraphics(fig, fullfile(fdir,'superlet_block_trajectories.png'),'Resolution',150); close(fig);
% ---- Figure 2: per-breath traces ----
fig=figure('Position',[30 30 1250 720],'Color','w','Visible','off'); hold on;
iB=find(strcmp(PB.cond,'baseline')); iA=find(strcmp(PB.cond,'ATB'));
aB=min(0.05,8/max(numel(iB),1)); aA=min(0.05,8/max(numel(iA),1));
for i=iB, plot(tt, PB.bz{i}(show),'-','Color',[GREY aB],'LineWidth',0.4); end
for i=iA, plot(tt, PB.bz{i}(show),'-','Color',[GREEN aA],'LineWidth',0.4); end
MB=mean(cell2mat(PB.bz(iB)'),1,'omitnan'); MA=mean(cell2mat(PB.bz(iA)'),1,'omitnan');
hB=plot(tt,MB(show),'-','Color',GREY*0.5,'LineWidth',3); hA=plot(tt,MA(show),'-','Color',GREEN,'LineWidth',3);
xline(0,'k-','LineWidth',1); yline(0,'k:');
xlabel('time from inhale onset (ms)','FontSize',15,'FontWeight','bold');
ylabel('gamma band-z (superlet, per breath)','FontSize',15,'FontWeight','bold');
legend([hB hA],{sprintf('baseline (%d breaths)',numel(iB)),sprintf('ATB (%d breaths)',numel(iA))},'FontSize',14,'Location','northeast');
title('Control: per-breath gamma band-z (same pipeline, not averaged) + bold means','FontSize',15,'FontWeight','bold');
set(gca,'FontSize',13,'LineWidth',1.2); xlim([-200 1500]); ylim([-6 8]);
exportgraphics(fig, fullfile(fdir,'superlet_perbreath.png'),'Resolution',150); close(fig);
fprintf('wrote superlet_block_trajectories.png, superlet_perbreath.png (+ .mat)\n');
end
function [m,n]=sess_mean(BLK,c,show), idx=find(strcmp(BLK.cond,c)); M=cell2mat(BLK.bz(idx)'); m=mean(M(:,show),1,'omitnan'); n=numel(idx); end
function y=ternIf(c,a,b), if c, y=a; else, y=b; end, end
function s=getcol(bd,nm), n=height(bd); s=repmat({''},n,1); if ~any(strcmp(bd.Properties.VariableNames,nm)), return; end
  v=bd.(nm); for i=1:n, try, if iscell(v), e=v{i}; else, e=v(i); end, t=string(e); if ~isscalar(t), if isempty(t), t=""; else, t=t(1); end, end, s{i}=char(t); catch, end, end, end
