function taskcmp_ridge_diff()
% Ridge-FREQUENCY sweep, canonical pipeline (extract_gamma_session): superlet
% (slt_power_cont, c1=3, ord [3 30], F=22:62) on best-macBP channel -> per breath,
% per-freq WHOLE-WINDOW z over the epoch (Zw), ridge_track(max(Zw,0)) on 25-58 Hz
% -> ridge frequency over time. Control only (CP excluded), baseline vs ATB.
% Per session: mean ridge-freq sweep per condition; within-session difference
% (ATB - baseline). Plots the two sweeps + the session-wise difference.
proj='C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2';
addpath(fullfile(proj,'code')); addpath('C:/Users/Adam/Documents/GitHub/ZelanoLabScripts');
addpath('C:/Users/Adam/Documents/GitHub/Superlets/matlab-pure');
tdir=fullfile(proj,'out','tables'); fdir=fullfile(proj,'out','figs','taskcmp');
F=22:1:62; C1=3; ORD=[3 30]; MULT=1; RB=[25 58]; PEN=1.0; RBW=2; fs=500; EXCL={'CP'};
bIdx=find(F>=RB(1)&F<=RB(2)); Fb=F(bIdx)';
GREY=[0.62 0.62 0.62]; GREEN=[0.133 0.545 0.133];
e0=round(-0.5*fs); e1=round(1.5*fs); tMs=(e0:e1)/fs*1000; nT=numel(tMs);
best=readtable(fullfile(tdir,'macbp_best.csv'),'TextType','string');
best=best(best.cohort=="OBE" & best.task=="breathingTask",:);
base='R:/Neurology/Zelano_Lab/Lab_Common/OBEControl/';
isAudio=@(s) strcmpi(s,'audio')|strcmpi(s,'audiobook');
isFocus=@(s) strcmpi(s,'focus')|strcmpi(s,'naturalFocus')|strcmpi(s,'slowFocus')|strcmpi(s,'focusedBreathing');
SB=[]; SF=[]; sess={};   % per-session mean ridge-freq sweeps
for r=1:height(best)
  id=char(best.sessID(r)); part=char(best.participant(r)); bl=char(best.bestLabel(r));
  if any(strcmpi(part,EXCL)), continue; end
  dd=dir(fullfile(base,id,'preProc','*breathing*.mat')); if isempty(dd), continue; end
  od=[]; for a=1:3, try, S=load(fullfile(dd(1).folder,dd(1).name)); f=fieldnames(S); od=S.(f{1}); clear S; break; catch e, pause(3); end, end
  if isempty(od), fprintf('%s LOAD FAIL\n',id); continue; end
  labs=cellfun(@(x)char(string(x)),od.labels,'uni',0); ci=find(strcmpi(labs,bl),1); if isempty(ci), continue; end
  x=double(od.data(ci,:)); x=fillmissing(x,'linear'); x=fillmissing(x,'nearest'); N=numel(x);
  P=slt_power_cont(x, fs, F, C1, ORD, MULT);   % [nF x N]
  bd=od.behDat; io=round(double(bd.finalOnset)); tc=strtrim(string(getcol(bd,'task')));
  RA=[]; RF=[];
  for k=1:numel(io)
    s0=io(k)+e0; s1=io(k)+e1; if isnan(io(k))||s0<1||s1>N, continue; end
    Pb=P(bIdx, s0:s1); Zw=(Pb-mean(Pb,2))./max(std(Pb,0,2),eps);
    fr=ridge_track(max(Zw,0), Fb, PEN, 1, RBW); fr=fr(:).';    % ridge freq (Hz) over time
    if isAudio(tc(k)), RA=[RA; fr]; elseif isFocus(tc(k)), RF=[RF; fr]; end %#ok<AGROW>
  end
  if isempty(RA)||isempty(RF), continue; end
  SB=[SB; mean(RA,1,'omitnan')]; SF=[SF; mean(RF,1,'omitnan')]; sess{end+1}=part; %#ok<AGROW>
  fprintf('%s (%s): %d baseline, %d ATB breaths\n', id, part, size(RA,1), size(RF,1));
end
save(fullfile(tdir,'taskcmp_ridge_diff.mat'),'SB','SF','sess','tMs');
D=SF-SB; nS=size(D,1);
show=tMs>=-100 & tMs<=1500; tt=tMs(show);

fig=figure('Position',[20 20 1400 620],'Color','w','Visible','off');
% left: the two sweeps
subplot(1,2,1); hold on;
for i=1:nS, plot(tt, SB(i,show),'-','Color',[GREY 0.45],'LineWidth',1); plot(tt, SF(i,show),'-','Color',[GREEN 0.45],'LineWidth',1); end
hB=plot(tt, mean(SB(:,show),1),'-','Color',GREY*0.5,'LineWidth',3.5);
hA=plot(tt, mean(SF(:,show),1),'-','Color',GREEN,'LineWidth',3.5);
xline(0,'k-'); xlabel('time from inhale onset (ms)','FontSize',14,'FontWeight','bold');
ylabel('ridge frequency (Hz)','FontSize',14,'FontWeight','bold');
legend([hB hA],{'baseline','ATB'},'FontSize',13,'Location','best');
title('Ridge-frequency sweep by session','FontSize',14,'FontWeight','bold'); set(gca,'FontSize',12,'LineWidth',1.1);
% right: session-wise difference
subplot(1,2,2); hold on;
for i=1:nS, plot(tt, D(i,show),'-','Color',[0.3 0.3 0.3 0.4],'LineWidth',1); end
mD=mean(D(:,show),1); sD=std(D(:,show),0,1)/sqrt(nS);
fill([tt fliplr(tt)],[mD+sD fliplr(mD-sD)],[0.2 0.5 0.2],'FaceAlpha',0.2,'EdgeColor','none');
plot(tt, mD,'-','Color',[0.13 0.45 0.13],'LineWidth',3.5); yline(0,'k:','LineWidth',1); xline(0,'k-');
xlabel('time from inhale onset (ms)','FontSize',14,'FontWeight','bold');
ylabel('\Delta ridge frequency (ATB - baseline, Hz)','FontSize',14,'FontWeight','bold');
title(sprintf('Session-wise sweep difference (n=%d, mean\\pmSEM)',nS),'FontSize',14,'FontWeight','bold'); set(gca,'FontSize',12,'LineWidth',1.1);
sgtitle('Control: gamma ridge-frequency sweep, baseline vs ATB (canonical ridge pipeline)','FontSize',15,'FontWeight','bold');
exportgraphics(fig, fullfile(fdir,'ridge_sweep_diff.png'),'Resolution',150); close(fig);
fprintf('wrote ridge_sweep_diff.png (%d sessions)\n', nS);
end
function s=getcol(bd,nm), n=height(bd); s=repmat({''},n,1); if ~any(strcmp(bd.Properties.VariableNames,nm)), return; end
  v=bd.(nm); for i=1:n, try, if iscell(v), e=v{i}; else, e=v(i); end, t=string(e); if ~isscalar(t), if isempty(t), t=""; else, t=t(1); end, end, s{i}=char(t); catch, end, end, end
