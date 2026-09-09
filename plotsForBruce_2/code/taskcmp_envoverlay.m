function taskcmp_envoverlay()
% Per-SESSION mean gamma (25-58 Hz) Hilbert-envelope traces, aligned to inhale
% onset, baseline-normalised (-500..-100 ms), overlaid: baseline (audiobook) grey,
% ATB (focused breathing) green; bold = grand mean across sessions. One line per
% control session (CP excluded). Best-macBP channel + envelope as gammaHRV.
proj='C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2';
tdir=fullfile(proj,'out','tables'); fdir=fullfile(proj,'out','figs','taskcmp');
GAMMA=[25 58]; EXCL={'CP'}; CLIP=8; fs=500;
GREY=[0.62 0.62 0.62]; GREEN=[0.133 0.545 0.133];
e0=round(-0.5*fs); e1=round(2.0*fs); tMs=(e0:e1)/fs*1000; nT=numel(tMs);
preM=tMs>=-500 & tMs<=-100;
best=readtable(fullfile(tdir,'macbp_best.csv'),'TextType','string');
best=best(best.cohort=="OBE" & best.task=="breathingTask",:);
base='R:/Neurology/Zelano_Lab/Lab_Common/OBEControl/';
isAudio=@(s) strcmpi(s,'audio')|strcmpi(s,'audiobook');
isFocus=@(s) strcmpi(s,'focus')|strcmpi(s,'naturalFocus')|strcmpi(s,'slowFocus')|strcmpi(s,'focusedBreathing');
sm=@(y) movmean(y,round(0.05*fs));   % light 50 ms smooth for display
SA=[]; SF=[]; sess={}; nAll=[0 0];
for r=1:height(best)
  id=char(best.sessID(r)); part=char(best.participant(r)); bl=char(best.bestLabel(r));
  if any(strcmpi(part,EXCL)), continue; end
  dd=dir(fullfile(base,id,'preProc','*breathing*.mat')); if isempty(dd), continue; end
  od=[]; for a=1:3, try, S=load(fullfile(dd(1).folder,dd(1).name)); f=fieldnames(S); od=S.(f{1}); clear S; break; catch e, pause(3); end, end
  if isempty(od), fprintf('%s LOAD FAIL\n',id); continue; end
  labs=cellfun(@(x)char(string(x)),od.labels,'uni',0); ci=find(strcmpi(labs,bl),1); if isempty(ci), continue; end
  x=double(od.data(ci,:)); x=fillmissing(x,'linear'); x=fillmissing(x,'nearest'); N=numel(x);
  [bB,aB]=butter(4,GAMMA/(fs/2),'bandpass'); env=abs(hilbert(filtfilt(bB,aB,x)));
  bd=od.behDat; io=round(double(bd.finalOnset)); tc=strtrim(string(getcol(bd,'task')));
  Ai=[]; Fi=[];
  for k=1:numel(io)
    s0=io(k)+e0; s1=io(k)+e1; if isnan(io(k))||s0<1||s1>N, continue; end
    seg=env(s0:s1); bmu=mean(seg(preM)); bsd=std(seg(preM)); if bsd<=0, bsd=eps; end
    z=max(min((seg-bmu)/bsd,CLIP),-CLIP);
    if isAudio(tc(k)), Ai=[Ai; z]; elseif isFocus(tc(k)), Fi=[Fi; z]; end %#ok<AGROW>
  end
  if isempty(Ai)||isempty(Fi), continue; end
  SA=[SA; sm(mean(Ai,1))]; SF=[SF; sm(mean(Fi,1))]; sess{end+1}=part; %#ok<AGROW>
  nAll=nAll+[size(Ai,1) size(Fi,1)];
  fprintf('%s (%s): %d baseline, %d ATB breaths\n', id, part, size(Ai,1), size(Fi,1));
end
save(fullfile(tdir,'taskcmp_envoverlay.mat'),'SA','SF','sess','tMs','nAll');
fprintf('%d sessions; %d baseline / %d ATB breaths total\n', size(SA,1), nAll(1), nAll(2));

show=tMs>=-200 & tMs<=1800; tt=tMs(show);
fig=figure('Position',[30 30 1250 720],'Color','w','Visible','off'); hold on;
for i=1:size(SA,1), plot(tt, SA(i,show), '-', 'Color',[GREY 0.5], 'LineWidth',1.3); end
for i=1:size(SF,1), plot(tt, SF(i,show), '-', 'Color',[GREEN 0.5], 'LineWidth',1.3); end
hB=plot(tt, mean(SA(:,show),1), '-', 'Color',GREY*0.5, 'LineWidth',3.5);
hA=plot(tt, mean(SF(:,show),1), '-', 'Color',GREEN, 'LineWidth',3.5);
xline(0,'k-','LineWidth',1); yline(0,'k:');
xlabel('time from inhale onset (ms)','FontSize',15,'FontWeight','bold');
ylabel('gamma (25-58 Hz) envelope, baseline-z','FontSize',15,'FontWeight','bold');
legend([hB hA], {sprintf('baseline (%d sessions)',size(SA,1)), sprintf('ATB (%d sessions)',size(SF,1))}, ...
  'FontSize',14,'Location','northeast');
title('Control: per-session mean gamma envelopes (baseline vs ATB) + bold grand mean','FontSize',15,'FontWeight','bold');
set(gca,'FontSize',13,'LineWidth',1.2); xlim([-200 1800]);
exportgraphics(fig, fullfile(fdir,'env_overlay_sessions.png'),'Resolution',150); close(fig);
fprintf('wrote env_overlay_sessions.png\n');
end
function s=getcol(bd,nm), n=height(bd); s=repmat({''},n,1); if ~any(strcmp(bd.Properties.VariableNames,nm)), return; end
  v=bd.(nm); for i=1:n, try, if iscell(v), e=v{i}; else, e=v(i); end, t=string(e); if ~isscalar(t), if isempty(t), t=""; else, t=t(1); end, end, s{i}=char(t); catch, end, end, end
