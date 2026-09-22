function taskcmp_olf_phase(testN)
% Respiratory-PHASE-matched olfactory vs non-olfactory gamma, control (OBE), CP dropped.
% 5 phases per event: (1) onset->inhale peak, (2) peak->exhale onset (inhale pause folded in),
% (3) exhale onset->trough, (4) trough->exhale-pause onset, (5) pause->breath end.
% Breathing (audiobook+focus combined): landmarks from bm_* columns (samples).
% Olfactory (cue/thresh/O15): landmarks derived per sniff from the respiration flow (SAME
%   detector is also run on breathing to VALIDATE it against bm_* ground truth).
% z = single-trial myChanZscore (per-freq, baseline -500..-100 ms, baseStd denom; n-invariant).
% Response per block x band x phase = median z over voxels pooled across events (freq x phase-time).
% Bands: high 40-58 Hz, low 25-40 Hz. Decomposition: superlet c1=3 ord[3 30] F=22:62.
if nargin<1, testN=Inf; end
proj='C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2';
addpath(fullfile(proj,'code')); addpath('C:/Users/Adam/Documents/GitHub/ZelanoLabScripts');
addpath('C:/Users/Adam/Documents/GitHub/Superlets/matlab-pure');
tdir=fullfile(proj,'out','tables');
F=22:1:62; C1=3; ORD=[3 30]; MULT=1; fs=500; EXCL={'CP'};
BANDS={'high',[40 58];'low',[25 40]};
PH={'inhaleRise','peakToExh','exhaleFall','troughToPause','pause'};
best=readtable(fullfile(tdir,'macbp_best.csv'),'TextType','string'); best=best(best.cohort=="OBE",:);
base='R:/Neurology/Zelano_Lab/Lab_Common/OBEControl/';
isAudio=@(s) strcmpi(s,'audio')|strcmpi(s,'audiobook');
isFocus=@(s) strcmpi(s,'focus')|strcmpi(s,'naturalFocus')|strcmpi(s,'slowFocus')|strcmpi(s,'focusedBreathing');
glob=containers.Map({'cueTask','threshTask','O15','breathingTask'},{'*cueTask*.mat','*threshold*.mat','*O15*.mat','*breathing*.mat'});

ROI=table(); DUR=table(); VAL=[];   % ROI medians, phase durations, validation errors
nrun=0;
for r=1:height(best)
  id=char(best.sessID(r)); part=char(best.participant(r)); bl=char(best.bestLabel(r)); tk=char(best.task(r));
  if any(strcmpi(part,EXCL)) || ~isKey(glob,tk), continue; end
  if nrun>=testN, break; end
  dd=dir(fullfile(base,id,'preProc',glob(tk))); if isempty(dd), continue; end
  fp=pick(dd); od=[]; for a=1:6, try, S=load(fp); f=fieldnames(S); od=S.(f{1}); clear S; break; catch, pause(6*a); end, end
  if isempty(od), fprintf('%s %s LOAD FAIL\n',id,tk); continue; end
  labs=cellfun(@(x)char(string(x)),od.labels,'uni',0); ci=find(strcmpi(labs,bl),1); if isempty(ci), continue; end
  x=double(od.data(ci,:)); x=fillmissing(x,'linear'); x=fillmissing(x,'nearest'); N=numel(x);
  P=slt_power_cont(x, fs, F, C1, ORD, MULT);
  % respiration (for flow-based landmarks)
  isR=cellfun(@(s)~isempty(regexpi(s,'rsp','once')),labs); rA=od.data(isR,:);
  ri=1; if isfield(od,'rspIDX'),ri=double(od.rspIDX); end; rf=1; if isfield(od,'rspFlip'),rf=double(od.rspFlip); end
  if ri>size(rA,1),ri=1; end
  rsp=double(rA(ri,:))*rf; rsp=fillmissing(rsp,'linear'); rsp=fillmissing(rsp,'nearest');
  sm=movmean(rsp,round(0.10*fs)); flow=gradient(sm);
  bd=od.behDat; fo=round(coerce(bd.finalOnset));

  % define blocks + per-event landmark matrix L (rows: [onset peak exhOn trough pauseOn end])
  blocks={};  % {name, category, L}
  if strcmp(tk,'breathingTask')
    tc=strtrim(string(getcol(bd,'task')));
    for cc={'audiobook','focusedBreathing'}
      sel = strcmp(cc{1},'audiobook') & isAudio(tc) | strcmp(cc{1},'focusedBreathing') & isFocus(tc);
      idx=find(sel(:).'); if isempty(idx), continue; end
      Lb=lm_bm(bd,idx); Lf=arrayfun(@(k) lm_flow(sm,flow,fo(k),nextOn(fo,k,N,fs),fs), idx,'uni',0);
      % validation: derived vs bm_* peak/exhOn/trough (ms)
      for j=1:numel(idx)
        if all(isfinite(Lb(j,[2 3 4]))) && ~isempty(Lf{j}) && all(isfinite(Lf{j}([2 3 4])))
          VAL=[VAL; (Lf{j}([2 3 4])-Lb(j,[2 3 4]))/fs*1000]; %#ok<AGROW>
        end
      end
      blocks(end+1,:)={cc{1},'nonolf',Lb}; %#ok<AGROW>
    end
  else
    io=fo(~isnan(fo)&fo>0); L=nan(numel(io),6);
    for k=1:numel(io), lk=lm_flow(sm,flow,io(k),nextOn(io,k,N,fs),fs); if ~isempty(lk), L(k,:)=lk; end, end
    blocks(end+1,:)={tk,'olf',L}; %#ok<AGROW>
  end

  for bb=1:size(blocks,1)
    bname=blocks{bb,1}; catg=blocks{bb,2}; L=blocks{bb,3};
    ok=all(isfinite(L(:,1:4)),2) & (L(:,2)>L(:,1)) & (L(:,3)>=L(:,2)) & (L(:,4)>L(:,3)) & (L(:,6)>L(:,4));
    L=L(ok,:); nB=size(L,1); if nB<5, fprintf('%s %s: %d valid events (skip)\n',id,bname,nB); continue; end
    % pooled baseline stats (-500..-100 ms before onset), per freq
    bs=[]; for k=1:nB, a=L(k,1)-round(0.5*fs); b=L(k,1)-round(0.1*fs); if a>=1&&b<=N, bs=[bs P(:,a:b)]; end, end %#ok<AGROW>
    baseMu=mean(bs,2); baseSd=std(bs,0,2); baseSd(baseSd<=0)=eps;
    for bi=1:size(BANDS,1)
      fm=F>=BANDS{bi,2}(1)&F<=BANDS{bi,2}(2);
      for p=1:5
        switch p
          case 1, s=L(:,1); e=L(:,2); case 2, s=L(:,2); e=L(:,3);
          case 3, s=L(:,3); e=L(:,4); case 4, s=L(:,4); e=L(:,5);
          case 5, s=L(:,5); e=L(:,6);
        end
        vox=[]; nev=0;
        for k=1:nB
          a=s(k); b=e(k); if ~isfinite(a)||~isfinite(b)||b<=a||a<1||b>N, continue; end
          z=(P(fm,a:b)-baseMu(fm))./baseSd(fm); vox=[vox; z(:)]; nev=nev+1; %#ok<AGROW>
        end
        if nev<3, med=NaN; else, med=median(vox,'omitnan'); end
        ROI=[ROI; table(string(id),string(part),string(bname),string(catg),string(BANDS{bi,1}), ...
             string(PH{p}),p,med,nev,nB,'VariableNames', ...
             {'sessID','participant','block','category','band','phase','phaseIdx','medZ','nEvents','nBreaths'})]; %#ok<AGROW>
      end
    end
    % phase durations (ms), medians per phase
    dur=[L(:,2)-L(:,1), L(:,3)-L(:,2), L(:,4)-L(:,3), L(:,5)-L(:,4), L(:,6)-L(:,5)]/fs*1000;
    DUR=[DUR; table(string(id),string(bname),string(catg), median(dur,1,'omitnan'), ...
         'VariableNames',{'sessID','block','category','medDurMs'})]; %#ok<AGROW>
    fprintf('%s | %-16s | %-6s | nB=%3d | durMs=[%s]\n', id,bname,catg,nB,num2str(round(median(dur,1,'omitnan'))));
  end
  nrun=nrun+1;
end
writetable(ROI, fullfile(tdir,'taskcmp_olf_phase_roi.csv'));
writetable(DUR, fullfile(tdir,'taskcmp_olf_phase_durations.csv'));
if ~isempty(VAL)
  fprintf('\nVALIDATION (flow detector - bm_* on breathing, ms): peak %+.0f/%.0f  exhOn %+.0f/%.0f  trough %+.0f/%.0f  (median/IQR)\n', ...
    median(VAL(:,1)), iqr(VAL(:,1)), median(VAL(:,2)), iqr(VAL(:,2)), median(VAL(:,3)), iqr(VAL(:,3)));
  fprintf('  median |err| ms: peak=%.0f exhOn=%.0f trough=%.0f\n', median(abs(VAL(:,1))),median(abs(VAL(:,2))),median(abs(VAL(:,3))));
end
fprintf('wrote taskcmp_olf_phase_roi.csv (%d rows), taskcmp_olf_phase_durations.csv\n', height(ROI));
end

% ---- landmark helpers ----
function L=lm_bm(bd,idx)   % [onset peak exhOn trough pauseOn end] in samples, from bm_*
g=@(nm) colv(bd,nm);
onset=g('finalOnset'); peak=g('bm_inhalePeaks'); exhOn=g('bm_exhaleOnsets');
trough=g('bm_exhaleTroughs'); pauseOn=g('bm_exhalePauseOnsets'); endt=g('endTim');
if all(isnan(peak)), peak=g('inMaxTim'); end
if all(isnan(trough)), trough=g('exMinTim'); end
L=[onset(idx) peak(idx) exhOn(idx) trough(idx) pauseOn(idx) endt(idx)];
% if no exhale pause, pauseOn := trough-recovery end approximated by end (phase5 empty)
np=~isfinite(L(:,5)); L(np,5)=L(np,6);
end
function L=lm_flow(sm,flow,t0,tCap,fs)  % derive [onset peak exhOn trough pauseOn end]
% breathmetrics-matched: exhale onset = downward crossing of the onset (baseline) level;
% end = next spontaneous inhale onset (so a cued sniff is segmented as ONE breath, not the
% whole inter-trial interval). tCap bounds the search (<=8 s past onset).
L=[]; if ~isfinite(t0)||t0<1, return; end
e=min(tCap, numel(sm)); if e<=t0+round(0.2*fs), return; end
base=sm(t0);
w1=t0:min(t0+round(2*fs),e); [~,rp]=max(sm(w1)); tIP=w1(rp);          % inhale peak
amp=sm(tIP)-base; if amp<=0, return; end
seg=tIP:e; dn=find(sm(seg)<=base,1,'first'); if isempty(dn), return; end; tEO=seg(dn);  % exhale onset = downward baseline crossing
w3=tEO:min(tEO+round(3*fs),e); [~,rt]=min(sm(w3)); tET=w3(rt);        % trough
w4=tET:e; up=find(sm(w4)>=base,1,'first'); if isempty(up), tPO=e; else, tPO=w4(up); end % pause onset = exhale recovery to baseline
w5=tPO:e; nx=find(sm(w5)>=base+0.10*amp,1,'first'); if isempty(nx), tEnd=e; else, tEnd=w5(nx); end % next inhale onset
L=[t0 tIP tEO tET tPO tEnd];
end
function t=nextOn(onsets,k,N,fs)   % search cap: next event onset, but never more than 8 s out
t=min(onsets(k)+round(8*fs), N);
if k<numel(onsets) && isfinite(onsets(k+1)), t=min(t, onsets(k+1)); end
end
function v=colv(bd,nm), if ~any(strcmp(bd.Properties.VariableNames,nm)), v=nan(height(bd),1); return; end
  c=bd.(nm); if iscell(c), v=cellfun(@(x) dbl(x), c); else, v=round(double(c)); end, end
function y=dbl(x), if isnumeric(x)&&~isempty(x), y=round(double(x(1))); else, y=round(double(string(x))); end, if isempty(y),y=NaN; end, end
function fp=pick(dd), n=lower({dd.name}); k=find(contains(n,'preproc'),1); if isempty(k),k=1; end, fp=fullfile(dd(k).folder,dd(k).name); end
function v=coerce(c)
if isnumeric(c), v=double(c(:)); return; end
if iscell(c), v=nan(numel(c),1); for i=1:numel(c), try, e=c{i}; if isnumeric(e)&&~isempty(e), v(i)=double(e(1)); else, v(i)=double(string(e)); end, catch, v(i)=NaN; end, end, return; end
try, v=double(string(c)); v=v(:); catch, v=nan(numel(c),1); end
end
function s=getcol(bd,nm), n=height(bd); s=repmat({''},n,1); if ~any(strcmp(bd.Properties.VariableNames,nm)), return; end
  v=bd.(nm); for i=1:n, try, if iscell(v), e=v{i}; else, e=v(i); end, t=string(e); if ~isscalar(t), if isempty(t), t=""; else, t=t(1); end, end, s{i}=char(t); catch, end, end, end
