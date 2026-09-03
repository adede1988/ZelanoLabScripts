function out = extract_gamma_session(fp, taskRow, bestLabel, cfg)
% EXTRACT_GAMMA_SESSION  Per-breath/sniff gamma characterization on one recording.
%   fp        : final .mat path
%   taskRow   : one of 'cueTask','threshTask','O15','audiobook','focusedBreathing'
%   bestLabel : macBP channel label to analyze (e.g. 'macBP2')
%   cfg       : struct of params (see defaults below)
%
% Returns out with: perBreath (table), meanZwhole/meanZbase (nF x nT sums),
%   nBreaths, tMs, freqs, sessID/cohort/group/participant/sessNum.
%
% Self-contained: uses slt_power_cont, ridge_track (this folder), myChanZscore
% & fooof_basic (ZelanoLabScripts repo root), all added to path by the runner.

if nargin<4, cfg = struct(); end
D = struct('F',22:1:62,'ridgeBand',[25 58],'c1',3,'ord',[3 30],'mult',1, ...
   'epochMs',[-1000 4000],'baseMs',[-1000 -250],'burstZ',3,'freqFloorZ',2, ...
   'peakWinMs',[0 2000],'chirpWinMs',[0 500],'chirpWin2Ms',[200 1000],'penalty',1.0,'ridgeBW',2, ...
   'peakSearchMs',2500,'fs',500);
fn=fieldnames(D); for i=1:numel(fn), if ~isfield(cfg,fn{i}), cfg.(fn{i})=D.(fn{i}); end, end

S=load(fp); vv=fieldnames(S); od=S.(vv{1}); clear S;
fs = cfg.fs; if isfield(od,'fs'), fs=double(od.fs); end
labs = cellfun(@(x)char(string(x)),od.labels,'uni',0);
ci = find(strcmpi(labs,bestLabel),1);
if isempty(ci), error('bestLabel %s not found', bestLabel); end
sig = double(od.data(ci,:)); sig = fillmissing(sig,'linear'); sig=fillmissing(sig,'nearest');
N = numel(sig);

% respiration trace (inhale up)
isRsp = cellfun(@(x) ~isempty(regexpi(x,'rsp','once')), labs);
rspAll = od.data(isRsp,:);
ri = 1; if isfield(od,'rspIDX'), ri=double(od.rspIDX); end
rf = 1; if isfield(od,'rspFlip'), rf=double(od.rspFlip); end
if ri>size(rspAll,1), ri=1; end
rsp = double(rspAll(ri,:))*rf; rsp=fillmissing(rsp,'linear'); rsp=fillmissing(rsp,'nearest');

% ---- events + landmarks ----
[ev, cov] = get_events(od, taskRow, rsp, fs, cfg);   % ev: struct array of landmark samples
nEv = numel(ev);
if nEv==0
    out = emptyOut(od, taskRow, bestLabel, cfg); return;
end

% ---- continuous superlet power TFR ----
P = slt_power_cont(sig, fs, cfg.F, cfg.c1, cfg.ord, cfg.mult);   % nF x N (power)
F = cfg.F(:); nF=numel(F);
bIdx = find(F>=cfg.ridgeBand(1) & F<=cfg.ridgeBand(2));
Fb = F(bIdx);

% ---- session-level respiration-gamma coupling (Tort MI, preferred phase, inhale/exhale) ----
sessCoup = session_coupling(P(bIdx,:), rsp, fs);

% epoch sample offsets
e0 = round(cfg.epochMs(1)/1000*fs); e1 = round(cfg.epochMs(2)/1000*fs);
tMs = (e0:e1)/fs*1000; nT=numel(tMs);
bMask = tMs>=cfg.baseMs(1) & tMs<=cfg.baseMs(2); bp0=find(bMask,1,'first'); bp1=find(bMask,1,'last');
dt_s = 1/fs;

% spectrogram accumulators over the FULL freq range:
%  sumZwhole  = sum over breaths of the per-breath single-trial within-frequency z
%               (== myChanZscore, no baseline, per breath) -> group map = sumZwhole/nBreaths
%  sumRawMap/poolSumsqRaw = raw-power sufficient stats (pooled myChanZscore variant / FOOOF-on-mean)
sumZwhole = zeros(nF, nT); sumRawMap = zeros(nF, nT); poolSumsqRaw = zeros(nF,1); nUsed=0;

rows = cell(nEv,1);
for k=1:nEv
    on = ev(k).onset;
    a = on+e0; b = on+e1;
    if a<1 || b>N, continue; end          % epoch off-recording
    sliceP = P(:, a:b);                     % nF x nT (raw power)
    Pb = sliceP(bIdx,:);                    % band raw power

    % --- per-frequency single-trial z (analytic = myChanZscore single-column limit) ---
    mu = mean(Pb,2); sdv = std(Pb,0,2); sdv(sdv<=0)=eps;
    Zw = (Pb - mu)./sdv;                                    % whole-window z
    muB = mean(Pb(:,bp0:bp1),2); sdB = std(Pb(:,bp0:bp1),0,2); sdB(sdB<=0)=eps;
    Zb = (Pb - muB)./sdB;                                   % pre-inhale baseline z

    % --- ridge on whole-window z (primary) --- (force row orientation to match tMs)
    [frW, idW, zpW] = ridge_track(max(Zw,0), Fb, cfg.penalty, 1, cfg.ridgeBW);
    frW=frW(:).'; idW=idW(:).'; zpW=zpW(:).';
    rpRaw = Pb(sub2ind(size(Pb), idW, 1:nT));   % 1 x nT raw power on ridge
    bandPow = mean(Pb,1);                         % 1 x nT mean 25-58 raw power
    % ridge on baseline z (robustness)
    [~, ~, zpB] = ridge_track(max(Zb,0), Fb, cfg.penalty, 1, cfg.ridgeBW);
    zpB=zpB(:).';

    % accumulate spectrogram stats (full freq range)
    muF=mean(sliceP,2); sdF=std(sliceP,0,2); sdF(sdF<=0)=eps;
    sumZwhole = sumZwhole + (sliceP - muF)./sdF;    % per-breath single-trial within-freq z
    sumRawMap = sumRawMap + sliceP;
    poolSumsqRaw = poolSumsqRaw + sum(sliceP.^2,2);
    nUsed=nUsed+1;

    % --- metrics ---
    M = struct();
    M.sessID=ev(k).sessID; M.task=taskRow; M.cohort=ev(k).cohort; M.group=ev(k).group;
    M.participant=ev(k).participant; M.sessNum=ev(k).sessNum;
    M.breathIdx=k; M.onsetSample=on; M.goodBreath=ev(k).goodBreath; M.nEv=nEv;
    % landmarks (ms rel onset)
    M.inhalePeakMs = smp2ms(ev(k).inhalePeak,on,fs);
    M.returnCrossMs= smp2ms(ev(k).returnCross,on,fs);
    M.exhaleStartMs= smp2ms(ev(k).exhaleStart,on,fs);
    M.exhaleTroughMs=smp2ms(ev(k).exhaleTrough,on,fs);
    M.endSymMs     = smp2ms(ev(k).endSym,on,fs);
    M.breathLenMs  = ev(k).breathLenMs;
    % covariates
    M.inhaleVolume=cov(k).inhaleVolume; M.inhaleDuration=cov(k).inhaleDuration;
    M.peakInspFlow=cov(k).peakInspFlow; M.breathLength=cov(k).breathLength;

    % fixed time windows
    W = {[-500 0],[0 500],[500 1000],[1000 1500],[1500 2000]};
    for wi=1:numel(W)
        idx = tMs>=W{wi}(1) & tMs<W{wi}(2);
        m = winmetrics(idx, frW, zpW, rpRaw, bandPow, dt_s, cfg.freqFloorZ);
        M = setwin(M, sprintf('w%d',wi), m);
    end
    % phase segments (samples rel onset -> ms windows)
    segs = phase_segs(ev(k), on, fs, cfg);
    for si=1:numel(segs)
        idx = tMs>=segs{si}(1) & tMs<segs{si}(2);
        m = winmetrics(idx, frW, zpW, rpRaw, bandPow, dt_s, cfg.freqFloorZ);
        M = setwin(M, sprintf('p%d',si), m);
    end

    % peak / bursts over response window [0 2000] using whole-window z
    rwin = tMs>=cfg.peakWinMs(1) & tMs<=cfg.peakWinMs(2);
    [M.peakZ, M.peakLatMs, M.peakFreq] = peakmetrics(rwin, tMs, zpW, frW);
    bm = burstmetrics(rwin, tMs, zpW, cfg.burstZ, dt_s);
    M.anyBurst=bm.any; M.burstLatMs=bm.lat; M.timeAboveMs=bm.above; M.nBursts=bm.n;
    M.dutyCycle=bm.duty; M.maxBurstMs=bm.maxdur;
    % baseline-z robustness
    [M.peakZ_base,~,~] = peakmetrics(rwin, tMs, zpB, frW);
    bmB = burstmetrics(rwin, tMs, zpB, cfg.burstZ, dt_s);
    M.anyBurst_base=bmB.any; M.timeAboveMs_base=bmB.above;

    % frequency dynamics
    cidx = tMs>=cfg.chirpWinMs(1) & tMs<=cfg.chirpWinMs(2);
    M.chirpSlope = linslope(tMs(cidx)/1000, frW(cidx));   % Hz/s over inhale window (0-500 ms)
    cidx2 = tMs>=cfg.chirpWin2Ms(1) & tMs<=cfg.chirpWin2Ms(2);
    M.chirpSlope2 = linslope(tMs(cidx2)/1000, frW(cidx2)); % Hz/s over 200-1000 ms
    gidx = rwin & (zpW>cfg.freqFloorZ);
    if nnz(gidx)>=3
        M.freqSpan = max(frW(gidx))-min(frW(gidx));
        M.freqJitter = std(frW(gidx));
    else, M.freqSpan=NaN; M.freqJitter=NaN; end

    % per-breath FOOOF-lite (broadband spectrum on the response window raw signal)
    rr = (on + round(cfg.peakWinMs(1)/1000*fs)) : (on + round(cfg.peakWinMs(2)/1000*fs));
    rr = rr(rr>=1 & rr<=N);
    fl = fooof_lite(sig(rr), fs, cfg.ridgeBand);
    M.apExp=fl.exp; M.apOffset=fl.offset; M.gammaBumpDb=fl.bumpDb; M.gammaPeakPresent=fl.peak;

    rows{k}=M;
end
rows = rows(~cellfun(@isempty,rows));
perBreath = struct2table([rows{:}]);

out = struct();
out.perBreath = perBreath;
out.sumZwhole = sumZwhole; out.sumRawMap = sumRawMap; out.poolSumsqRaw = poolSumsqRaw;
out.nBreaths = nUsed;
out.tMs = tMs; out.freqs = F; out.bandFreqs = Fb; out.nT = nT;
out.coupling = sessCoup;
if ~isempty(rows)
    out.sessID=rows{1}.sessID; out.cohort=rows{1}.cohort; out.group=rows{1}.group;
    out.participant=rows{1}.participant; out.sessNum=rows{1}.sessNum;
end
end

% ================= session respiration-gamma coupling =================
function s = session_coupling(Pband, rsp, fs)
s = struct('MI',NaN,'prefPhaseRad',NaN,'resultantLen',NaN,'rayleighP',NaN,'inhExhRatio',NaN);
try
    gammaAmp = mean(sqrt(max(Pband,0)),1);      % 1 x N gamma amplitude envelope
    rr = rsp - mean(rsp);
    [bb,aa] = butter(2,[0.05 1]/(fs/2),'bandpass');
    rf = filtfilt(bb,aa,double(rr));
    ph = angle(hilbert(rf));                      % respiration phase
    nb=18; edges=linspace(-pi,pi,nb+1); mAmp=zeros(1,nb);
    for b=1:nb, in=ph>=edges(b) & ph<edges(b+1); if any(in), mAmp(b)=mean(gammaAmp(in)); end, end
    p=mAmp/sum(mAmp); p(p<=0)=eps; H=-sum(p.*log(p));
    s.MI=(log(nb)-H)/log(nb);                     % Tort modulation index
    w=gammaAmp; R=sum(w.*exp(1i*ph))/sum(w);
    s.prefPhaseRad=angle(R); s.resultantLen=abs(R);
    neff=sum(w)^2/sum(w.^2); s.rayleighP=exp(-neff*s.resultantLen^2);
    dphi=[0 diff(rf)]; inh=dphi>0;
    s.inhExhRatio = mean(gammaAmp(inh))/(mean(gammaAmp(~inh))+eps);
catch
end
end

% ================= events + landmarks =================
function [ev, cov] = get_events(od, taskRow, rsp, fs, cfg)
sessID=char(string(od.sessID)); [coh,grp,part,sn]=cohort_of(sessID, od);
bd = od.behDat;
ev=struct([]); cov=struct([]);
switch taskRow
  case {'cueTask','threshTask','O15'}
    on = coerce(bd.finalOnset);
    keep = ~isnan(on) & on>0;
    on = round(on(keep));
    for k=1:numel(on)
        lm = derive_landmarks(rsp, on(k), fs, cfg);
        ev = append_ev(ev, sessID,coh,grp,part,sn, on(k), lm, 1);
        cov(end+1).inhaleVolume = lm.inhaleVolume; %#ok<AGROW>
        cov(end).inhaleDuration = lm.inhaleDuration;
        cov(end).peakInspFlow   = lm.peakInspFlow;
        cov(end).breathLength   = lm.breathLength;
    end
  case {'audiobook','focusedBreathing'}
    vn=bd.Properties.VariableNames;
    tcol = getstr(bd,'task');
    want = ternary(strcmp(taskRow,'audiobook'),'audio','focus');
    sel = strcmpi(strtrim(tcol), want);
    if ~any(sel), ev=struct([]); cov=struct([]); return; end
    on = round(coerce(bd.finalOnset));
    inMaxTim = round(coerce(getnum(bd,'inMaxTim')));
    endTim   = round(coerce(getnum(bd,'endTim')));
    exMinTim = round(coerce(getnum(bd,'exMinTim')));
    Yonset   = coerce(getnum(bd,'Yonset'));
    gb       = coerce(getnum(bd,'goodBreath'));
    len      = coerce(getnum(bd,'length'));
    iv=coerce(getnum(bd,'bm_inhaleVolumes')); idur=coerce(getnum(bd,'bm_inhaleDurations'));
    pif=coerce(getnum(bd,'bm_peakInspiratoryFlows'));
    idxs=find(sel(:).');
    for k=idxs
        lm=struct('inhalePeak',inMaxTim(k),'exhaleStart',endTim(k),'exhaleTrough',exMinTim(k), ...
            'onsetY',Yonset(k));
        lm.returnCross = return_cross(rsp, on(k), inMaxTim(k), Yonset(k));
        if isnan(lm.exhaleStart), lm.exhaleStart=lm.returnCross; end
        lm.endSym = exMinTim(k) + max(0,(exMinTim(k)-on(k)));
        lm.breathLenMs = len(k)*1000;
        gbk=gb(k); if isnan(gbk), gbk=1; end
        ev = append_ev(ev, sessID,coh,grp,part,sn, on(k), lm, gbk);
        cov(end+1).inhaleVolume=iv(k); %#ok<AGROW>
        cov(end).inhaleDuration=idur(k); cov(end).peakInspFlow=pif(k); cov(end).breathLength=len(k);
    end
  otherwise
    error('unknown taskRow %s', taskRow);
end
end

function ev = append_ev(ev, sessID,coh,grp,part,sn, on, lm, gb)
    i=numel(ev)+1;
    ev(i).sessID=sessID; ev(i).cohort=coh; ev(i).group=grp; ev(i).participant=part; ev(i).sessNum=sn;
    ev(i).onset=on; ev(i).goodBreath=gb;
    ev(i).inhalePeak=getf(lm,'inhalePeak'); ev(i).exhaleStart=getf(lm,'exhaleStart');
    ev(i).exhaleTrough=getf(lm,'exhaleTrough'); ev(i).returnCross=getf(lm,'returnCross');
    ev(i).endSym=getf(lm,'endSym');
    bl=getf(lm,'breathLenMs'); if isempty(bl)||isnan(bl), bl=NaN; end; ev(i).breathLenMs=bl;
end

function lm = derive_landmarks(rsp, on, fs, cfg)
% derive inhale peak / return crossing / exhale trough from rsp trace for a sniff
lm=struct('inhalePeak',NaN,'exhaleStart',NaN,'exhaleTrough',NaN,'returnCross',NaN,...
    'onsetY',NaN,'endSym',NaN,'breathLenMs',NaN,'inhaleVolume',NaN,'inhaleDuration',NaN,...
    'peakInspFlow',NaN,'breathLength',NaN);
N=numel(rsp); w=round(cfg.peakSearchMs/1000*fs);
if on<1||on>N, return; end
lm.onsetY = rsp(on);
seg_end = min(N, on+w);
seg = rsp(on:seg_end);
% inhale peak = max of trace in window (assumes inhale up)
[~,pk] = max(seg); lm.inhalePeak = on+pk-1;
% return crossing: first below onsetY after peak
lm.returnCross = return_cross(rsp, on, lm.inhalePeak, lm.onsetY);
lm.exhaleStart = lm.returnCross;
% exhale trough = min after exhaleStart within window
if ~isnan(lm.exhaleStart)
    s2 = lm.exhaleStart; e2 = min(N, s2 + w);
    if e2>s2, seg2=rsp(s2:e2); [~,tr]=min(seg2); lm.exhaleTrough=s2+tr-1; end
end
if ~isnan(lm.exhaleTrough), lm.endSym = lm.exhaleTrough + max(0,(lm.exhaleTrough-on)); end
% covariates (derived): inhaleDuration onset->peak; peakInspFlow ~ max diff; inhaleVolume ~ trapz
lm.inhaleDuration = (lm.inhalePeak-on)/fs;
if lm.inhalePeak>on
    d = diff(rsp(on:lm.inhalePeak)); lm.peakInspFlow = max(d)*fs;
    lm.inhaleVolume = trapz(rsp(on:lm.inhalePeak)-lm.onsetY)/fs;
end
if ~isnan(lm.exhaleTrough), lm.breathLength=(lm.exhaleTrough-on)/fs; lm.breathLenMs=lm.breathLength*1000; end
end

function rc = return_cross(rsp, on, pk, onsetY)
rc=NaN; N=numel(rsp);
if isnan(pk)||isnan(onsetY)||pk<1||pk>N, return; end
for t=pk+1:min(N, pk+round(3*500))    % search up to 3 s after peak
    if rsp(t) <= onsetY, rc=t; return; end
end
end

% ================= metric helpers =================
function m = winmetrics(idx, fr, zp, rp, bp, dt_s, floorZ)
m=struct('rfreqPW',NaN,'rfreqGated',NaN,'rfreqRaw',NaN,'rpowZ',NaN,'rpowDb',NaN,'bandDb',NaN,'aucZ',NaN);
if ~any(idx), return; end
f=fr(idx); z=zp(idx); r=rp(idx); b=bp(idx);
wpos=max(z,0);
if sum(wpos)>0, m.rfreqPW = sum(f.*wpos)/sum(wpos); end
g = z>floorZ; if any(g), m.rfreqGated=mean(f(g)); end
m.rfreqRaw=mean(f);
m.rpowZ=mean(z);
m.rpowDb=10*log10(max(mean(r),eps));
m.bandDb=10*log10(max(mean(b),eps));
m.aucZ=sum(wpos)*dt_s;
end

function M = setwin(M, pre, m)
fn=fieldnames(m);
for i=1:numel(fn), M.([pre '_' fn{i}]) = m.(fn{i}); end
end

function [pz,plat,pf] = peakmetrics(rwin, tMs, zp, fr)
pz=NaN;plat=NaN;pf=NaN;
if ~any(rwin), return; end
z=zp; z(~rwin)=-inf; [pz,pk]=max(z); if ~isfinite(pz), pz=NaN; return; end
plat=tMs(pk); pf=fr(pk);
end

function bm = burstmetrics(rwin, tMs, zp, thr, dt_s)
bm=struct('any',0,'lat',NaN,'above',0,'n',0,'duty',NaN,'maxdur',0);
z=zp(rwin); tt=tMs(rwin);
b = z>thr;
if ~any(b), bm.duty=0; return; end
bm.any=1; bm.above=sum(b)*dt_s*1000; bm.duty=mean(b);
fe=find(diff([0 b(:).'])==1); bm.n=numel(fe);
bm.lat = tt(find(b,1,'first')) - 0;   % ms rel onset (window starts at 0)
% max run
d=diff([0 b(:).' 0]); st=find(d==1); en=find(d==-1)-1; runs=en-st+1;
bm.maxdur = max(runs)*dt_s*1000;
end

function s = linslope(x,y)
x=x(:); y=y(:); ok=isfinite(x)&isfinite(y);
if nnz(ok)<3, s=NaN; return; end
p=polyfit(x(ok),y(ok),1); s=p(1);
end

function segs = phase_segs(e, on, fs, cfg)
% four phase segments as ms windows [start end] rel onset
g=@(smp) (smp-on)/fs*1000;
segs={};
p1s=0; p1e = valid(g(e.inhalePeak));
segs{1}=[p1s, def(p1e, 500)];
p2s=valid(g(e.inhalePeak)); p2e=valid(g(e.returnCross));
segs{2}=[def(p2s,segs{1}(2)), def(p2e, def(p2s,0)+500)];
p3s=valid(g(e.exhaleStart)); p3e=valid(g(e.exhaleTrough));
segs{3}=[def(p3s,segs{2}(2)), def(p3e, def(p3s,0)+500)];
p4s=valid(g(e.exhaleTrough)); p4e=valid(g(e.endSym));
segs{4}=[def(p4s,segs{3}(2)), def(p4e, def(p4s,0)+500)];
% clamp to epoch
for i=1:4
    segs{i}(1)=max(cfg.epochMs(1), segs{i}(1));
    segs{i}(2)=min(cfg.epochMs(2), segs{i}(2));
    if segs{i}(2)<=segs{i}(1), segs{i}=[NaN NaN]; end
end
end

% ================= FOOOF-lite (fast per-breath) =================
function fl = fooof_lite(x, fs, band)
fl=struct('exp',NaN,'offset',NaN,'bumpDb',NaN,'peak',0);
x=x(:); if numel(x)<fs*0.5, return; end
x=detrend(x,0);
w=min(numel(x), round(1.0*fs)); if mod(w,2)==1,w=w+1; end
[pxx,f]=pwelch(x, hann(w,'periodic'), round(w/2), max(1024,2^nextpow2(w)), fs);
in = f>=2 & f<=58;
lf=log10(f(in)); lp=log10(pxx(in)+eps);
% robust 1/f fit excluding the gamma band and the 55-62 line region
excl = (f(in)>=band(1) & f(in)<=band(2)) | (f(in)>=55);
mask=~excl;
if nnz(mask)<8, return; end
pcoef=polyfit(lf(mask), lp(mask),1);
fl.exp=-pcoef(1); fl.offset=pcoef(2);
flat = lp - (pcoef(1)*lf+pcoef(2));   % log10 over 1/f
gb = f(in)>=band(1) & f(in)<=band(2);
if any(gb)
    fl.bumpDb = 10*max(flat(gb));
    nz = 1.4826*mad(flat(mask),1);
    fl.peak = double(fl.bumpDb > 10*(2*nz));   % >2 SD bump
end
end

% ================= small utils =================
function v=coerce(c)
    if isnumeric(c), v=double(c); return; end
    if iscell(c)
        v=nan(numel(c),1);
        for i=1:numel(c)
            try, e=c{i}; if isnumeric(e)&&~isempty(e), v(i)=double(e(1)); else, v(i)=double(string(e)); end, catch, v(i)=NaN; end
        end
        return;
    end
    try, v=double(string(c)); catch, v=nan(numel(c),1); end
end
function s=getstr(bd,nm)
    n=height(bd); s=repmat({''},n,1);
    if ~any(strcmp(bd.Properties.VariableNames,nm)), return; end
    v=bd.(nm);
    for i=1:n
        try
            if iscell(v), e=v{i}; else, e=v(i); end
            if isempty(e), s{i}=''; else, s{i}=char(string(e)); end
        catch, s{i}=''; end
    end
end
function v=getnum(bd,nm), if any(strcmp(bd.Properties.VariableNames,nm)), v=bd.(nm); else, v=nan(height(bd),1); end, end
function y=ternary(c,a,b), if c, y=a; else, y=b; end, end
function v=getf(s,f), if isfield(s,f), v=s.(f); else, v=NaN; end, end
function ms=smp2ms(smp,on,fs), if isnan(smp), ms=NaN; else, ms=(smp-on)/fs*1000; end, end
function y=valid(x), if isempty(x)||~isfinite(x), y=NaN; else, y=x; end, end
function y=def(x,d), if isnan(x), y=d; else, y=x; end, end

function [coh,grp,part,sn]=cohort_of(sessID, od)
bits=strsplit(sessID,'_'); part=''; sn=1;
if numel(bits)>=4, part=bits{4}; end
if numel(bits)>=5, s=str2double(bits{5}); if ~isnan(s), sn=s; end, end
if strcmpi(part,'TPB'), part='TB'; end
typ=''; if numel(bits)>=2, typ=lower(bits{2}); end
if contains(typ,'dupi'), coh='Dupi'; elseif contains(typ,'eeg'), coh='EEG'; else, coh='OBE'; end
if strcmp(coh,'Dupi')
    if ismember(sn,[1 2 3]), grp=sprintf('DupiS%d',sn); else, grp='DupiSx'; end
else, grp='Control'; end
end

function out=emptyOut(od, taskRow, bestLabel, cfg)
out=struct('perBreath',table(),'sumZwhole',[],'sumRawMap',[],'poolSumsqRaw',[],'nBreaths',0,...
    'tMs',[],'freqs',[],'bandFreqs',[],'nT',0);
end
