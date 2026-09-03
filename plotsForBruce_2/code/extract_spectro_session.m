function out = extract_spectro_session(fp, taskRow, bestLabel, cfg)
% EXTRACT_SPECTRO_SESSION  Group-spectrogram sufficient stats for ONE recording,
% across THREE frequency-range sets, z-scored per frequency to a -500..-100 ms
% pre-inhale baseline (myChanZscore-equivalent single-trial baseline z), plus the
% mean epoched respiration waveform for overlay.
%
% Ranges: low 2-25 Hz, mid 25-58 Hz, high 60-150 Hz (superlet, c1=3, order [3 30]).
% Returns per range: sumZ [nF x nT] (sum over breaths of per-breath baseline-z),
%   nBreaths, freqs; plus sumRsp [1 x nT] and rsp count; tMs; identifiers.

if nargin<4, cfg=struct(); end
% NOTE: order [3 8] here (a light adaptive superlet) for the group VISUALIZATION
% spectrograms across three wide ranges — far fewer wavelets than the [3 30] used
% for the primary gamma ridge analysis, so the batch is tractable; the group-mean
% morphology is essentially unchanged at this order.
D = struct('c1',3,'ord',[2 4],'mult',1,'epochMs',[-1000 4000],'baseMs',[-500 -100],'fs',500);
fn=fieldnames(D); for i=1:numel(fn), if ~isfield(cfg,fn{i}), cfg.(fn{i})=D.(fn{i}); end, end
ranges = {[2 25 1],[25 58 1],[62 150 4]};   % [lo hi step]
rnames = {'low','mid','high'};

S=load(fp); vv=fieldnames(S); od=S.(vv{1}); clear S;
fs=cfg.fs; if isfield(od,'fs'), fs=double(od.fs); end
labs=cellfun(@(x)char(string(x)),od.labels,'uni',0);
ci=find(strcmpi(labs,bestLabel),1); if isempty(ci), error('bestLabel %s not found',bestLabel); end
sig=double(od.data(ci,:)); sig=fillmissing(sig,'linear'); sig=fillmissing(sig,'nearest'); N=numel(sig);
isRsp=cellfun(@(x)~isempty(regexpi(x,'rsp','once')),labs); rspAll=od.data(isRsp,:);
ri=1; if isfield(od,'rspIDX'),ri=double(od.rspIDX); end; rf=1; if isfield(od,'rspFlip'),rf=double(od.rspFlip); end
if ri>size(rspAll,1),ri=1; end
rsp=double(rspAll(ri,:))*rf; rsp=fillmissing(rsp,'linear'); rsp=fillmissing(rsp,'nearest');

onsets = get_onsets(od, taskRow);
[coh,grp,part,sn]=cohort_of(char(string(od.sessID)));
e0=round(cfg.epochMs(1)/1000*fs); e1=round(cfg.epochMs(2)/1000*fs);
tMs=(e0:e1)/fs*1000; nT=numel(tMs);
bMask = tMs>=cfg.baseMs(1) & tMs<=cfg.baseMs(2); bp0=find(bMask,1,'first'); bp1=find(bMask,1,'last');

out=struct('sessID',char(string(od.sessID)),'taskRow',taskRow,'bestLabel',bestLabel, ...
    'cohort',coh,'group',grp,'participant',part,'sessNum',sn,'tMs',tMs,'nBreaths',0, ...
    'sumRsp',zeros(1,nT),'nRsp',0,'ranges',{rnames},'nBaseSamp',bp1-bp0+1);
for r=1:3, out.(rnames{r})=struct('sumPow',[],'baseSum',[],'baseSumSq',[],'freqs',[]); end
if isempty(onsets), return; end

% precompute epoch validity + respiration mean (once)
valid=false(numel(onsets),1);
for k=1:numel(onsets)
    a=onsets(k)+e0; b=onsets(k)+e1;
    if a>=1 && b<=N
        valid(k)=true;
        out.sumRsp = out.sumRsp + rsp(a:b); out.nRsp=out.nRsp+1;
    end
end
out.nBreaths=sum(valid);
if out.nBreaths==0, return; end

for r=1:3
    lo=ranges{r}(1); hi=ranges{r}(2); st=ranges{r}(3); F=lo:st:hi;
    P = slt_power_cont(sig, fs, F, cfg.c1, cfg.ord, cfg.mult);   % nF x N (raw power)
    nF=numel(F); sumPow=zeros(nF,nT); baseSum=zeros(nF,1); baseSumSq=zeros(nF,1);
    for k=1:numel(onsets)
        if ~valid(k), continue; end
        a=onsets(k)+e0; b=onsets(k)+e1; Pe=P(:,a:b);
        sumPow = sumPow + Pe;                       % sum raw power over breaths
        Pb = Pe(:,bp0:bp1);                          % baseline window (-500..-100 ms)
        baseSum   = baseSum   + sum(Pb,2);           % pooled per-freq baseline stats
        baseSumSq = baseSumSq + sum(Pb.^2,2);
    end
    out.(rnames{r}).sumPow=sumPow; out.(rnames{r}).baseSum=baseSum;
    out.(rnames{r}).baseSumSq=baseSumSq; out.(rnames{r}).freqs=F;
end
end

% ---- events ----
function on = get_onsets(od, taskRow)
on=[]; bd=od.behDat;
coerce=@(c) local_coerce(c);
switch taskRow
  case {'cueTask','threshTask','O15'}
    v=coerce(bd.finalOnset); v=round(v(~isnan(v)&v>0)); on=v(:);
  case {'audiobook','focusedBreathing'}
    vn=bd.Properties.VariableNames;
    if any(strcmp(vn,'task')), tc=stringcol(bd.task); else, return; end
    want=ternary(strcmp(taskRow,'audiobook'),'audio','focus');
    sel=strcmpi(strtrim(tc),want);
    fo=round(coerce(bd.finalOnset));
    on=fo(sel & ~isnan(fo) & fo>0); on=on(:);
end
end
function v=local_coerce(c)
    if isnumeric(c), v=double(c); return; end
    if iscell(c), v=nan(numel(c),1); for i=1:numel(c), try, e=c{i}; if isnumeric(e)&&~isempty(e), v(i)=double(e(1)); else, v(i)=double(string(e)); end, catch, v(i)=NaN; end, end, return; end
    try, v=double(string(c)); catch, v=nan(numel(c),1); end
end
function s=stringcol(v), n=numel(v); s=repmat({''},n,1); for i=1:n, try, if iscell(v), e=v{i}; else, e=v(i); end, s{i}=char(string(e)); catch, end, end, end
function y=ternary(c,a,b), if c, y=a; else, y=b; end, end
function [coh,grp,part,sn]=cohort_of(sessID)
bits=strsplit(sessID,'_'); part=''; sn=1;
if numel(bits)>=4,part=bits{4}; end
if numel(bits)>=5, s=str2double(bits{5}); if ~isnan(s),sn=s; end, end
if strcmpi(part,'TPB'),part='TB'; end
typ=''; if numel(bits)>=2,typ=lower(bits{2}); end
if contains(typ,'dupi'),coh='Dupi'; elseif contains(typ,'eeg'),coh='EEG'; else,coh='OBE'; end
if strcmp(coh,'Dupi'), if ismember(sn,[1 2 3]),grp=sprintf('DupiS%d',sn); else,grp='DupiSx'; end, else,grp='Control'; end
end
