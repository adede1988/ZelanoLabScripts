function scratch_ridge_bandmax()
% SCRATCH: for each control breath, first 1000 ms after inhale onset, compute
% (1) the ridge frequency progression on the 25-58 Hz superlet ridge, MEAN-CENTERED
%     per breath (so we see the sweep shape, not the level), and
% (2) the whole-band max power = max over 25-58 Hz of the per-breath baseline-z at
%     each time point (peak gamma z anywhere in band).
% Aggregate to participant (mean over breaths), then grand mean over participants.
% Plot 2x2 (rows: ridge-freq-centered, band-max; cols: audiobook, focus), each
% participant a faint line, grand mean bold. Baseline for z = -500..-100 ms
% (matches the group spectrograms). Scratch only.
proj='C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2'; code=fullfile(proj,'code');
addpath(code); addpath('C:/Users/Adam/Documents/GitHub/ZelanoLabScripts');
addpath('C:/Users/Adam/Documents/GitHub/Superlets/matlab-pure');
assert(exist('slt_power_cont','file')==2 && exist('ridge_track','file')==2,'need superlet+ridge on path');

idx=readtable(fullfile(proj,'out','tables','session_index.csv'),'TextType','string');
best=readtable(fullfile(proj,'out','tables','macbp_best.csv'),'TextType','string');
idx=idx(idx.cohort=="OBE" & idx.task=="breathingTask" & idx.onDisk==1, :);

F=22:1:62; ridgeBand=[25 58]; c1=3; ord=[3 30]; mult=1; penalty=1.0; ridgeBW=2; fs=500;
bIdx=find(F>=ridgeBand(1)&F<=ridgeBand(2)); Fb=F(bIdx)';
preMs=[-500 -100]; winMs=[0 1000];
e0=round(-600/1000*fs); e1=round(1000/1000*fs); tMs=(e0:e1)/fs*1000; nT=numel(tMs);
preM=tMs>=preMs(1)&tMs<=preMs(2); wM=tMs>=winMs(1)&tMs<=winMs(2); tW=tMs(wM);
isFocus=@(s) strcmpi(s,'focus')|strcmpi(s,'naturalFocus')|strcmpi(s,'slowFocus');

parts={}; A_rf={}; F_rf={}; A_bm={}; F_bm={};
for r=1:height(idx)
    id=char(idx.sessID(r)); fp=char(idx.finalPath(r));
    m=strcmp(best.sessID,id)&strcmp(best.task,"breathingTask"); if ~any(m), continue; end
    bl=char(string(best.bestLabel(find(m,1))));
    p=strsplit(id,'_'); subj=p{4};
    fprintf('%s (%s) bl=%s ... ', id, subj, bl); t0=tic;
    od=[]; for a=1:3, try, S=load(fp); fn=fieldnames(S); od=S.(fn{1}); clear S; break; catch e, fprintf('[retry %d] ',a); pause(3); end, end
    if isempty(od), fprintf('LOAD FAIL\n'); continue; end
    labs=cellfun(@(x)char(string(x)),od.labels,'uni',0); ci=find(strcmpi(labs,bl),1);
    if isempty(ci), fprintf('no bl\n'); continue; end
    sig=double(od.data(ci,:)); sig=fillmissing(sig,'linear'); sig=fillmissing(sig,'nearest'); N=numel(sig);
    bd=od.behDat; tc=strtrim(string(getcol(bd,'task'))); on=round(coerce(bd.finalOnset));
    selA=strcmpi(tc,'audio'); selF=isFocus(tc); keep=~isnan(on)&on>0;
    P=slt_power_cont(sig, fs, F, c1, ord, mult);                 % nF x N
    [arf,abm]=deal([]); [frf,fbm]=deal([]);
    for kk=find((selA|selF)' & keep')
        a=on(kk)+e0; b=on(kk)+e1; if a<1||b>N, continue; end
        Pb=P(bIdx, a:b);                                        % band raw power (nFb x nT)
        muB=mean(Pb(:,preM),2); sdB=std(Pb(:,preM),0,2); sdB(sdB<=0)=eps;
        Zb=(Pb-muB)./sdB;                                       % baseline z
        [frW,~,~]=ridge_track(max(Zb,0), Fb, penalty, 1, ridgeBW); frW=frW(:).';
        bmax=max(Zb,[],1);                                      % peak z across band per time
        rfc=frW(wM)-mean(frW(wM),'omitnan');                   % mean-centered ridge freq over the window
        bmw=bmax(wM);
        if selA(kk), arf=[arf; rfc]; abm=[abm; bmw]; else, frf=[frf; rfc]; fbm=[fbm; bmw]; end %#ok<AGROW>
    end
    if isempty(arf)||isempty(frf), fprintf('no breaths one cond\n'); continue; end
    parts{end+1}=subj; %#ok<AGROW>
    A_rf{end+1}=mean(arf,1,'omitnan'); F_rf{end+1}=mean(frf,1,'omitnan');
    A_bm{end+1}=mean(abm,1,'omitnan'); F_bm{end+1}=mean(fbm,1,'omitnan');
    fprintf('nA=%d nF=%d (%.0fs)\n', size(arf,1), size(frf,1), toc(t0));
end
save(fullfile(proj,'out','tables','scratch_ridge_bandmax.mat'),'parts','A_rf','F_rf','A_bm','F_bm','tW');
plot_2x2(tW, parts, A_rf, F_rf, A_bm, F_bm, fullfile(proj,'out','figs','taskcmp','scratch_ridge_bandmax.png'));
fprintf('DONE scratch_ridge_bandmax (%d participants)\n', numel(parts));
end

function plot_2x2(tW, parts, A_rf, F_rf, A_bm, F_bm, outpng)
nS=numel(parts); cmap=lines(max(nS,3));
fig=figure('Position',[30 30 1300 900],'Color','w','Visible','off');
pan={A_rf,'audiobook',1,'ridge freq (mean-centered, Hz)'; F_rf,'focus',2,'ridge freq (mean-centered, Hz)';
     A_bm,'audiobook',3,'whole-band max z'; F_bm,'focus',4,'whole-band max z'};
% shared y per row
rf_all=cell2mat([A_rf(:);F_rf(:)]); bm_all=cell2mat([A_bm(:);F_bm(:)]);
yl_rf=[min(rf_all(:)) max(rf_all(:))]; yl_bm=[min(bm_all(:)) max(bm_all(:))];
for q=1:4
    subplot(2,2,pan{q,3}); hold on; C=pan{q,1};
    for i=1:nS, plot(tW, C{i}, '-','Color',[cmap(i,:) 0.35],'LineWidth',1); end
    M=mean(cell2mat(C(:)),1,'omitnan'); plot(tW, M, '-','Color',[0 0 0],'LineWidth',2.5);
    yline(0,'k:'); xlim([0 1000]);
    if pan{q,3}<=2, ylim(yl_rf); else, ylim(yl_bm); end
    title(sprintf('%s — %s', regexprep(pan{q,4},' \(.*',''), pan{q,2}),'FontSize',12,'FontWeight','bold');
    xlabel('ms after inhale onset'); ylabel(pan{q,4}); set(gca,'FontSize',10);
    if pan{q,3}==2, legend([repmat({''},1,nS), {'grand mean'}],'Location','best'); end
end
sgtitle('Control: first 1000 ms — ridge-frequency sweep & whole-band max z (faint = participant)','FontSize',13,'FontWeight','bold');
exportgraphics(fig,outpng,'Resolution',140); close(fig);
end

function v=coerce(c)
    if isnumeric(c), v=double(c); return; end
    if iscell(c), v=nan(numel(c),1); for i=1:numel(c), try, e=c{i}; if isnumeric(e)&&~isempty(e), v(i)=double(e(1)); else, v(i)=double(string(e)); end, catch, v(i)=NaN; end, end, return; end
    try, v=double(string(c)); catch, v=nan(numel(c),1); end
end
function s=getcol(bd,nm)
    n=height(bd); s=repmat({''},n,1); if ~any(strcmp(bd.Properties.VariableNames,nm)), return; end
    v=bd.(nm); for i=1:n, try, if iscell(v), e=v{i}; else, e=v(i); end, t=string(e); if ~isscalar(t), if isempty(t), t=""; else, t=t(1); end, end, s{i}=char(t); catch, end, end
end
