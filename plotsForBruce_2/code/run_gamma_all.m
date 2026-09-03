function run_gamma_all(useParfor, rootOverride)
% RUN_GAMMA_ALL  Batch per-breath/sniff gamma extraction over all Dupi+OBE
% recordings x 5 task-rows (cueTask, threshTask, O15, audiobook, focusedBreathing).
% Parfor-safe: each work unit writes its own per-unit CSV + spectrogram .mat.
%
%   run_gamma_all()                 % serial
%   run_gamma_all(true)             % parfor
%   run_gamma_all(true,'E:\Lab_Common')  % rebase R:\...\Lab_Common -> E:\Lab_Common
%
% Self-contained; adds repo root + Superlets to path. Reads:
%   out/tables/session_index.csv, out/tables/macbp_best.csv
% Writes:
%   out/gamma/perbreath/<sessID>__<taskRow>.csv
%   out/gamma/spectro/<sessID>__<taskRow>.mat  (sumRawMap,poolSumsqRaw,nBreaths,tMs,freqs,meta)

if nargin<1||isempty(useParfor), useParfor=false; end
if nargin<2, rootOverride=''; end

here = fileparts(mfilename('fullpath'));
proj = fileparts(here); repo = fileparts(proj);
addpath(repo); addpath(here);
% ensure external tools on path (ZelanoLabScripts repo for fooof_basic/myChanZscore; Superlets)
zls = {repo,'C:\Users\Adam\Documents\GitHub\ZelanoLabScripts','E:\GitHub\ZelanoLabScripts','G:\My Drive\GitHub\ZelanoLabScripts'};
for c=zls, if exist(fullfile(c{1},'fooof_basic.m'),'file'), addpath(c{1}); end, end
for c = {'C:\Users\Adam\Documents\GitHub\Superlets\matlab-pure','E:\GitHub\Superlets\matlab-pure','G:\My Drive\GitHub\Superlets\matlab-pure'}
    if exist(c{1},'dir'), addpath(c{1}); end
end
assert(exist('fooof_basic','file')==2, 'fooof_basic not on path');
assert(exist('myChanZscore','file')==2, 'myChanZscore not on path');
assert(exist('faslt','file')==2 || exist('slt_power_cont','file')==2, 'superlet not on path');
tdir = fullfile(proj,'out','tables');
gdir = fullfile(proj,'out','gamma'); pbdir=fullfile(gdir,'perbreath'); spdir=fullfile(gdir,'spectro'); cpdir=fullfile(gdir,'coupling');
if ~exist(pbdir,'dir'), mkdir(pbdir); end
if ~exist(spdir,'dir'), mkdir(spdir); end
if ~exist(cpdir,'dir'), mkdir(cpdir); end

idx = readtable(fullfile(tdir,'session_index.csv'),'TextType','string');
best= readtable(fullfile(tdir,'macbp_best.csv'),'TextType','string');
idx = idx( (idx.cohort=="Dupi"|idx.cohort=="OBE") & idx.onDisk==1, :);

% build work units
U = {};  % {fp, taskRow, bestLabel, sessID}
for r=1:height(idx)
    task=char(idx.task(r)); id=char(idx.sessID(r)); fp=char(idx.finalPath(r));
    fp = resolve_path(fp, rootOverride);
    bl = lookup_best(best, id, task);
    if isempty(bl), fprintf('no bestMac for %s %s -- skip\n', id, task); continue; end
    switch task
        case 'cueTask',    U(end+1,:)={fp,'cueTask',bl,id}; %#ok<AGROW>
        case 'threshTask', U(end+1,:)={fp,'threshTask',bl,id}; %#ok<AGROW>
        case 'O15',        U(end+1,:)={fp,'O15',bl,id}; %#ok<AGROW>
        case 'breathingTask'
            U(end+1,:)={fp,'audiobook',bl,id}; %#ok<AGROW>
            U(end+1,:)={fp,'focusedBreathing',bl,id}; %#ok<AGROW>
    end
end
nU = size(U,1);
fprintf('gamma batch: %d work units, parfor=%d\n', nU, useParfor);

if useParfor
    nw = max(1, min(6, feature('numcores')-2));
    delete(gcp('nocreate'));
    parpool('Processes', nw);
    fprintf('parpool up: %d workers\n', nw);
    parfor i=1:nU
        process_unit(U(i,:), pbdir, spdir, cpdir);
    end
    delete(gcp('nocreate'));
else
    for i=1:nU, fprintf('[%3d/%3d] %s %s\n', i, nU, U{i,4}, U{i,2}); process_unit(U(i,:), pbdir, spdir, cpdir); end
end
fprintf('DONE run_gamma_all\n');
end

function process_unit(u, pbdir, spdir, cpdir)
fp=u{1}; taskRow=u{2}; bl=u{3}; id=u{4};
outcsv = fullfile(pbdir, sprintf('%s__%s.csv', id, taskRow));
outmat = fullfile(spdir, sprintf('%s__%s.mat', id, taskRow));
outcp  = fullfile(cpdir, sprintf('%s__%s.csv', id, taskRow));
if exist(outcsv,'file') && exist(outmat,'file') && exist(outcp,'file'), return; end   % resume
try
    out = extract_gamma_session(fp, taskRow, bl, struct());
    if height(out.perBreath)>0
        writetable(out.perBreath, outcsv);
    end
    if isfield(out,'coupling')
        cp=out.coupling; [coh,grp,part,sn]=deal("","","",NaN);
        if isfield(out,'cohort'), coh=string(out.cohort); grp=string(out.group); part=string(out.participant); sn=out.sessNum; end
        crow=table(string(id),string(taskRow),string(bl),coh,grp,part,sn, ...
            cp.MI,cp.prefPhaseRad,cp.resultantLen,cp.rayleighP,cp.inhExhRatio, ...
            'VariableNames',{'sessID','taskRow','bestLabel','cohort','group','participant','sessNum', ...
            'coup_MI','coup_prefPhaseRad','coup_resultantLen','coup_rayleighP','coup_inhExhRatio'});
        writetable(crow, outcp);
    end
    meta = struct('sessID',id,'taskRow',taskRow,'bestLabel',bl, ...
        'nBreaths',out.nBreaths,'tMs',out.tMs,'freqs',out.freqs, ...
        'bandFreqs',out.bandFreqs,'nT',out.nT);
    if isfield(out,'cohort'), meta.cohort=out.cohort; meta.group=out.group; ...
        meta.participant=out.participant; meta.sessNum=out.sessNum; end
    sumZwhole=out.sumZwhole; sumRawMap=out.sumRawMap; poolSumsqRaw=out.poolSumsqRaw; %#ok<NASGU>
    save(outmat,'sumZwhole','sumRawMap','poolSumsqRaw','meta','-v7');
catch e
    fprintf('  FAIL %s %s: %s\n', id, taskRow, e.message);
end
end

function bl = lookup_best(best, id, task)
bl='';
m = strcmp(best.sessID,id) & strcmp(best.task,task);
if any(m)
    v = best.bestLabel(find(m,1));
    bl = char(string(v));
end
end

function p = resolve_path(fp, rootOverride)
p = fp;
if exist(p,'file')==2, return; end
if ~isempty(rootOverride)
    % rebase '...\Lab_Common\...' onto rootOverride
    k = regexpi(fp,'Lab_Common','once');
    if ~isempty(k)
        tail = fp(k+numel('Lab_Common'):end);
        cand = [rootOverride tail];
        if exist(cand,'file')==2, p=cand; return; end
    end
end
end
