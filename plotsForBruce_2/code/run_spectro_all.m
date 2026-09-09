function run_spectro_all(useParfor, rootOverride)
% RUN_SPECTRO_ALL  Batch group-spectrogram sufficient stats (3 freq ranges,
% -500..-100 ms baseline, respiration overlay) over all Dupi+OBE recordings x
% 5 task-rows. Writes out/spectro2/<sessID>__<taskRow>.mat.
if nargin<1||isempty(useParfor), useParfor=false; end
if nargin<2, rootOverride=''; end
here=fileparts(mfilename('fullpath')); proj=fileparts(here); repo=fileparts(proj);
addpath(repo); addpath(here);
for c={repo,'C:\Users\Adam\Documents\GitHub\ZelanoLabScripts','E:\GitHub\ZelanoLabScripts','G:\My Drive\GitHub\ZelanoLabScripts'}
    if exist(fullfile(c{1},'fooof_basic.m'),'file'), addpath(c{1}); end, end
for c={'C:\Users\Adam\Documents\GitHub\Superlets\matlab-pure','E:\GitHub\Superlets\matlab-pure','G:\My Drive\GitHub\Superlets\matlab-pure'}
    if exist(c{1},'dir'), addpath(c{1}); end, end
assert(exist('slt_power_cont','file')==2,'superlet missing');
tdir=fullfile(proj,'out','tables'); sdir=fullfile(proj,'out','spectro2');
if ~exist(sdir,'dir'), mkdir(sdir); end
idx=readtable(fullfile(tdir,'session_index.csv'),'TextType','string');
best=readtable(fullfile(tdir,'macbp_best.csv'),'TextType','string');
idx=idx((idx.cohort=="Dupi"|idx.cohort=="OBE") & idx.onDisk==1,:);
U={};
for r=1:height(idx)
    task=char(idx.task(r)); id=char(idx.sessID(r)); fp=resolve_path(char(idx.finalPath(r)),rootOverride);
    m=strcmp(best.sessID,id)&strcmp(best.task,task); if ~any(m), continue; end
    bl=char(string(best.bestLabel(find(m,1))));
    switch task
        case 'cueTask', U(end+1,:)={fp,'cueTask',bl,id};
        case 'threshTask', U(end+1,:)={fp,'threshTask',bl,id};
        case 'O15', U(end+1,:)={fp,'O15',bl,id};
        case 'breathingTask'
            U(end+1,:)={fp,'audiobook',bl,id}; U(end+1,:)={fp,'focusedBreathing',bl,id};
    end
end
nU=size(U,1); fprintf('spectro batch: %d units parfor=%d\n', nU, useParfor);
if useParfor
    nw=3; delete(gcp('nocreate')); parpool('Processes',nw);   % 3 workers: avoid parallel-FFT memory-bandwidth wall
    fprintf('parpool up: %d workers\n', nw);
    parfor i=1:nU, proc(U(i,:), sdir); end
    delete(gcp('nocreate'));
else
    for i=1:nU, fprintf('[%3d/%3d] %s %s\n',i,nU,U{i,4},U{i,2}); proc(U(i,:), sdir); end
end
fprintf('DONE run_spectro_all\n');
end

function proc(u, sdir)
fp=u{1}; taskRow=u{2}; bl=u{3}; id=u{4};
outm=fullfile(sdir,sprintf('%s__%s.mat',id,taskRow)); if exist(outm,'file'), return; end
try
    out=extract_spectro_session(fp,taskRow,bl,struct());
    save(outm,'out','-v7');
catch e
    fprintf('  FAIL %s %s: %s\n', id, taskRow, e.message);
end
end
function p=resolve_path(fp,ro)
p=fp; if exist(p,'file')==2, return; end
if ~isempty(ro), k=regexpi(fp,'Lab_Common','once'); if ~isempty(k), c=[ro fp(k+numel('Lab_Common'):end)]; if exist(c,'file')==2, p=c; end, end, end
end
