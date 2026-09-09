function probe_structs()
% Probe the four task finals for one Dupi session: structure, labels, behDat cols, timing.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));      % ...\ZelanoLabScripts
addpath(repo);
try, L = labPaths(); catch e, fprintf('labPaths err: %s\n', e.message); L = struct('rootDupi','R:\Neurology\Zelano_Lab\Lab_Common\Dupi\'); end

base = 'R:\Neurology\Zelano_Lab\Lab_Common\Dupi\250623_Dupi_NMH_KS_1\preProc\';
files = struct( ...
  'breathing','250623_Dupi_NMH_KS_1_breathingPreProc.mat', ...
  'cue',      '250623_Dupi_NMH_KS_1_cueTaskPreproc.mat', ...
  'O15',      '250623_Dupi_NMH_KS_1_O15preproc.mat', ...
  'thresh',   '250623_Dupi_NMH_KS_1_PEA_threshold_preproc.mat');

tasks = fieldnames(files);
for ti = 1:numel(tasks)
    tk = tasks{ti};
    fp = fullfile(base, files.(tk));
    fprintf('\n================ %s ================\n%s\n', tk, fp);
    if exist(fp,'file')~=2, fprintf('  MISSING\n'); continue; end
    d = dir(fp); fprintf('  size = %.1f MB\n', d.bytes/1e6);
    t0 = tic;
    S = load(fp);
    fn = fieldnames(S);
    od = S.(fn{1});
    fprintf('  load time = %.1f s ; topVar = %s\n', toc(t0), fn{1});
    fprintf('  outDat fields: %s\n', strjoin(fieldnames(od)', ', '));
    if isfield(od,'data'),   fprintf('  data size = [%d x %d]\n', size(od.data,1), size(od.data,2)); end
    if isfield(od,'fs'),     fprintf('  fs = %g\n', od.fs); end
    if isfield(od,'task'),   fprintf('  task = %s\n', char(string(od.task))); end
    if isfield(od,'type'),   fprintf('  type = %s\n', char(string(od.type))); end
    if isfield(od,'sessID'), fprintf('  sessID = %s\n', char(string(od.sessID))); end
    if isfield(od,'moreThan1'), fprintf('  moreThan1 = %g\n', od.moreThan1); end
    if isfield(od,'bestMac'), fprintf('  bestMac(cached) = %s\n', char(string(od.bestMac))); end
    if isfield(od,'labels')
        labs = od.labels;
        isMac = cellfun(@(x) contains(lower(char(string(x))),'macbp'), labs);
        fprintf('  nLabels = %d ; macBP labels = %s\n', numel(labs), strjoin(cellfun(@(x)char(string(x)),labs(isMac),'uni',0),', '));
        isRsp = cellfun(@(x) contains(lower(char(string(x))),'rsp'), labs);
        fprintf('  rsp labels = %s ; rspIDX=%g rspFlip=%g\n', strjoin(cellfun(@(x)char(string(x)),labs(isRsp),'uni',0),', '), ...
            getfielddef(od,'rspIDX',NaN), getfielddef(od,'rspFlip',NaN));
    end
    if isfield(od,'behDat')
        bd = od.behDat;
        if istable(bd)
            fprintf('  behDat is TABLE [%d x %d]; vars: %s\n', height(bd), width(bd), strjoin(bd.Properties.VariableNames, ', '));
        else
            fprintf('  behDat is %s size [%s]\n', class(bd), num2str(size(bd)));
        end
    end
    if isfield(od,'bmObj'), fprintf('  bmObj size = [%d x %d]\n', size(od.bmObj,1), size(od.bmObj,2)); end
    if isfield(od,'TTL')
        if istable(od.TTL), fprintf('  TTL table [%d x %d]: %s\n', height(od.TTL), width(od.TTL), strjoin(od.TTL.Properties.VariableNames,', '));
        else, fprintf('  TTL %s size [%s]\n', class(od.TTL), num2str(size(od.TTL))); end
    end
    clear S od;
end
fprintf('\nDONE probe\n');
end

function v = getfielddef(s,f,d)
    if isfield(s,f), v = s.(f); else, v = d; end
    if ~isnumeric(v)||~isscalar(v), v = NaN; end
end
