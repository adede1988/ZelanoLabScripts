function diag_thresh()
% What scale is the intensity rating on, and is it stable across sessions?
cfg  = ohrv_config();
root = cfg.dataRoot;
subs = ["AB","DB","GH","JH","JL","KS","PC"];
d = dir(fullfile(root,'*','preProc','*PEA_threshold*.mat'));
fprintf('%-26s %5s | %-28s | %-28s\n','session','nTr','intensity by odor (mean)','pleasantness range');
fprintf('%-26s %5s | %8s %8s %8s | %8s %8s\n','','','air','odor2','odor3','min','max');
for ii = 1:numel(d)
    nm = d(ii).name; tok = split(string(nm),'_');
    if numel(tok) < 4 || ~ismember(tok(4), subs), continue; end
    try
        S = load(fullfile(d(ii).folder, nm));
        fn = fieldnames(S); pick = fn{1};
        for k=1:numel(fn), if ismember(fn{k},{'outDat','chanDat','out'}), pick=fn{k}; break; end, end
        b = S.(pick).behDat; clear S
        if ~all(ismember({'odor','intensity'}, b.Properties.VariableNames)), continue; end
        if ismember('n', b.Properties.VariableNames)
            [~,ia] = unique(double(b.n),'stable'); b = b(ia,:);
        end
        od = double(b.odor); it = double(b.intensity);
        m = @(v) mean(it(od==v & isfinite(it)),'omitnan');
        pl = NaN(1,2);
        if ismember('pleasantness', b.Properties.VariableNames)
            p = double(b.pleasantness); pl = [min(p) max(p)];
        end
        base = erase(string(nm), "_PEA_threshold_preproc.mat");
        fprintf('%-26s %5d | %8.1f %8.1f %8.1f | %8.1f %8.1f   [it %.0f..%.0f]\n', ...
            base, height(b), m(1), m(2), m(3), pl(1), pl(2), min(it), max(it));
    catch ME
        fprintf('%-26s  ERROR %s\n', nm, ME.message(1:min(40,end)));
    end
end
end
