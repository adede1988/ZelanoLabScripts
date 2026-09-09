function v2_regroup(task)
% V2_REGROUP  Relabel existing per-subject aggregates (Dupi -> DupiS1/S2/S3; OBE stays) and
%   re-run the group step for a task, WITHOUT re-running the per-session batch. Cleans up the
%   orphaned old two-group "_Dupi" group figures. Assumes chirp paths are already on the path.
    C = v2_config();
    outDir = getenv('CHIRP_V2OUT'); if isempty(outDir), outDir = 'E:\chirpV2out'; end
    aggDir = fullfile(outDir, ['agg_' task]);
    d = dir(fullfile(aggDir,'*_agg.mat'));
    for k = 1:numel(d)
        S = load(fullfile(d(k).folder,d(k).name)); a = S.agg;
        a.group = v2_grouplabel(a.group, a.sessID);
        agg = a; save(fullfile(d(k).folder,d(k).name),'agg','-v7.3'); %#ok<NASGU>
    end
    gsf = fullfile(C.figRootE,'groupStatFigs');
    old = dir(fullfile(gsf, ['group*_' task '_*_Dupi.png']));   % orphaned 2-group Dupi figs
    for k = 1:numel(old), try, delete(fullfile(old(k).folder,old(k).name)); catch, end, end
    tc = v2_taskconfig(task);
    v2_group(aggDir, gsf, outDir, C, task, tc.cats);
    fprintf('v2_regroup(%s): %d agg relabeled; removed %d orphan figs; group refreshed\n', task, numel(d), numel(old));
end
