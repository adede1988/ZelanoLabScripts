function [idx, cfg] = reproc_curatedIdx(task)
% Session indices (into applyParams(task,'main')) to reprocess: everything whose
% paramSource is not 'guess' (curated + blank-treated-as-curated).
    cfg = applyParams(task, 'main');
    idx = find(~strcmpi(cfg.paramSource, 'guess'));
end
