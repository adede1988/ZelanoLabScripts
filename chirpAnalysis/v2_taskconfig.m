function tc = v2_taskconfig(task)
% V2_TASKCONFIG  Per-task config: final-file suffix, behDat category column, categories, O15 guard.
    switch task
        case 'cueTask',    tc = struct('suffix','cueTaskpreproc',       'catCol','type',      'cats',{{'hit','cr'}},            'guard',false);
        case 'threshTask', tc = struct('suffix','PEA_threshold_preproc','catCol','type',      'cats',{{'air','low','med'}},     'guard',false);
        case 'O15',        tc = struct('suffix','O15preproc',           'catCol','sniffLabel','cats',{{'start','free','confirm'}},'guard',true);
        otherwise, error('unknown task %s', task);
    end
end
