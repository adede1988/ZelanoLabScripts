function run_o15_task_group(groupDir)
% RUN_O15_TASK_GROUP  Group-mean baseline-z spectrograms (+ respiration overlay)
%   for the O15 task, from the per-subject z-TFRs. Averages the per-subject maps
%   within each group (DupiS1/S2/S3/Control) for each of the THREE sniff-type
%   lockings (start/free/confirm), on one shared (data-driven) color scale, and
%   overlays the group-mean respiration. Mirrors run_cue_task4_group.

    repo = fileparts(fileparts(mfilename('fullpath')));
    addpath(repo); addpath(fullfile(repo,'cueAnalysis')); addpath(fullfile(repo,'O15Analysis'));
    L = labPaths();
    if nargin < 1 || isempty(groupDir)
        groupDir = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\groupStatFigs';
    end
    if ~isfolder(groupDir), mkdir(groupDir); end

    dispWin = [-1000 3000];
    groups = {'DupiS1','DupiS2','DupiS3','Control'};
    lockings = {'start','free','confirm'};

    T = o15_session_table(false); T = T(T.onDisk & T.fresh,:);

    G = struct();
    for gi = 1:numel(groups)
        for li = 1:numel(lockings)
            G.(groups{gi}).(lockings{li}) = struct('maps',{{}},'resp',{{}}, ...
                'times',[],'freqs',[],'respT',[],'ids',{{}});
        end
    end

    for i = 1:height(T)
        id = char(T.sessID(i)); grp = char(T.group(i));
        if ~ismember(grp, groups), continue; end
        f = fullfile(L.figPath, id, 'O15', [id '_O15_bestMac_TFR.mat']);
        if ~isfile(f), continue; end
        Sv = load(f); to = Sv.tfrOut;
        if ~isfield(to,'freqs'), continue; end
        for li = 1:numel(lockings)
            lk = lockings{li};
            if isfield(to,lk) && isstruct(to.(lk)) && ~isempty(to.(lk)) && isfield(to.(lk),'map')
                S = G.(grp).(lk);
                S.maps{end+1} = to.(lk).map;
                S.resp{end+1} = to.(lk).resp;
                S.times = to.(lk).times; S.freqs = to.freqs; S.respT = to.(lk).respT;
                S.ids{end+1} = id;
                G.(grp).(lk) = S;
            end
        end
    end

    % ---- means + data-driven shared color scale across all group maps ----
    allvals = []; M = struct();
    for gi = 1:numel(groups)
        for li = 1:numel(lockings)
            S = G.(groups{gi}).(lockings{li});
            if isempty(S.maps), M.(groups{gi}).(lockings{li}) = []; continue; end
            mp = mean(cat(3, S.maps{:}), 3, 'omitnan');
            rp = mean(cat(1, S.resp{:}), 1, 'omitnan');
            M.(groups{gi}).(lockings{li}) = struct('map',mp,'resp',rp, ...
                'times',S.times,'freqs',S.freqs,'respT',S.respT,'n',numel(S.maps));
            allvals = [allvals; abs(mp(:))]; %#ok<AGROW>
        end
    end
    a = prctile(allvals, 99); if isempty(a) || ~isfinite(a) || a <= 0, a = 3; end
    clim = [-a a];

    for gi = 1:numel(groups)
        g = groups{gi};
        for li = 1:numel(lockings)
            lk = lockings{li};
            mm = M.(g).(lk);
            if isempty(mm), fprintf('%-8s %-8s : (none)\n', g, lk); continue; end
            cue_plot_ztfr(mm.map, mm.times, mm.freqs, clim, mm.resp, mm.respT, dispWin, ...
                sprintf('GROUP %s  %s-sniff  (z, n=%d)', g, lk, mm.n), ...
                fullfile(groupDir, ['group_' g '_O15TFR_' lk '.png']));
            fprintf('%-8s %-8s : n=%d\n', g, lk, mm.n);
        end
    end

    save(fullfile(groupDir, 'O15Task_group_means.mat'), 'M', 'clim', 'groups', 'lockings', 'dispWin', '-v7');
    fprintf('Saved O15 group means + figures to %s (clim=[%.2f %.2f])\n', groupDir, clim(1), clim(2));

    o15_make_manifest();
end
