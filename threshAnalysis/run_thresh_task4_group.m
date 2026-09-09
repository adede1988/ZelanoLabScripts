function run_thresh_task4_group(groupDir)
% RUN_THRESH_TASK4_GROUP  Group-mean z-score spectrograms (+ respiration overlay)
%   from the per-subject condition-split z-TFRs. Averages the per-subject
%   bootstrap-z maps within each group (DupiS1/S2/S3/Control) and odor condition
%   (low/med/air), on one shared color scale, and overlays the group-mean
%   respiration trace. 4 groups x 3 conditions = up to 12 PNGs.

    repo = fileparts(fileparts(mfilename('fullpath')));
    addpath(repo); addpath(fullfile(repo,'threshAnalysis'));
    L = labPaths();
    if nargin < 1 || isempty(groupDir)
        groupDir = getenv('THRESH_GROUPDIR');
        if isempty(groupDir), groupDir = fullfile(L.figPath, 'groupStatFigs'); end
    end
    if ~isfolder(groupDir), mkdir(groupDir); end

    dispWin = [-1000 3000];
    groups = {'DupiS1','DupiS2','DupiS3','Control'};
    conds  = {'low','med','air'};

    T = thresh_session_table(false); T = T(T.onDisk,:);

    % collect per-subject maps + resp traces in cells, per group x condition
    G = struct();
    for gi = 1:numel(groups)
        for li = 1:numel(conds)
            G.(groups{gi}).(conds{li}) = struct('maps',{{}},'resp',{{}}, ...
                'times',[],'freqs',[],'respT',[],'ids',{{}});
        end
    end

    for i = 1:height(T)
        id = char(T.sessID(i)); grp = char(T.group(i));
        if ~ismember(grp, groups), continue; end
        f = fullfile(L.figPath, id, 'threshTask', [id '_thresh_bestMac_TFR.mat']);
        if ~isfile(f), continue; end
        Sv = load(f); to = Sv.tfrOut;
        if ~isfield(to,'freqs'), continue; end
        for li = 1:numel(conds)
            ck = conds{li};
            if isfield(to,ck) && isstruct(to.(ck)) && ~isempty(to.(ck)) && isfield(to.(ck),'map')
                S = G.(grp).(ck);
                S.maps{end+1} = to.(ck).map;
                S.resp{end+1} = to.(ck).resp;
                S.times = to.(ck).times; S.freqs = to.freqs; S.respT = to.(ck).respT;
                S.ids{end+1} = id;
                G.(grp).(ck) = S;
            end
        end
    end

    % ---- means (fixed shared color scale +-5 z) ----
    clim = [-5 5];
    M = struct();
    for gi = 1:numel(groups)
        for li = 1:numel(conds)
            S = G.(groups{gi}).(conds{li});
            if isempty(S.maps)
                M.(groups{gi}).(conds{li}) = []; continue;
            end
            mp = mean(cat(3, S.maps{:}), 3, 'omitnan');
            rp = mean(cat(1, S.resp{:}), 1, 'omitnan');
            M.(groups{gi}).(conds{li}) = struct('map',mp,'resp',rp, ...
                'times',S.times,'freqs',S.freqs,'respT',S.respT,'n',numel(S.maps));
        end
    end

    % ---- plots ----
    for gi = 1:numel(groups)
        g = groups{gi};
        for li = 1:numel(conds)
            ck = conds{li};
            mm = M.(g).(ck);
            if isempty(mm), fprintf('%-8s %-4s : (none)\n', g, ck); continue; end
            thresh_plot_ztfr(mm.map, mm.times, mm.freqs, clim, mm.resp, mm.respT, dispWin, ...
                sprintf('GROUP %s  %s-odor finalOnset-locked  (z, n=%d)', g, ck, mm.n), ...
                fullfile(groupDir, ['group_' g '_TFR_' ck '.png']));
            fprintf('%-8s %-4s : n=%d\n', g, ck, mm.n);
        end
    end

    save(fullfile(groupDir, 'threshTask_group_means.mat'), 'M', 'clim', 'groups', 'conds', 'dispWin', '-v7');
    fprintf('Saved group means + figures to %s (clim=[%.2f %.2f])\n', groupDir, clim(1), clim(2));

    % keep the data-availability manifest in sync with the live FOOOF+QC CSV
    thresh_make_manifest();
end
