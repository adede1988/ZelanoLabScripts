function run_cue_task4_group(groupDir)
% RUN_CUE_TASK4_GROUP  Group-mean z-score spectrograms (+ respiration overlay)
%   from the per-subject z-TFRs (analysis1 format). Averages the per-subject
%   bootstrap-z maps within each group (DupiS1/S2/S3/Control), on one shared
%   color scale, and overlays the group-mean respiration trace.

    repo = fileparts(fileparts(mfilename('fullpath')));
    addpath(repo); addpath(fullfile(repo,'cueAnalysis'));
    L = labPaths();
    if nargin < 1 || isempty(groupDir)
        groupDir = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\groupStatFigs';
    end
    if ~isfolder(groupDir), mkdir(groupDir); end

    dispWin = [-1000 3000];
    groups = {'DupiS1','DupiS2','DupiS3','Control'};
    lockings = {'trialStart','finalOnset'};

    T = cue_session_table(false); T = T(T.onDisk,:);

    % collect per-subject maps + resp traces in cells, per group x locking
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
        f = fullfile(L.figPath, id, 'cueTask', [id '_cue_bestMac_TFR.mat']);
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

    % ---- means + shared color scale across all group maps ----
    allvals = [];
    M = struct();
    for gi = 1:numel(groups)
        for li = 1:numel(lockings)
            S = G.(groups{gi}).(lockings{li});
            if isempty(S.maps)
                M.(groups{gi}).(lockings{li}) = []; continue;
            end
            mp = mean(cat(3, S.maps{:}), 3, 'omitnan');
            rp = mean(cat(1, S.resp{:}), 1, 'omitnan');
            M.(groups{gi}).(lockings{li}) = struct('map',mp,'resp',rp, ...
                'times',S.times,'freqs',S.freqs,'respT',S.respT,'n',numel(S.maps));
            allvals = [allvals; abs(mp(:))]; %#ok<AGROW>
        end
    end
    % Fixed color axis (+-5 z) so smaller-but-meaningful effects are emphasized.
    clim = [-5 5]; %#ok<NASGU>  (allvals computed above is no longer used for clim)

    % ---- plots ----
    for gi = 1:numel(groups)
        g = groups{gi};
        for li = 1:numel(lockings)
            lk = lockings{li};
            mm = M.(g).(lk);
            if isempty(mm), fprintf('%-8s %-11s : (none)\n', g, lk); continue; end
            cue_plot_ztfr(mm.map, mm.times, mm.freqs, clim, mm.resp, mm.respT, dispWin, ...
                sprintf('GROUP %s  %s-locked  (z, n=%d)', g, lk, mm.n), ...
                fullfile(groupDir, ['group_' g '_TFR_' lk '.png']));
            fprintf('%-8s %-11s : n=%d\n', g, lk, mm.n);
        end
    end

    save(fullfile(groupDir, 'cueTask_group_means.mat'), 'M', 'clim', 'groups', 'lockings', 'dispWin', '-v7');
    fprintf('Saved group means + figures to %s (clim=[%.2f %.2f])\n', groupDir, clim(1), clim(2));

    % keep the data-availability manifest in sync with the live FOOOF+QC CSV
    cue_make_manifest();
end
