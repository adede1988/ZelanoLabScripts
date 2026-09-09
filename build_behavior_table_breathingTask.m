function outDat = build_behavior_table_breathingTask(outDat, bmObj)
% constructBehDat
%   Rebuilds outDat.behDat from breath-mark matrix bmObj and the
%   original behavioral table in outDat.behDat, and also constructs
%   outDat.baseEmotion.
%
% Inputs
%   outDat : struct with fields
%              .fs      : sampling rate
%              .data    : [channels x time]
%              .behDat  : original behavior table (tmpBehDat)
%   bmObj  : [nBreaths x 14+] matrix of breath markers:
%            col  1: onset Y value
%            col  2: onset time (s)
%            col  3: peak Y value
%            col  4: peak time (s)
%            col  5: end Y value
%            col  6: end time (s)
%            col  7: length (end - onset)
%            col  8: amplitude
%            col  9: idx of peak in rspSig2
%            col 10: exhale peak Y
%            col 11: exhale peak time (s)
%            col 12: condition
%            col 14: index
%
% Output
%   outDat : same struct with updated
%              .behDat      : per-breath table
%              .baseEmotion : baseline emotion row

    % Preserve original behavior table for emotion info
    tmpBehDat = outDat.behDat;

    % New behDat table based on bmObj
    outDat.behDat = table();

    % Build time axis in seconds for index conversion
    tim = (1:size(outDat.data, 2)) / outDat.fs;

    % Column 2: onset time → sniffOnset and finalOnset (sample indices)
    idx = arrayfun(@(x) find(x <= tim, 1), bmObj(:, 2));
    outDat.behDat.sniffOnset = idx;
    outDat.behDat.finalOnset = idx;

    % Manual-onset placeholder column: NaN here, filled in by hand during later QC.
    outDat.behDat.manOnset = nan(size(idx));

    % Column 12: condition
    outDat.behDat.condition = bmObj(:, 12);

    % Column 1: onset Y value
    outDat.behDat.Yonset = bmObj(:, 1);

    % Column 3: peak Y value
    outDat.behDat.inhaleMax = bmObj(:, 3);

    % Column 4: peak time → inMaxTim (sample indices)
    idx = arrayfun(@(x) find(x <= tim, 1), bmObj(:, 4));
    outDat.behDat.inMaxTim = idx;

    % Column 5: end Y value
    outDat.behDat.Yend = bmObj(:, 5);

    % Column 6: end time → endTim (sample indices)
    idx = arrayfun(@(x) find(x <= tim, 1), bmObj(:, 6));
    outDat.behDat.endTim = idx;

    % Column 7: length (end - onset)
    outDat.behDat.length = bmObj(:, 7);

    % Column 8: amp (peak Y - avg of two ends)
    outDat.behDat.amp = bmObj(:, 8);

    % Column 10: exhale peak Y value
    outDat.behDat.exhaleMin = bmObj(:, 10);

    % Column 11: exhale peak time → exMinTim (sample indices)
    idx = arrayfun(@(x) find(x <= tim, 1), bmObj(:, 11));
    outDat.behDat.exMinTim = idx;

    % Column 14: index
    outDat.behDat.index = bmObj(:, 14);

    % ---------------- integrate emotion data into respiration ----------------
    Qs = unique(tmpBehDat.Q_short);  % e.g., emotion questions

    for cndi = 1:max(outDat.behDat.condition)
        idx = find(outDat.behDat.condition == cndi);
        if isempty(idx)
            continue;
        end

        tmp = tmpBehDat(tmpBehDat.order == cndi, :);
        if isempty(tmp)
            continue;
        end

        outDat.behDat.task(idx)       = tmp.task(1);
        outDat.behDat.noseMouth(idx)  = tmp.noseMouth(1);
        outDat.behDat.shadowFile(idx) = tmp.shadowFile(1);
        outDat.behDat.warp(idx)       = tmp.warp(1);

        % Add one column per Q/type combination: e.g., "Happy_pre" etc.
        for q = 1:length(Qs)
            ii = find(cellfun(@(x) strcmp(Qs{q}, x), tmp.Q_short));
            if isempty(ii)
                continue;
            end
            varName = [Qs{q} '_' tmp.type{ii}];
            outDat.behDat.(varName)(idx) = tmp.rsp(ii);
        end
    end

    % ---------------- baseline emotion row (order == 0) ----------------
    cndi = 0;
    baseEmotion = table;
    tmp = tmpBehDat(tmpBehDat.order == cndi, :);

    if ~isempty(tmp)
        baseEmotion.task       = tmp.task(1);
        baseEmotion.noseMouth  = tmp.noseMouth(1);
        baseEmotion.shadowFile = tmp.shadowFile(1);
        baseEmotion.warp       = tmp.warp(1);

        for q = 1:length(Qs)
            ii = find(cellfun(@(x) strcmp(Qs{q}, x), tmp.Q_short));
            if isempty(ii)
                continue;
            end
            varName = [Qs{q} '_' tmp.type{ii}];
            baseEmotion.(varName) = tmp.rsp(ii);
        end
    end

    outDat.baseEmotion = baseEmotion;

end
