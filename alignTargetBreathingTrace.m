function [outDat, targTraces] = alignTargetBreathingTrace(outDat, targTraceDir)
% alignTargetBreathingTrace
%   Build and align target breathing traces for the shadow conditions,
%   then append them into outDat as a flattened channel 'targTrace'.
%
% Inputs
%   outDat    : struct with .fs, .data, .labels, .behDat
%   codePre   : base path to code repo (prefix for 'closed-loop-respiration\data')
%
% Outputs
%   outDat     : updated struct with new channel 'targTrace'
%   targTraces : [time x condition] matrix of aligned target traces
    
    sessionID = outDat.sessID; 
    tmpBehDat = outDat.behDat; 
    %extract breathing data: 
    idx = cellfun(@(x) contains(x, 'rsp'), outDat.labels);
    rspDat = outDat.data(idx,:); 
    rspDat = rspDat(outDat.rspIDX,:);
    rspDat = rspDat .* outDat.rspFlip;

    % Number of conditions (rows) and samples (columns)
    nCond      = max(outDat.behDat.order);
    nSamples   = outDat.fs * 300;  % 300 s window at fs
    targTraces = zeros(nCond, nSamples);

    dataDir = fullfile(targTraceDir);
    if ismember('cndName', outDat.behDat.Properties.VariableNames)
        orderIdx = arrayfun(@(x) find(outDat.behDat.order == x, 1),...
                                unique(outDat.behDat.order));
        OrderCndNames = outDat.behDat.cndName(orderIdx);
        searchOn = true; 
        startFrom = 1; 
        while searchOn
            if strcmp(OrderCndNames{startFrom}, 'pre') || ...
               strcmp(OrderCndNames{startFrom}, 'audio')
                startFrom = startFrom + 1; 
            else
                searchOn = false; 
            end
        end
    else
        startFrom = 3; 
    end

    % Get the target files for the shadow conditions:
    for cndi = startFrom:nCond
        idx = find(tmpBehDat.order == cndi, 1);

        if isempty(idx)
            % No rows for this condition in tmpBehDat; skip
            continue;
        end

        % Determine which shadow file to use
        sfName = tmpBehDat.shadowFile{idx};
        if strcmp(sfName, 'NA')
            sfName = 'audioResp';
        end

        % Build full CSV path
        try
            csvName = sprintf('%s_%s_recording.csv', sessionID, sfName);
            csvPath = fullfile(dataDir, csvName);
    
            targTbl = readtable(csvPath);
        catch
            csvName = sprintf('%s%s_recording.csv', sessionID, sfName);
            csvPath = fullfile(dataDir, csvName);
    
            targTbl = readtable(csvPath);

        end

        % Remove any zero-voltage rows
        targTbl(targTbl.voltage == 0, :) = [];

        % Tempo scale can be a string or numeric
        try
            tempo_scale = str2num(tmpBehDat.warp{idx}); %#ok<ST2NM>
        catch
            tempo_scale = tmpBehDat.warp(idx);
        end

        % Target length in "stim samples" (before mapping to ephys fs)
        targ_len = round(tmpBehDat.trialTim(idx) * tmpBehDat.FPS(idx) * 2);
        new_len  = round(length(targTbl.voltage) / tempo_scale);

        % Resample original voltage trace to tempo-scaled length
        L        = length(targTbl.voltage);
        timRec   = linspace(1/L, 1, L);
        timGoal  = linspace(1/new_len, 1, new_len);
        voltages_resampled = interp1(timRec, ...
                                     targTbl.voltage, ...
                                     timGoal, 'linear');

        % Take loop segment, trimming 180 s (scaled) from each end
        loop_start   = round(180 / tempo_scale);
        loop_end     = length(voltages_resampled) - round(180 / tempo_scale);
        loop_segment = voltages_resampled(loop_start:loop_end);
        loopLen      = length(loop_segment);

        % Match loop length to targ_len by truncation or repetition
        if loopLen > targ_len
            voltages = loop_segment(1:targ_len);
        elseif loopLen < targ_len
            repeats  = ceil(targ_len / loopLen);
            voltages = repmat(loop_segment, 1, repeats);
            voltages = voltages(1:targ_len);
        else
            voltages = loop_segment;
        end

        % Cut to trialTim length using recording timestamps
        timStp = mean(diff(targTbl.timestamp));
        tmpTim = timStp:timStp:tmpBehDat.trialTim(idx);
        voltages = voltages(1:length(tmpTim));

        % Resample to match ephys data (fs) over 300 s
        targTime = 1/outDat.fs : 1/outDat.fs : 300;
        voltages = interp1(tmpTim, voltages, targTime, 'linear');

        targTraces(cndi, :) = voltages;
    end

    % Time x condition
    targTraces = targTraces';

    % Append flattened target trace into outDat
    try
        outDat.data(end+1, :) = targTraces(:);
    catch
        outDat.data(end+1, 1:length(targTraces(:))) = targTraces(:); 
    end
    outDat.labels{end+1}  = 'targTrace';

     %% Per-condition plots: respiration vs target trace
    % One figure per condition (cndi = 3:nCond → nCond-2 figures)
    if ~isempty(rspDat)
        nSamples  = outDat.fs * 300;       % already used above
        [~, nCond] = size(targTraces);     % time x condition

        for cndi = 3:nCond
            % Segment indices for this condition in the full-session rspDat
            startIdx = (cndi-1) * nSamples + 1;
            endIdx   = cndi * nSamples;

            if startIdx > numel(rspDat)
                warning('alignTargetBreathingTrace:RespTooShort', ...
                        'rspDat too short for condition %d (startIdx=%d). Skipping plot.', ...
                        cndi, startIdx);
                continue;
            end

            endIdx = min(endIdx, numel(rspDat));
            segLen = endIdx - startIdx + 1;

            segRsp  = rspDat(startIdx:endIdx);
            segTarg = targTraces(1:segLen, cndi);

            t = (0:segLen-1) / outDat.fs;

            figure('visible', false, 'position', [0,0,1000,500]);
            yyaxis left;
            plot(t, segRsp);
            ylabel('Respiration');

            yyaxis right;
            plot(t, segTarg);
            ylabel('Target trace (a.u.)');
            xlim([100 130])
            xlabel('Time (s)');
            title(sprintf('Condition %d: Respiration vs Target Trace (%s)', ...
                          cndi, sessionID), ...
                  'Interpreter','none');
            saveas(gcf,fullfile(outDat.figs, ['shadowResp' num2str(cndi) '.jpg']));
        end
    end

end
