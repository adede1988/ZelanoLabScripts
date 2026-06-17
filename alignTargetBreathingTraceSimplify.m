function [outDat, targTraces] = alignTargetBreathingTraceSimplify(outDat, targTraceDir)
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
    nSamples   = diff([outDat.TTL size(outDat.data,2)]); 

    targTraces = cell(length(nSamples), 1); 

    % nSamples   = outDat.fs * 300;  % 300 s window at fs
    % targTraces = zeros(nCond, nSamples);

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

                targTraces{startFrom} = zeros(1, nSamples(startFrom)); 
                startFrom = startFrom + 1; 
            else
                searchOn = false; 
            end
        end
        startFrom = outDat.behDat.order(orderIdx(startFrom)); 
    else
        startFrom = 3; 
        targTraces{1} = zeros(1, nSamples(1)); 
        targTraces{2} = zeros(1, nSamples(2)); 
    end
    
    

    % Get the target files for the shadow conditions:
    for cndi = startFrom:nCond
        

        idx = find(tmpBehDat.order == cndi, 1);

        if isempty(idx)
            % No rows for this condition in tmpBehDat; skip
            error('unexpected missing condition!')
            continue;
        end
        if ismember(tmpBehDat.task(idx), {'fastFocus', 'slowFocus', 'naturalFocus'})

                targTraces{cndi} = zeros(1, nSamples(cndi)); 
            continue
        end
        if ~ismember(tmpBehDat.shadowFile{idx}, {'audioResp', 'focusedResp'})
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
            
            if cndi == length(outDat.TTL)
                targ_len = length(rspDat) - outDat.TTL(end); 
                seg = rspDat(outDat.TTL(cndi): end);
            else
                targ_len = outDat.TTL(cndi+1) - outDat.TTL(cndi); 
                seg = rspDat(outDat.TTL(cndi)+1: outDat.TTL(cndi+1));
            end
            L_goal   = length(seg); 
            timGoal  = 1/outDat.fs:1/outDat.fs:L_goal/outDat.fs; 
            L        = length(targTbl.voltage);
            timRec   = linspace(1/L, timGoal(end), L);
            target_resampled = interp1(timRec, ...
                                         targTbl.target, ...
                                         timGoal, 'linear');
            target_resampled = (target_resampled - ...
                                  mean(target_resampled, 'omitnan')) ./...
                                        std(target_resampled, [], 'omitnan'); 
            voltages_resampled = interp1(timRec, ...
                                         targTbl.voltage, ...
                                         timGoal, 'linear');
            voltages_resampled = (voltages_resampled - ...
                                  mean(voltages_resampled, 'omitnan')) ./...
                                        std(voltages_resampled, [], 'omitnan'); 
               % z-score each trace
            z_neuralynx = (seg - mean(seg, 'omitnan')) ./ std(seg, 0, 'omitnan');
            z_target    = (target_resampled - mean(target_resampled, 'omitnan')) ./ std(target_resampled, 0, 'omitnan');
            z_psychopy  = (voltages_resampled - mean(voltages_resampled, 'omitnan')) ./ std(voltages_resampled, 0, 'omitnan');
            
            figure('visible', true, 'position', [0,0,1000,500]);
            hold on
            
            % better separated colors
            plot(timGoal, z_neuralynx, 'LineWidth', 1.8, 'Color', [0.00 0.45 0.74]); % blue
            plot(timGoal, z_target,    'LineWidth', 1.8, 'Color', [0.85 0.33 0.10]); % orange
            plot(timGoal, z_psychopy,  'LineWidth', 1.8, 'Color', [0.20 0.60 0.20]); % green
            
            ylabel('Z-scored amplitude');
            xlim([30 70])
            xlabel('Time (s)');
            legend({'neuralynx', 'target', 'psychopy'}, 'Location', 'best')
            
            title(sprintf('Condition %d: Respiration vs Target Trace (%s : %s)', ...
                          cndi, sessionID, tmpBehDat.task{idx}), ...
                  'Interpreter','none');
            box off
            set(gca, 'LineWidth', 1.2, 'FontSize', 12)
            saveas(gcf,fullfile(outDat.figs, ['shadowResp' num2str(cndi) '.jpg']));
            targTraces{cndi} = target_resampled;
        else
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
                tempo_scale = str2num(tmpBehDat.warp{idx}); 
            catch
                tempo_scale = tmpBehDat.warp(idx);
            end
    
            % Target length in "stim samples" (before mapping to ephys fs)
            targ_len = round(nSamples(cndi) /outDat.fs * tmpBehDat.FPS(idx) * 2);
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
            voltages(isnan(voltages)) = mean(voltages, 'omitnan'); 


             if cndi == length(outDat.TTL)
                seg = rspDat(outDat.TTL(cndi)+1: end);
             else
                seg = rspDat(outDat.TTL(cndi)+1: outDat.TTL(cndi+1));
            end


            z_neuralynx = (seg - mean(seg, 'omitnan')) ./ std(seg, 0, 'omitnan');
            z_target    = (voltages - mean(voltages, 'omitnan')) ./ std(voltages, 0, 'omitnan');



            targTraces{cndi} = voltages;
            
            
            
            figure('visible', true, 'position', [0,0,1000,500]);
            hold on
            if length(z_neuralynx) < length(targTime)
                z_neuralynx = [z_neuralynx(:)' zeros(length(targTime)-length(z_neuralynx), 1)'];
                z_target = [z_target(:)' zeros(length(targTime)-length(z_target), 1)'];
            end
            if length(z_neuralynx) > length(targTime)
                z_neuralynx = z_neuralynx(1:length(targTime));
                z_target = z_target(1:length(targTime));
            end
            % better separated colors
            plot(targTime, z_neuralynx, 'LineWidth', 1.8, 'Color', [0.00 0.45 0.74]); % blue
            plot(targTime, z_target,    'LineWidth', 1.8, 'Color', [0.85 0.33 0.10]); % orange
            
            ylabel('Z-scored amplitude');
            xlim([30 70])
            xlabel('Time (s)');
            legend({'neuralynx', 'target', 'psychopy'}, 'Location', 'best')
            
            title(sprintf('Condition %d: Respiration vs Target Trace (%s : %s)', ...
                          cndi, sessionID, tmpBehDat.task{idx}), ...
                  'Interpreter','none');
            box off
            set(gca, 'LineWidth', 1.2, 'FontSize', 12)
            saveas(gcf,fullfile(outDat.figs, ['shadowResp' num2str(cndi) '.jpg']));
        end

    

        % figure
        % hold on 
        % plot(timGoal, seg)
        % plot(timGoal,target_resampled.*100)
        % plot(timGoal,voltages_resampled.*100)
        % targTraces{cndi} = target_resampled;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%start here






        % 
        % 
        % 
       
    end
    
    % Time x condition
    targTraces = cat(2, targTraces{:});

    % Append flattened target trace into outDat
    try
        outDat.data(end+1, :) = targTraces(:);
    catch
        if size(outDat.data, 2) > length(targTraces)
            outDat.data(end+1, :) = ...
                [targTraces(:)' ...
                zeros(size(outDat.data, 2)-length(targTraces), 1)'];
        else

            outDat.data(end+1, :) = targTraces(1:size(outDat.data, 2)); 
        end
        warning('wrong length target trace problem!')
    end
    outDat.labels{end+1}  = 'targTrace';

     %% Per-condition plots: respiration vs target trace
    % One figure per condition (cndi = 3:nCond → nCond-2 figures)
    % if ~isempty(rspDat)
    % 
    % 
    % 
    %     for cndi = 3:nCond
    %         % Segment indices for this condition in the full-session rspDat
    %         startIdx = outDat.TTL(cndi)+1; 
    %         endIdx   = outDat.TTL(cndi+1);
    % 
    %         if startIdx > numel(rspDat)
    %             warning('alignTargetBreathingTrace:RespTooShort', ...
    %                     'rspDat too short for condition %d (startIdx=%d). Skipping plot.', ...
    %                     cndi, startIdx);
    %             continue;
    %         end
    % 
    %         endIdx = min(endIdx, numel(rspDat));
    %         segLen = endIdx - startIdx + 1;
    % 
    %         segRsp  = rspDat(startIdx:endIdx);
    %         segTarg = targTraces(startIdx:endIdx);
    % 
    %         t = (0:segLen-1) / outDat.fs;
    % 
    %         figure('visible', false, 'position', [0,0,1000,500]);
    %         yyaxis left;
    %         plot(t, segRsp);
    %         ylabel('Respiration');
    % 
    %         yyaxis right;
    %         plot(t, segTarg);
    %         ylabel('Target trace (a.u.)');
    %         xlim([100 130])
    %         xlabel('Time (s)');
    %         title(sprintf('Condition %d: Respiration vs Target Trace (%s)', ...
    %                       cndi, sessionID), ...
    %               'Interpreter','none');
    %         saveas(gcf,fullfile(outDat.figs, ['shadowResp' num2str(cndi) '.jpg']));
    %     end
    % end

end
