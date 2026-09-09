function [out, badChan, blinkIndicator] = blinkRemoveWrapper(outDat, ...
                                                            chanIDX, fs, interpChan)

    data = outDat.data(chanIDX, :); 
    origIs2D = ismatrix(data);
    if origIs2D
        data    = reshape(data,    size(data,1), size(data,2), 1);
        
    end
    [C,T,N] = size(data);
    data = reshape(data, C, []);
            
    %hard coded index of blink channel! 
    blinkChan = 1;      
    
    % data: [nChannels x nSamples], units = microvolts
    % outDat.fs: sampling rate (Hz)
    
    fs = outDat.fs;
    
    %% ---------- Criterion 1: high-voltage deviations in 2 s epochs ----------
    
    epochDurSec  = 2;                   % 2-second epochs
    epochSamples = round(epochDurSec * fs);
    nEpochs      = floor(T / epochSamples);
    
    % Trim to full epochs
    dataTrim   = data(:, 1:nEpochs*epochSamples);
    dataEpochs = reshape(dataTrim, C, epochSamples, nEpochs);  % [ch x samples x epochs]
    
    highRangeThresh = 700;   % µV: max - min >= 700 counts as a "bad" epoch
    epochFracThresh = 0.15;  % 15% of epochs
    
    % Range per channel per epoch
    epochMax   = squeeze(max(dataEpochs, [], 2));   % [nCh x nEpochs]
    epochMin   = squeeze(min(dataEpochs, [], 2));   % [nCh x nEpochs]
    epochRange = epochMax - epochMin;               % [nCh x nEpochs]
    
    % Epoch is "bad" if range >= 700 µV
    epochBad = epochRange >= highRangeThresh;       % logical [nCh x nEpochs]
    
    % Channel is bad if >= 15% of its epochs are bad
    badHighDev = (sum(epochBad, 2) ./ nEpochs) >= epochFracThresh;   % [nCh x 1]
    
    
    %% ---------- Criterion 2: flat periods (> 100 ms within ±1 µV) ----------
    
    flatAmpThresh      = 1;                % µV: between -1 and +1
    flatLenSec         = 0.100;            % 100 ms
    flatLenSamples     = round(flatLenSec * fs);
    flatTotalSecThresh = 30;               % 30 seconds total
    flatTotalSamples   = round(flatTotalSecThresh * fs);
    
    badFlat = false(C,1);
    
    for ch = 1:C
        x = data(ch,:);                         % 1 x nSamp
        isFlat = abs(x) <= flatAmpThresh;       % logical
    
        if ~any(isFlat)
            continue;                           % no flat segments at all
        end
    
        % Find contiguous flat segments via run-length encoding
        d = diff([0 isFlat 0]);                % edges
        starts = find(d == 1);                 % start indices of flat segments
        ends   = find(d == -1) - 1;            % end indices
    
        segLen = ends - starts + 1;            % length of each flat segment in samples
    
        % Only count segments that are at least 100 ms long
        longSegIdx = segLen >= flatLenSamples;
        totalFlatSamples = sum(segLen(longSegIdx));
    
        % Mark channel as bad if total flat time >= 30 s
        badFlat(ch) = totalFlatSamples >= flatTotalSamples;
    end
    
    
    %% ---------- Combine criteria ----------
    
    badChan = find(badHighDev | badFlat);% channels failing either criterion
    badChannelsHighDev = find(badHighDev);
    badChannelsFlat    = find(badFlat);
    
    if ismember(blinkChan, badChan)
        if ismember(C, badChan)
            error('no available blink chan!')
        else
            blinkChan = C; 
        end
    end


    % Optional: display
    % fprintf('Bad channels (any criterion): %s\n', mat2str(badChan));
    % fprintf('Bad channels (high deviation): %s\n', mat2str(badChannelsHighDev));
    % fprintf('Bad channels (flat): %s\n', mat2str(badChannelsFlat));


    % figure; 
    % imagesc(data) %check for bad channels overall 
    % caxis([-200, 200])
    % 
    % badChan = input(sprintf(...
    % 'Enter the index of badChans (1..%d), or [] to skip: '...
    %                                         ,size(data,1)));
    % 
    % figure; 
    % imagesc(data) %check for bad channels overall 
    % caxis([-1, 1])
    % 
    % badChan = [badChan input(sprintf(...
    % 'Enter the index of badChans (1..%d), or [] to skip: '...
    %                                         ,size(data,1)))];
   
    chanIDX(badChan) = []; 
    trainDat = data; 
    trainDat(badChan, :) = []; 
    data(badChan, :) = []; 
    if blinkChan == C
        blinkChan = size(data, 1); 
    end
    
    
    eyeBlinkDat = data(blinkChan,:) - mean(trainDat,1); 
    [b,a] = butter(4, [3,10]/(fs/2), 'bandpass');
    eyeBlinkDat = filtfilt(b,a, eyeBlinkDat); 
    eyeBlinkDat = (eyeBlinkDat - mean(eyeBlinkDat,2)) ./ ...
                        std(eyeBlinkDat,[],2);
        
    % %check that eyeblink dat looks as expected
    % figure; 
    % plot(eyeBlinkDat')

    test = eyeBlinkDat>2;
    blinkIndicator = reshape(test, 1, T, N); 


    %% ---------- Select an ICA training window ----------
    % Prefer a window with a moderate blink rate (enough blinks to define the
    % component, but not so many the data are dominated by artifact) and the
    % fewest interpolated samples. Falls back gracefully for short recordings
    % or when no window lands in the target blink-density band.
    winLen = 100000; 
    if length(test) <= winLen + 1
        % recording too short for windowed selection: train on all of it
        selStart = 1; 
        selEnd   = size(trainDat, 2); 
    else
        cand         = 1:500:(length(test) - winLen); 
        blinkCounts  = arrayfun(@(x) sum(test(x:x+winLen)), cand); 
        interpCounts = arrayfun(@(x) mean(interpChan(x:x+winLen)), cand); 
        sel = find(blinkCounts > prctile(blinkCounts, 75) & ...
                   blinkCounts < prctile(blinkCounts, 90)); 
        if isempty(sel)
            % nothing in the target band: fall back to the busiest window
            [~, sel] = max(blinkCounts); 
        end
        interpSel = interpCounts(sel); 
        sel = sel(min(interpSel) == interpSel);   % among those, fewest interpolated
        if numel(sel) > 1
            sel = sel(round(numel(sel)/2));        % middle of the ties
        end
        selStart = cand(sel);                      % FIX: sample index, not list position
        selEnd   = selStart + winLen; 
    end

    %% ---------- Learn the blink topography on the >2 Hz band ----------
    % ICA is trained on a 2 Hz high-passed copy so slow drift does not smear
    % or split the blink component. The resulting spatial filter is applied
    % to the broadband data below (blink topography is frequency-stable, so a
    % filter learned >2 Hz is valid across the full band). Fs is deliberately
    % NOT passed to ica_blinks, so it applies no additional high-pass/notch.
    [bHP, aHP] = butter(4, 2/(fs/2), 'high'); 
    trainDatHP = filtfilt(bHP, aHP, trainDat.').'; 
    trainDatHP = trainDatHP(:, selStart:selEnd); 

    out = ica_blinks(trainDatHP, 'blinkChan', blinkChan);

    if isfield(out, 'ambiguous') && out.ambiguous
        % Batch run could not auto-select a blink IC. Save a candidate-vs-blink
        % plot for manual review, then abort THIS session (the caller's
        % try/catch logs it as skipped). To reprocess after inspecting, set
        % ZLP_BLINK_IC to the chosen IC index and re-run just that session.
        try
            fA = figure('visible','off','Color','w'); nC = size(out.candSact,1);
            for ic = 1:nC
                subplot(nC,1,ic);
                zic = (out.candSact(ic,:)-mean(out.candSact(ic,:)))/max(std(out.candSact(ic,:)),eps);
                plot(zic); hold on; plot(out.blinkSig);
                ylabel(sprintf('IC %d', out.candIdx(ic)));
                if ic==1, title(sprintf('%s: choose blink IC, set ZLP\\_BLINK\\_IC', outDat.sessID)); end
            end
            if ~exist(outDat.figs,'dir'), mkdir(outDat.figs); end
            saveas(fA, fullfile(outDat.figs, ['blinkAmbiguous_' outDat.sessID '.jpg']));
            close(fA);
        catch, end
        error('ZLP:blinkAmbiguous', 'blink IC ambiguous for %s (candidate ICs %s)', ...
              outDat.sessID, mat2str(out.candIdx));
    end

    if ~isempty(out.badICs)
        ax = figure; 
        miniTopo(out.A(:,out.badICs(1)), outDat.eegLocs.X_flat(chanIDX), outDat.eegLocs.Y_flat(chanIDX)); 
      
        saveas(ax,fullfile(outDat.figs, 'removedBlink.jpg'));
    end

    %% ---------- Apply the spatial filter to the BROADBAND data ----------
    % Mean-centre first so that zeroing the blink IC does not subtract a
    % blink-shaped fraction of the channel means; add the means back after.
    dataMean = mean(data, 2); 
    Sclean   = out.W * (data - dataMean); 
    Sclean(out.badICs,:) = 0;                 % removal of blink IC entirely
    data_clean = out.A * Sclean + dataMean;   % back to channel space, means restored

    newEphys = outDat.data(1:32,:); 
    newEphys(chanIDX,:) = data_clean; 
    if ~origIs2D
        out = reshape(newEphys, C, T, N);
    else
        out = newEphys; 
    end
end
