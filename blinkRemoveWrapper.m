function [out, badChan, blinkIndicator] = blinkRemoveWrapper(outDat, ...
                                                            chanIDX, fs)

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
    
    highRangeThresh = 700;   % µV: max - min ≥ 700 counts as a "bad" epoch
    epochFracThresh = 0.15;  % 25% of epochs
    
    % Range per channel per epoch
    epochMax   = squeeze(max(dataEpochs, [], 2));   % [nCh x nEpochs]
    epochMin   = squeeze(min(dataEpochs, [], 2));   % [nCh x nEpochs]
    epochRange = epochMax - epochMin;               % [nCh x nEpochs]
    
    % Epoch is "bad" if range ≥ 200 µV
    epochBad = epochRange >= highRangeThresh;       % logical [nCh x nEpochs]
    
    % Channel is bad if ≥ 25% of its epochs are bad
    badHighDev = (sum(epochBad, 2) ./ nEpochs) >= epochFracThresh;   % [nCh x 1]
    
    
    %% ---------- Criterion 2: flat periods (> 100 ms within ±1 µV) ----------
    
    flatAmpThresh      = 1;                % µV: between -1 and +1
    flatLenSec         = 0.100;            % 100 ms
    flatLenSamples     = round(flatLenSec * fs);
    flatTotalSecThresh = 30;                % 5 seconds total
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
    
        % Mark channel as bad if total flat time ≥ 5 s
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
       

       
    startIdx = [1:500:length(test)-100000]; 
    blinkCounts = arrayfun(@(x) sum(test(x:x+100000)), startIdx); 
    idx = find(blinkCounts>median(blinkCounts)-100 & ...
        blinkCounts<median(blinkCounts)+100);


    startIdx = startIdx(idx(1)); 




    trainDat = trainDat(:,startIdx:startIdx+100000);

    out = ica_blinks(trainDat, 'blinkChan', ...
        find(chanIDX == blinkChan));
    if ~isempty(out.badICs)
        ax = plotICATopo(out, out.badICs(1), outDat.eegLocs.theta(chanIDX), ...
            outDat.eegLocs.phi(chanIDX));
        saveas(ax,fullfile(outDat.figs, 'removedBlink.jpg'));
    end
    Sclean = out.W * data; 
    Sclean(out.badICs,:) = 0; % removal of blink IC entirely 
    data_clean = out.A * Sclean;                     % back to channel space
    newEphys = outDat.data(1:32,:); 
    newEphys(chanIDX,:) = data_clean; 
    if ~origIs2D
        out = reshape(newEphys, C, T, N);
    else
        out = newEphys; 
    end
end