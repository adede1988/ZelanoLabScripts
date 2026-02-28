function outDat = preprocess_macros(outDat, P)

    idx = cellfun(@(x) contains(x, 'macro'), outDat.labels);
    figure('visible', false, 'position', [0,0,1000,500])
    macroDat = outDat.data(idx, :); 
    plot(macroDat(1,:))
    hold on 
    for ii = 2:6
        plot(macroDat(ii,:)+(ii-1)*50)
    end
    legend()
    if ~isempty(P.macroRemove)
        plot(macroDat(P.macroRemove(1), :) + (P.macroRemove(1)-1)*50, ...
            'color', 'red')
    end
    title([outDat.sessID ' macros raw'], 'Interpreter','none')
    saveas(gcf,fullfile(outDat.figs, 'macrosRaw.jpg'));


    figure('visible', false, 'position', [0,0,1000,500])
    macroDat = outDat.data(idx, :); 
    plot(macroDat(1,50000:55000))
    hold on 
    for ii = 2:6
        plot(macroDat(ii,50000:55000)+(ii-1)*50)
    end
    legend()
    if ~isempty(P.macroRemove)
        plot(macroDat(P.macroRemove(1), 50000:55000) + ...
            (P.macroRemove(1)-1)*50, ...
            'color', 'red')
    end
    title([outDat.sessID ' macros raw 10s'], 'Interpreter','none')
    saveas(gcf,fullfile(outDat.figs, 'macrosRaw_zoom.jpg'));

    idx = cellfun(@(x) contains(x, 'macro'), outDat.labels);
    macroDat = outDat.data(idx, :); 
    if ~isempty(P.macroRemove)
        macroDat(P.macroRemove,:) = []; 
    end
    macOut = zeros([size(macroDat,1)-1, size(macroDat, [2,3])]);
    %do bipolar rereferencing 
    for chani = 1:size(macroDat,1)-1
        macOut(chani, :, :) = squeeze(macroDat(chani, :) -...
                                      macroDat(chani+1,:)); 
    end

if P.spikeClean
    % [b,a] = butter(4, [5,150]/(outDat.fs/2), 'bandpass');
    % gammaSig = filtfilt(b,a, macOut')'; 

    macOut = double(macOut); 

    splitFreq = 10;          % Hz cutoff between "low" and "high" components
    hpOrder   = 4;           % 4th-order Butterworth for high-pass
    
    % Design high-pass for the spike-y part (> splitFreq)
    [b_hp, a_hp] = butter(hpOrder, splitFreq/(outDat.fs/2), 'high');
    
    % High-frequency component of the IC
    x_high = filtfilt(b_hp, a_hp, macOut.').';   % column
    % Low-frequency residual (everything not captured by high-pass)
    x_low  = macOut - x_high;             

    z_high = (x_high - mean(x_high,2)) ./ std(x_high, [], 2); 
            
    [test, prominence] = detect_spikes(z_high, 2,...
        P.spikeWin); 
        %ICA is on the macro data after bipolar rereferencing because ICs
        %are more stable this way.
    out = ica_flag_spikes_targeted(x_high, test, prominence, ...
        'Fs', outDat.fs);
    %add low frequency component back into the clean data: 
    out.data_clean = out.data_clean + x_low; 
    whereSpikes = movmean(out.mixVector, 10*outDat.fs);
    [~, idx] = min(whereSpikes(20000:end-20000)); 
    idx = idx + 10000; 
    idx = min(size(macOut,2)-outDat.fs*5-1, idx);
    x = figure('visible', false, 'position', [0,0,1000,500]);
    plot(macOut(1,idx-outDat.fs*5:idx+outDat.fs*5))
    hold on 
    plot(out.data_clean(1,idx-outDat.fs*5:idx+outDat.fs*5))
    for ii = 2:size(macOut,1)
    plot(macOut(ii,idx-outDat.fs*5:idx+outDat.fs*5) + ii*30)
    hold on 
    plot(out.data_clean(ii,idx-outDat.fs*5:idx+outDat.fs*5)+ ii*30)


    end
    
    title([outDat.sessID ' spike removal'], 'interpreter', 'none')
    saveas(x, fullfile(outDat.figs, 'macroSpikeRemoval.jpg'));

    C = size(macOut,1);
    newLabs = {'macBP1', 'macBP2', 'macBP3','macBP4', 'macBP5'};
    newLabs = newLabs(1:C);  
    
    outDat.data(end+1:end+C, :) = out.data_clean; 
    outDat.labels(end+1:end+C) = newLabs; 
    outDat.data(end+1, :) = out.mixVector; 
    outDat.labels{end+1} = "spikeCleanVec";
    outDat.spikeRemoval = 1; 
    

    

else
    C = size(macOut,1);
    newLabs = {'macBP1', 'macBP2', 'macBP3','macBP4', 'macBP5'};
    newLabs = newLabs(1:C); 
   


    outDat.data(end+1:end+C, :) = macOut;
    outDat.labels(end+1:end+C) = newLabs;

    outDat.data(end+1, :) = ones(size(outDat.data,2),1); 
    outDat.labels{end+1} = "spikeCleanVec";
    outDat.spikeRemoval = 1; 


end





end