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
    title([outDat.sessID ' macros raw 10s'], 'Interpreter','none')
    saveas(gcf,fullfile(outDat.figs, 'macrosRaw_zoom.jpg'));

    idx = cellfun(@(x) contains(x, 'macro'), outDat.labels);
    macroDat = outDat.data(idx, :); 
    macOut = zeros([5, size(macroDat, [2,3])]);
    %do bipolar rereferencing 
    for chani = 1:5
        macOut(chani, :, :) = squeeze(macroDat(chani, :) -...
                                      macroDat(chani+1,:)); 
    end

if P.spikeClean
    [b,a] = butter(4, [5,150]/(outDat.fs/2), 'bandpass');
    gammaSig = filtfilt(b,a, macOut')'; 
    [test, prominence] = detect_spikes(macOut,P.spikeThresh,...
        P.spikeWin,...
        false, gammaSig); 
        %ICA is on the macro data without bipolar rereference 
        %this allows later rereferencing at will
    out = ica_flag_spikes_targeted(macOut, test, prominence, ...
        'Fs', outDat.fs);

    whereSpikes = movmean(out.mixVector, 10*outDat.fs);
    [~, idx] = min(whereSpikes); 
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

    outDat.data(end+1:end+5, :) = out.data_clean; 
    outDat.labels(end+1:end+5) = {'macBP1', 'macBP2', 'macBP3', ...
                                                'macBP4', 'macBP5'};
    outDat.data(end+1, :) = out.mixVector; 
    outDat.labels{end+1} = "spikeCleanVec";
    outDat.spikeRemoval = 1; 

    

else
    outDat.data(end+1:end+5, :) = macOut;
    outDat.labels(end+1:end+5) = {'macBP1', 'macBP2', 'macBP3', ...
                                                'macBP4', 'macBP5'};
    outDat.data(end+1, :) = ones(size(outDat.data,2),1); 
    outDat.labels{end+1} = "spikeCleanVec";
    outDat.spikeRemoval = 1; 


end





end