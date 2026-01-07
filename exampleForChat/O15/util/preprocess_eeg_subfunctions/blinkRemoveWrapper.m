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
    
    figure; 
    imagesc(data) %check for bad channels overall 
    caxis([-200, 200])
    
    badChan = input(sprintf(...
    'Enter the index of badChans (1..%d), or [] to skip: '...
                                            ,size(data,1)));
    
    figure; 
    imagesc(data) %check for bad channels overall 
    caxis([-1, 1])
    
    badChan = [badChan input(sprintf(...
    'Enter the index of badChans (1..%d), or [] to skip: '...
                                            ,size(data,1)))];
    chanidx = 1:size(data,1);
    trainDat = data; 
    trainDat(badChan, :) = []; 
    chanidx(badChan) = []; 
    
    
    eyeBlinkDat = data(blinkChan,:) - mean(trainDat,1); 
    [b,a] = butter(4, [3,10]/(fs/2), 'bandpass');
    eyeBlinkDat = filtfilt(b,a, eyeBlinkDat); 
    eyeBlinkDat = (eyeBlinkDat - mean(eyeBlinkDat,2)) ./ ...
                        std(eyeBlinkDat,[],2);
        
    %check that eyeblink dat looks as expected
    figure; 
    plot(eyeBlinkDat')

    test = eyeBlinkDat>2;
    blinkIndicator = reshape(test, 1, T, N); 
       

       
    startIdx = [1:500:length(test)-100000]; 
    blinkCounts = arrayfun(@(x) sum(test(x:x+100000)), startIdx); 
    idx = find(blinkCounts>median(blinkCounts)-5 & ...
        blinkCounts<median(blinkCounts)+5);


    startIdx = startIdx(idx(1)); 




    trainDat = trainDat(:,startIdx:startIdx+100000);

    out = ica_blinks(trainDat, 'blinkChan', ...
        find(chanidx == blinkChan));

    ax = plotICATopo(out, out.badICs(1), outDat.eegLocs.theta(chanidx), ...
        outDat.eegLocs.phi(chanidx));
    saveas(ax,fullfile(outDat.figs, 'removedBlink.jpg'));

    Sclean = out.W * data(chanidx,:); 
    Sclean(out.badICs,:) = 0; % removal of blink IC entirely 
    data_clean = out.A * Sclean;                     % back to channel space
    newEphys = data; 
    newEphys(chanidx,:) = data_clean; 
    if ~origIs2D
        out = reshape(newEphys, C, T, N);
    else
        out = newEphys; 
    end
end