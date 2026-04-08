function outDat = preprocess_eeg(outDat, standardEEGlocs, P)



    %check electrode names
    for chan = 1:32
        if ~strcmp(standardEEGlocs.Label(chan), outDat.labels(chan))
            error('EEG channels not labeled as expected')
        end
    end

    phi = standardEEGlocs.phi_deg .*pi ./ 180; 
    theta = standardEEGlocs.theta_deg .*pi ./ 180; 
    X = standardEEGlocs.X_nose_m;
    Y = standardEEGlocs.Y_left_m; 
    Z = standardEEGlocs.Z_up_m; 

    lbls = string(outDat.labels(1:32));
    lbls = reshape(lbls, [length(lbls), 1]);

    outDat.eegLocs = table(lbls, 'VariableNames', {'labels'}); 
    outDat.eegLocs.X = X; 
    outDat.eegLocs.Y = Y; 
    outDat.eegLocs.Z = Z; 
    outDat.eegLocs.theta = theta; 
    outDat.eegLocs.phi = phi; 
    outDat.eegLocs.X_flat = standardEEGlocs.X2D_right;
    outDat.eegLocs.Y_flat = standardEEGlocs.Y2D_front; 
    
    [badTS, badChans, outDat.data(1:32,:), interpChan] = ...
            removeNoiseChansVolt(outDat.data(1:32,:), outDat.fs, [1 32], outDat.eegLocs);
    title(outDat.sessID, 'interpreter', 'none')
    EXEMPT = [1 32];
     % only non-exempt get flagged
    badChans = setdiff(unique(badChans(:)), EXEMPT);         

    % 3) Build good-channel index for blink removal
    chanIDX = setdiff(1:32, badChans);% includes 1 & 32 even if detected bad

    % 4) Blink removal on good channels 
   
    [out, badChan2, blinkIndicator] = blinkRemoveWrapper(outDat,...
                                                chanIDX,...
                                                outDat.fs, interpChan);

    tmp = chanIDX(badChan2); 
    badChans = [badChans(:); tmp(:)]; 
  
    
    if ~isempty(badChans)
        ephysDat = interpolate_perrinX(out,X,Y,Z,badChans);
    else
        ephysDat = out;
    end
    ephysDat = ephysDat - mean(ephysDat,1); 

    dataLap = laplacian_perrinX(ephysDat, X, Y, Z); 

   
    outDat.badChans = outDat.labels(badChans); 
    outDat.dataLap = dataLap; 
    outDat.data(1:32,:) = ephysDat; 
    outDat.data(end+1, :,:) = blinkIndicator; 
    outDat.labels{end+1} = "blinkIndicator";
    outDat.data(end+1, :,:) = badTS; 
    outDat.labels{end+1} = "badTS";
    outDat.data(end+1, :,:) = interpChan; 
    outDat.labels{end+1} = "interpChan";
    outDat.EEGInterpolation = 1;
    outDat.EEGCleaning = 1;
    outDat.blinkRemoval = 1; 








end