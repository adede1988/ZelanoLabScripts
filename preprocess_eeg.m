function outDat = preprocess_eeg(outDat, standardEEGlocs, P)



    %check electrode names
    for chan = 1:32
        if ~strcmp(standardEEGlocs.Label(chan), outDat.labels(chan))
            error('EEG channels not labeled as expected')
        end
    end

    phi = standardEEGlocs.Phi .*pi ./ 180; 
    theta = standardEEGlocs.Theta .*pi ./ 180; 
    X = cos(phi) .* sin(theta);
    Y = sin(phi) .* sin(theta); 
    Z = cos(theta); 

    lbls = string(outDat.labels(1:32));
    lbls = reshape(lbls, [length(lbls), 1]);

    outDat.eegLocs = table(lbls, 'VariableNames', {'labels'}); 
    outDat.eegLocs.X = X; 
    outDat.eegLocs.Y = Y; 
    outDat.eegLocs.Z = Z; 
    outDat.eegLocs.theta = theta; 
    outDat.eegLocs.phi = phi; 
    
    [badTS, badChans] = ...
            removeNoiseChansVolt(outDat.data(1:32,:), outDat.fs);
    chanIDX = 1:32; 
    if ismember(1, badChans)
      chanIDX(badChans(2:end)) = [];
    else
      chanIDX(badChans) = [];
      
    end
    ephysDat = outDat.data(1:32,:);
    [out, badChan2, blinkIndicator] = blinkRemoveWrapper(outDat,...
                                                chanIDX,...
                                                outDat.fs);

    
    badChans = [badChans; chanIDX(badChan2)]; 
    ephysDat(chanIDX,:) = out; 
    

    ephysDat = interpolate_perrinX(ephysDat,X,Y,Z,badChans);

    ephysDat = ephysDat - mean(ephysDat,1); 

    dataLap = laplacian_perrinX(ephysDat, X, Y, Z); 

   
    outDat.badChans = outDat.labels(badChans); 
    outDat.dataLap = dataLap; 
    outDat.data(1:32,:) = ephysDat; 
    outDat.data(end+1, :,:) = blinkIndicator; 
    outDat.labels{end+1} = "blinkIndicator";
    outDat.data(end+1, :,:) = badTS; 
    outDat.labels{end+1} = "badTS";
    outDat.EEGInterpolation = 1;
    outDat.EEGCleaning = 1;
    outDat.blinkRemoval = 1; 








end