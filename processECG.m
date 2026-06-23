function outDat = processECG(outDat, P)

    % ECG channels -> band-passed, z-scored, with the min beat separation.
    % (Shared with paramCheckECG so the verification figure matches detection.)
    [ECGz, beatSep] = buildECGz(outDat);



    heartBeats = P.getBeats(ECGz, beatSep);



    figure('visible', false, 'position', [0,0,1000,500]);
    plot(ECGz(1,100000:110000), 'color', 'k')
    hold on 
    plot(ECGz(2,100000:110000), 'color', 'red')
    plot(ECGz(3,100000:110000), 'color', 'green')
    xlim([0 10000])
    xticks([0:1000:10000])
    xticklabels(0:2:20)
    xlabel('Time (s)')
    xline(heartBeats(heartBeats>100000 & heartBeats<110000)-100000, 'color', 'magenta', 'linestyle', '--')
    maxVal = max(ECGz(:,100000:110000), [], 'all');
    minVal = min(ECGz(:,100000:110000), [], 'all'); 
    ylim([minVal*1.5 maxVal * 1.5])
    title(sprintf('ECG beat detection (%s)', ...
                  outDat.sessID), ...
          'Interpreter','none');
    saveas(gcf,fullfile(outDat.figs, ['ECG_beatDetect' '.jpg']));


    outDat.heartBeats = heartBeats; 


    %% with beats detected, establish RR interval timeseries
    
    idx = cellfun(@(x) contains(x, 'rsp'), outDat.labels);
    rspDat = outDat.data(idx, :); 
    
    %choose which one looks right
    % rspDatz = (rspDat - mean(rspDat, [2,3])) ./ std(rspDat, [], [2,3]); 
    
    idx = P.rspIDX; 
    rspDat = squeeze(rspDat(idx, :));  
    
    %flip signal
    rspDat = rspDat .* P.rspFlip;

    %convert beats indecies to beat times
    tim = 1/outDat.fs:1/outDat.fs:size(outDat.data,2)/outDat.fs; 
    beatTims = tim(outDat.heartBeats); 
    beatDiffs = diff(beatTims);
    beatTims(end) = [];
       
   
    %account for accidental misses and double counts and interpolate
    %across
    breakVals = [.2:.05:10]; 
    counts = arrayfun(@(x,y) sum(beatDiffs>x & beatDiffs<y), ...
    breakVals(1:end-1), breakVals(2:end));
    %locate the mode and find zeros around it to define central dist.
    [~, idx] = max(counts); 
    minVal = breakVals(idx - find(flip(counts(1:idx))<5, 1) + 1); 
    if isempty(minVal)
        minVal = .6; %interbeat interval of less than .6 is unrealistic for generally healthy people sitting still
    end
    maxVal = breakVals(idx + find(counts(idx:end)<5, 1) - 1); 

    %remove bad vals: 
    beatTims(beatDiffs < minVal) = []; 
    beatDiffs(beatDiffs < minVal) = []; 

    beatTims(beatDiffs > maxVal) = []; 
    beatDiffs(beatDiffs > maxVal) = []; 
    figure('visible', false, 'position', [0,0,1000,500]);
    histogram(beatDiffs)
    xlabel('interbeat interval (s)')  
    title(sprintf('Interbeat heart interval (%s)', ...
                  outDat.sessID), ...
          'Interpreter','none');
    saveas(gcf,fullfile(outDat.figs, ['interbeatHist' '.jpg']));

    %interpolate the interbeat interval into a full timeseries
    RRint = interp1(beatTims,beatDiffs, tim, 'linear');
    
    
    figure('visible', false, 'position', [0,0,1000,500]);
    plot(rspDat(100000:150000))
    ylabel('Respiration velocity')
    yyaxis right
    plot(RRint(100000:150000))
    ylabel('Interbeat Interval (s)')
    xlim([0 50000])
    xticks([0:5000:50000])
    xticklabels([0:10:100])
    title(sprintf('Respiration and HRV (%s)', ...
                  outDat.sessID), ...
          'Interpreter','none');
    saveas(gcf,fullfile(outDat.figs, ['RespirationHeart' '.jpg']));

    outDat.data(end+1, :) = RRint; 
    outDat.labels{end+1} = 'RRint'; 




end