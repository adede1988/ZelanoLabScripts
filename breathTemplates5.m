function [outMap]  = breathTemplates5(rspSig, fs)


[peaks, troughs] = findRespiratoryExtrema(rspSig, fs);

%I want to search from a trough to a peak to find an onset, but peaks are
%always before troughs in the current output. To adjust, throw away the
%first peak and the last trough
troughs = [max(1, peaks(1)-fs*2), troughs]; 

troughs(end) = []; 
smthRsp = smoothdata(rspSig, 'gaussian', round(fs/10));

[onsets2, diag] = findInhalationOnsets(smthRsp, fs, peaks, troughs); 


%col 1: onset Y value
%col 2: onset tim
%col 3: peak Y value
%col 4: peak tim
%col 5: end Y value
%col 6: end tim
%col 7: length (end tim - onset tim)
%col 8: amp (peak Y - avg of two ends)
%col 9: idx of peak in rspSig2
%col10: exhale peak Y value
%col11: exhale peak tim
%col12: condition
%col13: empty
%col14: index
tim = 1/fs:1/fs:length(rspSig)/fs;
outMap = zeros(length(peaks), 14); 
outMap(:,1) = rspSig(onsets2); 
outMap(:,2) = tim(onsets2); 
outMap(:,3) = rspSig(peaks); 
outMap(:,4) = tim(peaks); 
outMap(1:end-1,5) = rspSig(onsets2(2:end)); 
outMap(end,5) = mean(outMap(1:end-1,5)); 
outMap(1:end-1,6) = outMap(2:end,2);
outMap(end, 6) = min(tim(onsets2(end)) + 2, max(tim));
outMap(:,7) = outMap(:,6) - outMap(:,2); 
outMap(:,8) = outMap(:,3) - (outMap(:,1) + outMap(:,5)) / 2; 
outMap(:,9) = peaks; 
outMap(:,10) =  rspSig(troughs); 
outMap(:,11) = tim(troughs);
outMap(:,14) = 1:length(troughs); 







end
