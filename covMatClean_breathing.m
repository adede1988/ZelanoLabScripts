function [badIDX] = covMatClean_breathing(data, starts, stops)

%data should be in Channel X time X trial configuration


[chanCount, tLength, Tcount] = size(data);  
%preallocate for the covariance matrices
covMats = nan(chanCount, chanCount, Tcount);  

%get cov matrices in every window
for win = 1:Tcount
    if stops(win)<10000
   covMats(:,:,win) = (squeeze(data(:,starts(win):stops(win),win)) * ...
                       squeeze(data(:,starts(win):stops(win),win))') ./ ...
                       (stops(win) - starts(win));
    end
end

%get an average covMat
covMatAvg = mean(covMats,3, 'omitnan'); 

%find the sum of squared difference of each epoch from the average
diffs = zeros(size(data,3),1);
for win = 1:size(data,3)
   cur = covMats(:,:,win);  
   diffs(win) = sum((cur(:) - covMatAvg(:)).^2, 'omitnan');   
end

% z-score the SS differences 
cutDif = diffs(diffs<prctile(diffs, 95)); 
diffsZ = (diffs - mean(cutDif)) / std(cutDif); 

%1=window that starts at the corresponding windowStarts index is good
%0=window that starts at the corresponding windowStarts index is bad
goodWindows = ones(Tcount, 1);
%z-score threshold 
goodWindows(diffsZ>2) = 0; 
goodWindows(diffsZ<-2) = 0; 
badIDX = find(goodWindows == 0); 



end