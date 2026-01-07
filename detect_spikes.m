function [flags, prominence] = detect_spikes(data, threshold, winSize)
% DETECT_SPIKES detect sharp spike-like deflections in 3D data
% 
%   Inputs:
%     data            - C x T x N  (Channels x Time x Trials)
%     threshold       - scalar or Cx1 vector giving the threshold for (max-min)
%     winSize         - odd integer window length (default = 5)

%
%   Outputs:
%     flags    - binary matrix C x T x N, 1 = spike detected (at that time sample)
%     prominence - per-spike local SNR-like prominence; zeros where no spike
%
% Example:
%   [f,r] = detect_spikes(myData, 50, 5);
%

% Defaults

% Accept 2D (C x T) or 3D (C x T x N)
origIs2D = ismatrix(data);
if origIs2D
    data    = reshape(data,    size(data,1), size(data,2), 1);
end

[C, T, N] = size(data);


% ---- compute moving max - min with centered window ----
half = floor((winSize-1)/2);
left = half;
right = winSize - 1 - left; % equals half for odd winSize

% Use movmax - movmin with centered window [left right] along time (dim=2)
rangeVals = movmax(data, [left right], 2) - movmin(data, [left right], 2);

% ---- detection ----
detections = rangeVals > threshold;
prominence = zeros(size(detections));  
outDetections = zeros(size(detections)); 

chGlobalMAD = zeros(size(data,1),1);
for ch = 1:size(data,1)
    x = reshape(data(ch,:,:), 1, []);
    chGlobalMAD(ch) = mad(x, 0);  % MATLAB's MAD; already ~1.4826*median|...|
end
flankHalf = 50;       % your ~100-sample total flank width
madFloorFrac = 0.3;   % 30% of global MAD as a safety floor




for ch = 1:C
   
    limVals = prctile(real(data(ch,:,:)), [10, 90], [2, 3]);
    lower = limVals(1); upper = limVals(2);
    
    for tr = 1:N

        sig = squeeze(data(ch,:,tr));  % 1 x T
        curIdx = find(detections(ch,:,tr)); 
        

        %detections of spikes: 
        posIdx = curIdx(sig(curIdx) > upper);
        negIdx = curIdx(sig(curIdx) < lower); 
        
        %detections of below to above crossings for upper threshold: 
        midToHigh = find(sig(1:end-1) < upper & sig(2:end)>upper); 
        highToMid = find(sig(1:end-1) > upper & sig(2:end)<upper); 

        %% ----- map of nearest midToHigh to the LEFT of each sample -----
        leftMap = zeros(1, T);              % 0 means "no crossing yet"
        leftMap(midToHigh) = midToHigh;     % put crossing indices at their locations
        leftMap(1) = 1; 
   
        leftMap = cummax(leftMap); %get prior max to each point

        % For each spike, left shoulder (mid→high crossing) is:
        leftShoulder = leftMap(posIdx);     % same length as posIdx


        %% ----- map of nearest highToMid to the RIGHT of each sample -----
        rightMap = zeros(1, T);
        rightMap(highToMid) = highToMid;    % mark crossings
        
        % Treat zeros as end of vector so cummin ignores them
        rightMap(rightMap == 0) = T;
        
        % cumulative min from the right: at each t, store the first crossing ≥ t
        rightMap = fliplr( cummin(fliplr(rightMap)) );   % 1 x T
        
        % For each spike, right shoulder (high→mid crossing) is:
        rightShoulder = rightMap(posIdx);    

        % eliminate shoulders too narrow: 
        badIdx = leftShoulder == posIdx | rightShoulder == posIdx; 
    
        leftShoulder(badIdx) = []; 
        rightShoulder(badIdx) = []; 
        posIdx(badIdx) = []; 

        prominence(ch,posIdx,tr) = arrayfun(@(x,y) local_prominence_snr( ...
                                    sig, x, y, flankHalf, 'pos', ...
                                    chGlobalMAD(ch), madFloorFrac), ...
                                    leftShoulder, rightShoulder);

        outDetections(ch,posIdx,tr) = true; 
        %now repeat for negative spikes: 


        %detections of below to above crossings for upper threshold: 
        midToLow = find(sig(1:end-1) > lower & sig(2:end)<lower); 
        lowToMid = find(sig(1:end-1) < lower & sig(2:end)>lower); 

        %% ----- map of nearest midToLow to the LEFT of each sample -----
        leftMap = zeros(1, T);              % 0 means "no crossing yet"
        leftMap(midToLow) = midToLow;     % put crossing indices at their locations
        leftMap(1) = 1; 
   
        leftMap = cummax(leftMap); %get prior max to each point

        % For each spike, left shoulder (mid→high crossing) is:
        leftShoulder = leftMap(negIdx);     % same length as negIdx


        %% ----- map of nearest lowToMid to the RIGHT of each sample -----
        rightMap = zeros(1, T);
        rightMap(lowToMid) = lowToMid;    % mark crossings
        
        % Treat zeros as end of vector so cummin ignores them
        rightMap(rightMap == 0) = T;
        
        % cumulative min from the right: at each t, store the first crossing ≥ t
        rightMap = fliplr( cummin(fliplr(rightMap)) );   % 1 x T
        
        % For each spike, right shoulder (high→mid crossing) is:
        rightShoulder = rightMap(negIdx);    

        % eliminate shoulders too narrow: 
        badIdx = leftShoulder == negIdx | rightShoulder == negIdx; 
    
        leftShoulder(badIdx) = []; 
        rightShoulder(badIdx) = []; 
        negIdx(badIdx) = []; 

        prominence(ch,negIdx,tr) = arrayfun(@(x,y) local_prominence_snr( ...
                                    sig, x, y, flankHalf, 'neg', ...
                                    chGlobalMAD(ch), madFloorFrac), ...
                                    leftShoulder, rightShoulder);



        outDetections(ch,negIdx,tr) = true; 

    end
end







% If original input was 2D, squeeze outputs back to C x T
if origIs2D
    flags      = squeeze(outDetections);
    prominence = squeeze(prominence);
else
    flags = outDetections; 
end

end




% 
% 
%         for k = 1:T
%             if detections(ch,k,tr) && ...
%                (data(ch,k,tr) < lower || data(ch,k,tr) > upper)
%                 % --------- POSITIVE spike ---------
%                 if sig(k) > upper
%                     % first index to the RIGHT that falls back below 'upper'
%                     right = find(sig(k:end) < upper, 1, 'first');
%                     % make absolute
%                     if ~isempty(right), right = k + right - 1; end  
%                     % last index to the LEFT that is below 'upper'
%                     left  = find(sig(1:k) < upper, 1, 'last');
% 
% 
%                     if isempty(left)
%                         left = 1;
% 
%                     end
%                     if isempty(right)
%                         right = T; 
% 
%                     end
%                     if ~isempty(right) && ~isempty(left)
%                         prominence(ch,k,tr) = local_prominence_snr( ...
%                                     sig, left, right, flankHalf, 'pos', ...
%                                     chGlobalMAD(ch), madFloorFrac);
%                         % prominence(ch,k,tr) = cenMax / max(Rmax, Lmax, eps,...
%                         %                         'omitnan');
%                     else
%                         prominence(ch,k,tr) = NaN; % no context available
%                     end
%                 % --------- NEGATIVE spike ---------
%                 elseif sig(k) < lower
%                     % first index to the RIGHT that rises back above 'lower'
%                     right = find(sig(k:end) > lower, 1, 'first');
%                     if ~isempty(right), right = k + right - 1; end
% 
%                     % last index to the LEFT that is above 'lower'
%                     left  = find(sig(1:k) > lower, 1, 'last');
% 
%                     if isempty(right)
%                         right = T; 
% 
%                     end
%                     if isempty(left)
%                         left = 1; 
% 
%                     end
%                     if ~isempty(right) && ~isempty(left)
%                         prominence(ch,k,tr) = local_prominence_snr( ...
%                                     sig, left, right, flankHalf, 'neg', ...
%                                     chGlobalMAD(ch), madFloorFrac);
%                         % prominence(ch,k,tr) = cenMax / max(Rmax, Lmax, eps,...
%                         %                         'omitnan');
%                     else
%                         prominence(ch,k,tr) = NaN; % no context available
%                     end
%                 end
%             end
%         end
%     end
% end
