function zDat = myChanZscore(varargin) 
   
%data:              time X trials
%basePeriod:        [start, stop] indices from which to draw the baseline


%zDat:              trials X time 




switch nargin
    case 1
        data = varargin{1}; 
        basePeriod = [1, size(data,1)];
        
    case 2
        data = varargin{1}; 
        basePeriod = varargin{2}; 
    otherwise
        warning('Error: at least one input is needed')
        return
end






nsamps = 1000; 


ntrials = size(data,2); 


tmp = data(basePeriod(1):basePeriod(2), :);
tmp = tmp(:); 

dist = zeros(nsamps,1); 
for ii = 1:nsamps
    test = randsample(tmp, ntrials); 
    dist(ii) = mean(test,1, 'omitnan'); 
end

zDat = (data - mean(dist)) ./ std(dist);








end




%permute it to channels X time X trials
% data = permute(data, [2,3,1]); 