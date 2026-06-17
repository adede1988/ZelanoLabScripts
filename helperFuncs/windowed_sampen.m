function [sampen_vals, t_idx] = windowed_sampen(x, win_size, step, m, r)
% x         : input signal (vector)
% win_size  : window size (samples)
% step      : step size (samples)
% m         : embedding dimension (e.g., 2)
% r         : tolerance (absolute, not scaled internally)

x = x(:); % ensure column vector
N = length(x);

n_windows = floor((N - win_size) / step) + 1;
sampen_vals = nan(n_windows,1);
t_idx = nan(n_windows,1);

parfor i = 1:n_windows
    disp([num2str(i) ' of ' num2str(n_windows)])
    start_idx = (i-1)*step + 1;
    end_idx = start_idx + win_size - 1;
    
    segment = x(start_idx:end_idx);
    
    % optional detrend (recommended for respiration)
    segment = detrend(segment);
    curR = r; %0.2 * std(segment); 
    sampen_vals(i) = sampen(segment, m, curR);
    
    % center index for time reference
    t_idx(i) = start_idx + floor(win_size/2);
end
end

function se = sampen(x, m, r)
% Sample entropy (SampEn)
% x : signal (vector)
% m : embedding dimension
% r : tolerance

x = x(:);
N = length(x);

% build template vectors
count_m = 0;
count_m1 = 0;

for i = 1:N - m
    template_m = x(i:i+m-1);
    template_m1 = x(i:i+m);
    
    for j = i+1:N - m
        comp_m = x(j:j+m-1);
        
        if max(abs(template_m - comp_m)) <= r
            count_m = count_m + 1;
            
            % check m+1 match
            if j <= N - m - 1
                comp_m1 = x(j:j+m);
                if max(abs(template_m1 - comp_m1)) <= r
                    count_m1 = count_m1 + 1;
                end
            end
        end
    end
end

% avoid division by zero
if count_m == 0 || count_m1 == 0
    se = NaN;
else
    se = -log(count_m1 / count_m);
end
end