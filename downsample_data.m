function outDat = downsample_data(outDat, fs_new)
% outDat: struct containing .data (ch x time OR ch x time x trials), .fs
% fs_new: desired sample rate (e.g., 500)

    fs_old = outDat.fs;

    % Compute resampling factors
    [p, q] = rat(fs_new/fs_old);  
    
    % Get data and size
    data = outDat.data;
    dims = ndims(data);
    sz = size(data);
    
    % Handle 2D or 3D input
    if dims == 2
        % channels x time
        ch = sz(1);
        t = sz(2);
        tr = 1;
        data = reshape(data, ch, t, tr);
    elseif dims == 3
        % channels x time x trials
        [ch, t, tr] = size(data);
    else
        error('Input data must be 2D or 3D: channels x time OR channels x time x trials.');
    end
    
    % Preallocate
    new_t = ceil(t * p / q);
    data_ds = zeros(ch, new_t, tr);
    
    % Design filters
    highd = designfilt('highpassiir','FilterOrder',4, ...
        'HalfPowerFrequency',0.03,'SampleRate',fs_new);
    lowd = designfilt('lowpassiir','FilterOrder',4, ...
        'HalfPowerFrequency',fs_new/2 - 0.1,'SampleRate',fs_new); % Nyquist-safe
    
    % Process each channel & trial
    for ii = 1:ch
        

        for jj = 1:tr
            sig = double(squeeze(data(ii,:,jj)));
            sig = resample(sig, p, q);
            sig = filtfilt(highd, sig);
            sig = filtfilt(lowd, sig);

            % NOTE: line-noise notch filtering (60/120/180 Hz) intentionally
            % removed -- downstream analyses handle line noise themselves.

            data_ds(ii,:,jj) = sig;
        end
    end
    
   


    % Restore original shape
    if dims == 2
        data_ds = squeeze(data_ds); % back to ch x time
    end
    
    % Update struct
    outDat.data = data_ds;
    outDat.fs = fs_new;
    outDat.origFS = fs_old;
    outDat.downsampled = 1;
    
    % Optional: update time vector if available
    if isfield(outDat,'tim')
        outDat.tim = (1:size(data_ds,2)) / fs_new;
    end
end
