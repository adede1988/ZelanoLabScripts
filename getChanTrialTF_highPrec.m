function [trialTF] = getChanTrialTF_highPrec(trialDat, frex, numfrex, stds, srate)

%trialDat:            timepoints X trials data points


%trialTF is going to be the hilbert transformed complex time series of all
%trials across a range of frequencies
%time X trials X freq  
trialTF = zeros([size(trialDat), numfrex]); 

padDat = mirrorPad(trialDat); 

time  = -10:1/(srate/10):10; % time, from -1 to 1 second in steps of 1/sampling-rate

n_wavelet            = length(time);
n_data               = prod(size(padDat)); 
n_convolution        = n_wavelet+n_data-1;
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;

% hz = linspace(0, srate, size(n_convolution,1));
    
  
fftDat = fft(reshape(padDat,1,numel(padDat)),n_conv_pow2);
    
    parfor fi = 1:numfrex
        fi
        tic
        f  = frex(fi)/10; % frequency of wavelet in Hz
        s  = stds(fi)*10; % / (2*pi*f); 
        

        %gaussian: 
        GaussWin = exp(-(time.^2)/(2*s^2));
        SinWave = exp(1i*2*pi*f.*time);
        %scaling factor (only important if you're not going to normalize or
        %z score: 
        A = 1 / sqrt(s*sqrt(pi)); 
        % and together they make a wavelet
        wavelet = GaussWin .* SinWave * A; 
        
        wavelet = fft(wavelet, n_conv_pow2); 

        eegconv = ifft(wavelet.*fftDat);
        eegconv = eegconv(1:n_convolution);
        eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
        temppower = reshape(eegconv,size(padDat));
        outTF = temppower(size(trialDat,1)+1:size(trialDat,1)*2, :); 

%         % create Gaussian
%         s  = stds(fi)*(2*pi-1)/(4*pi); % normalized width
%         x  = hz-frex(fi);                 % shifted frequencies (pre insert the mean for the gaussian)
%         fx = exp(-.5*(x/s).^2);    % gaussian
%         fx = fx./abs(max(fx));     % gain-normalized

       
        trialTF(:,:,fi) = outTF; 
       toc
    end
        



end