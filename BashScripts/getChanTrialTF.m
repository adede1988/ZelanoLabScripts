function [trialTF] = getChanTrialTF(trialDat, frex, numfrex, stds, srate)

%trialDat:            timepoints X trials data points


%trialTF is going to be the hilbert transformed complex time series of all
%trials across a range of frequencies
%time X trials X freq  
trialTF = zeros([size(trialDat), numfrex]); 

padDat = mirrorPad(trialDat); 

time  = -1:1/srate:1; % time, from -1 to 1 second in steps of 1/sampling-rate

n_wavelet            = length(time);
n_data               = prod(size(padDat)); 
n_convolution        = n_wavelet+n_data-1;
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;

% hz = linspace(0, srate, size(n_convolution,1));
    
  
fftDat = fft(reshape(padDat,1,numel(padDat)),n_conv_pow2);
    
    for fi = 1:numfrex

        f  = frex(fi); % frequency of wavelet in Hz
        s  = stds(fi); 
        % and together they make a wavelet
        wavelet = sqrt(1/(s*sqrt(pi))) * ...
            exp(2*1i*pi*f.*time) .* exp(-time.^2./(2*(s^2))); 
        
        wavelet = fft(wavelet, n_conv_pow2); 

        eegconv = ifft(wavelet.*fftDat);
        eegconv = eegconv(1:n_convolution);
        eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
        temppower = reshape(eegconv,size(padDat));
        outTF = temppower(size(trialDat,1)+1:size(trialDat,1)*2, :); 


       
        trialTF(:,:,fi) = outTF; 
       
    end
        



end