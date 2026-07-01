function y = chirp_bbfilt(x, fs, band, C)
% CHIRP_BBFILT  Zero-phase FIR bandpass (firws if FieldTrip present, else fir1+filtfilt).
%   y = chirp_bbfilt(x, fs, band, C)   % x [nRow x nSamp]; filtered row-wise; zero phase.
%   Used for ALL hilbert-based phase/envelope substrates (F4): broadband 25-58 and narrow
%   band-limited [f+-band]. Operate on the PADDED epoch then trim downstream so FIR edge
%   transients live in the pad.

    if nargin < 4, C = struct(); end
    if isrow(x) || iscolumn(x), x = x(:)'; end
    use_ft = ~isempty(which('ft_preproc_bandpassfilter'));
    if isfield(C,'bbFiltType'), ftype = C.bbFiltType; else, ftype = 'firws'; end

    if use_ft
        try
            y = ft_preproc_bandpassfilter(x, fs, band, [], ftype);
            return;
        catch
            % fall through to fir1
        end
    end
    % fir1 + filtfilt fallback (zero-phase). Order ~ 3 cycles of the low edge.
    ord = 2*round(1.5*fs/band(1));            % even order
    ord = min(ord, floor((size(x,2)-1)/3));   % filtfilt needs 3*order < nSamp
    b = fir1(ord, band/(fs/2), 'bandpass');
    y = filtfilt(b, 1, x')';
end
