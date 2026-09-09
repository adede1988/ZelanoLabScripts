function validate_slt()
% Validate slt_power_cont against Superlets/matlab-pure/faslt.m and sanity-check ridge_track.
here = fileparts(mfilename('fullpath')); addpath(here);
sup = 'C:\Users\Adam\Documents\GitHub\Superlets\matlab-pure';
if exist(sup,'dir'), addpath(sup); else
    sup2 = 'E:\GitHub\Superlets\matlab-pure'; if exist(sup2,'dir'), addpath(sup2); end
end

fs=500; T=4; t=(0:1/fs:T-1/fs);
rng(0);
x = 0.5*sin(2*pi*40*t) .* (t>1 & t<2) + 0.3*sin(2*pi*30*t) + 0.2*randn(1,numel(t));
F = 22:2:62; c1=3; ord=[3 30]; mult=1;

Pme = slt_power_cont(x, fs, F, c1, ord, mult);          % nF x N
if exist('faslt','file')==2
    Pref = faslt(x, fs, F, c1, ord, mult);              % nF x N
    % compare interior (avoid edges within longest wavelet)
    edge = 400;
    a = Pme(:, edge:end-edge); b = Pref(:, edge:end-edge);
    relerr = abs(a-b) ./ (abs(b)+eps);
    fprintf('SLT vs faslt: max rel err (interior) = %.3e ; median = %.3e ; corr = %.6f\n', ...
        max(relerr(:)), median(relerr(:)), corr(a(:), b(:)));
else
    fprintf('faslt.m not found on path -- skipping parity check (Superlets repo missing here).\n');
end

% ridge sanity: track on z-scored power of the 40 Hz burst region
addpath(fileparts(here));  % repo root not needed; myChanZscore is in repo
try
    Z = Pme;
    for f=1:size(Pme,1)
        col = Pme(f,:).';
        Z(f,:) = ( (col-mean(col))/std(col) ).';   % simple within-freq z for the test
    end
    [fr, ~, rp] = ridge_track(max(Z,0), F(:), 1.0, 1, 2);
    burstIdx = t>1.2 & t<1.8;
    fprintf('ridge median freq in 40Hz burst = %.1f Hz (expect ~40); mean ridge z there = %.2f\n', ...
        median(fr(burstIdx)), mean(rp(burstIdx)));
catch e
    fprintf('ridge sanity failed: %s\n', e.message);
end
fprintf('DONE validate_slt\n');
end
