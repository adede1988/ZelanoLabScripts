% run_chirpV2_job_lab.m  -- WMI-detached chirpAnalysisV2 batch launcher for the lab E: session.
% Launch via:
%   ([wmiclass]'Win32_Process').Create('"matlab.exe" -batch "run(''E:/GitHub/ZelanoLabScripts/chirpAnalysis/run_chirpV2_job_lab.m'')"')
% Log E:\chirpV2_job.log; sentinel E:\_CHIRPV2_DONE.txt (or _CHIRPV2_ERROR.txt).
% Writes ridgeInfo + powerPhaseContinuity into each cue final on E: (overwrite) and figures under
% E:\Lab_Common\Adam\Dupi_processing\<id>\singleTrialSpectrograms. Copyback E:->R: is done from home.

diary('E:\chirpV2_job.log');
try
    fprintf('=== chirpV2 job start %s ===\n', datestr(now));
    addpath('E:\GitHub\ZelanoLabScripts');                 % labPaths.m
    addpath('E:\');                                        % labPaths_local (codePre=E:\GitHub, fieldtrip=E:)
    addpath('E:\GitHub\ZelanoLabScripts\chirpAnalysis');   % V2 code
    setup_chirpAnalysis_paths(true);
    assert(~isempty(which('tfridge')),          'tfridge (Signal Toolbox) not on path');
    assert(~isempty(which('ft_preproc_bandpassfilter')), 'FieldTrip firws not on path');
    assert(~isempty(which('cue_noise_trials')) || true, 'cue_noise_trials resolved in driver');

    setenv('CHIRP_DATAROOT','E:\Lab_Common');
    setenv('CHIRP_FIGROOT_E','E:\Lab_Common\Adam\Dupi_processing');  % E: mirror of subject figs
    setenv('CHIRP_V2OUT','E:\chirpV2out');

    maxNumCompThreads(6);
    chirpAnalysisV2([], struct('saveFig',true,'saveMat',true));

    fprintf('=== chirpV2 job DONE %s ===\n', datestr(now));
    fid=fopen('E:\_CHIRPV2_DONE.txt','w'); fprintf(fid,'done\n'); fclose(fid);
catch ME
    fprintf('ERROR: %s\n%s\n', ME.message, getReport(ME,'extended'));
    fid=fopen('E:\_CHIRPV2_ERROR.txt','w'); fprintf(fid,'%s\n',ME.message); fclose(fid);
end
diary off
