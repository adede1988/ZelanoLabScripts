% run_chirpV2_thresh_lab.m  -- WMI-detached chirpAnalysisV2 threshTask batch for the lab E: session.
% Log E:\chirpV2_thresh.log; sentinel E:\_CHIRPV2_THRESH_DONE.txt / _CHIRPV2_THRESH_ERROR.txt.
diary('E:\chirpV2_thresh.log');
try
    fprintf('=== chirpV2 THRESH job start %s ===\n', datestr(now));
    addpath('E:\GitHub\ZelanoLabScripts'); addpath('E:\'); addpath('E:\GitHub\ZelanoLabScripts\chirpAnalysis');
    setup_chirpAnalysis_paths(true);
    setenv('CHIRP_DATAROOT','E:\Lab_Common');
    setenv('CHIRP_FIGROOT_E','E:\Lab_Common\Adam\Dupi_processing');
    setenv('CHIRP_V2OUT','E:\chirpV2out');
    maxNumCompThreads(6);
    chirpAnalysisV2('threshTask', [], struct('saveFig',true,'saveMat',true));
    fprintf('=== chirpV2 THRESH job DONE %s ===\n', datestr(now));
    fid=fopen('E:\_CHIRPV2_THRESH_DONE.txt','w'); fprintf(fid,'done\n'); fclose(fid);
catch ME
    fprintf('ERROR: %s\n%s\n', ME.message, getReport(ME,'extended'));
    fid=fopen('E:\_CHIRPV2_THRESH_ERROR.txt','w'); fprintf(fid,'%s\n',ME.message); fclose(fid);
end
diary off
