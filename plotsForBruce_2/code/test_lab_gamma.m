function test_lab_gamma()
cd('E:\GitHub\ZelanoLabScripts\plotsForBruce_2\code');
addpath('E:\GitHub\ZelanoLabScripts'); addpath('E:\GitHub\Superlets\matlab-pure');
fid=fopen('E:\gammatest.txt','w');
try
  fp='R:\Neurology\Zelano_Lab\Lab_Common\Dupi\250623_Dupi_NMH_KS_1\preProc\250623_Dupi_NMH_KS_1_cueTaskPreproc.mat';
  fprintf(fid,'exist=%d\n', exist(fp,'file'));
  t0=tic; out=extract_gamma_session(fp,'cueTask','macBP1',struct());
  fprintf(fid,'LAB_OK time=%.1fs nBreaths=%d cols=%d coupMI=%.4f\n', toc(t0), out.nBreaths, width(out.perBreath), out.coupling.MI);
catch e
  fprintf(fid,'ERR %s\n', e.message);
end
fclose(fid);
end
