function integrate_new_controls()
% Integrate the two additional OBE control breathing sessions that use the
% 'breathingTasks_separatepreproc.mat' naming AND label blocks 'audiobook' /
% 'focusedBreathing' (HM_2, SP_2) — both were missed by the pipeline. Computes
% their best macBP gamma channel (same FOOOF-flattened-peak rule as run_scores)
% and writes their spectro2 band-z maps for both conditions, plus a macbp_best row.
proj='C:/Users/Adam/Documents/GitHub/ZelanoLabScripts/plotsForBruce_2'; code=fullfile(proj,'code');
addpath(code); addpath('C:/Users/Adam/Documents/GitHub/ZelanoLabScripts');
addpath('C:/Users/Adam/Documents/GitHub/Superlets/matlab-pure');
H=lib_scores(); BAND=[25 58]; s2dir=fullfile(proj,'out','spectro2');
base='R:/Neurology/Zelano_Lab/Lab_Common/OBEControl/';
sess={'260625_OBE_NWU_HM_2','260702_OBE_NWU_SP_2'};
bestRows={};
for i=1:numel(sess)
  id=sess{i}; fp=fullfile(base,id,'preProc',[id '_breathingTasks_separatepreproc.mat']);
  od=[]; for a=1:3, try, X=load(fp); vv=fieldnames(X); od=X.(vv{1}); clear X; break; catch e, fprintf('[retry %s] ',e.message); pause(3); end, end
  if isempty(od), fprintf('%s LOAD FAIL\n', id); continue; end
  fs=500; if isfield(od,'fs'), fs=double(od.fs); end
  labs=cellfun(@(x)char(string(x)),od.labels,'uni',0);
  isMac=cellfun(@(x)~isempty(regexpi(x,'macbp','once')),labs); macIdx=find(isMac); nMac=numel(macIdx);
  pf=nan(nMac,1); det=false(nMac,1);
  for m=1:nMac, Rm=H.macbp_gamma(double(od.data(macIdx(m),:)), fs, BAND); pf(m)=Rm.peakFlatDb; det(m)=Rm.gammaDetected; end
  [~,bi]=max(pf); bl=labs{macIdx(bi)};
  fprintf('%s: bestMac=%s (%.2f dB, det=%d of %d)\n', id, bl, pf(bi), sum(det), nMac);
  bestRows{end+1}=sprintf('%s,breathingTask,OBE,Control,%s,%d,%d,%d,%s,%.6g,%d,%.6g,%d,%d', ...
     id, subj_of(id), sn_of(id), nMac, bi, bl, pf(bi), det(bi), 0, sum(det)-double(det(bi)), sum(det)); %#ok<AGROW>
  for tr={'audiobook','focusedBreathing'}
    try
      out=extract_spectro_session(fp, tr{1}, bl, struct()); %#ok<NASGU>
      save(fullfile(s2dir,[id '__' tr{1} '.mat']),'out','-v7');
      fprintf('   %s nBreaths=%d saved\n', tr{1}, out.nBreaths);
    catch e, fprintf('   %s FAIL: %s\n', tr{1}, e.message); end
  end
end
% append to macbp_best.csv
bf=fullfile(proj,'out','tables','macbp_best.csv');
fid=fopen(bf,'a'); for i=1:numel(bestRows), fprintf(fid,'%s\n', bestRows{i}); end; fclose(fid);
fprintf('DONE integrate_new_controls (appended %d macbp_best rows)\n', numel(bestRows));
end
function s=subj_of(id), p=strsplit(id,'_'); s=p{4}; end
function n=sn_of(id), p=strsplit(id,'_'); n=str2double(p{end}); if isnan(n), n=1; end, end
