function pctcheck()
fid=fopen('E:\pct.txt','w');
fprintf(fid,'PCT_license=%d numcores=%d\n', license('test','Distrib_Computing_Toolbox'), feature('numcores'));
try
  p=parpool('Processes',4);
  fprintf(fid,'pool started: %d workers\n', p.NumWorkers);
  r=zeros(1,8); parfor i=1:8, r(i)=i^2; end
  fprintf(fid,'parfor OK sum=%d\n', sum(r));
  delete(p);
catch e
  fprintf(fid,'parpool ERR: %s\n', e.message);
end
fclose(fid);
end
