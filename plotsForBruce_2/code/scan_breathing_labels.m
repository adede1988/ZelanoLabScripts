function scan_breathing_labels()
% Enumerate the behDat task/condition label variants across all breathing finals
% (Dupi + OBE), to find labels the 'audio'/'focus' exact-match filter misses.
here=fileparts(mfilename('fullpath')); proj=fileparts(here);
idx=readtable(fullfile(proj,'out','tables','session_index.csv'),'TextType','string');
isB = contains(idx.task,'breathing','IgnoreCase',true) & ...
      (idx.cohort=="Dupi"|idx.cohort=="OBE") & idx.onDisk==1;
rows=idx(isB,:);
fprintf('scanning %d breathing finals\n\n', height(rows));
fid=fopen(fullfile(proj,'out','tables','breathing_label_scan.csv'),'w');
fprintf(fid,'sessID,cohort,taskColName,uniqueLabels,nAudioExact,nFocusExact,nAudioLike,nFocusLike,nTotal\n');
for i=1:height(rows)
    fp=char(rows.finalPath(i)); sid=char(rows.sessID(i)); coh=char(rows.cohort(i));
    if ~isfile(fp), fprintf('MISSING %s\n', fp); continue; end
    try
        s=load(fp); fn=fieldnames(s); od=s.(fn{1}); clear s;
        if ~isfield(od,'behDat'), fprintf('%s: no behDat\n',sid); clear od; continue; end
        bd=od.behDat; clear od;
        vn=bd.Properties.VariableNames;
        % candidate task columns (case-insensitive), prefer exact 'task'
        cand=vn(~cellfun(@isempty, regexpi(vn,'^(task|condition|cond|block|stim|blockType|type)$','once')));
        if isempty(cand), cand=vn(contains(lower(vn),'task')); end
        tcName='<none>'; tc=strings(0);
        if ~isempty(cand)
            tcName=cand{1};
            col=bd.(tcName);
            if iscell(col), tc=string(col); elseif iscategorical(col), tc=string(col);
            elseif isstring(col), tc=col; elseif ischar(col), tc=string(cellstr(col));
            else, tc=string(col); end
        end
        tc=strtrim(tc);
        u=unique(tc); u(u=="")=[];
        ulist=strjoin(cellstr(u),' | ');
        nAudioX=sum(strcmpi(tc,'audio')); nFocusX=sum(strcmpi(tc,'focus'));
        nAudioL=sum(contains(lower(tc),'audio')); nFocusL=sum(contains(lower(tc),'focus'));
        fprintf('%-28s %-5s col=%-12s | %s\n', sid, coh, tcName, ulist);
        fprintf(fid,'%s,%s,%s,"%s",%d,%d,%d,%d,%d\n', sid,coh,tcName,ulist,nAudioX,nFocusX,nAudioL,nFocusL,numel(tc));
        clear bd tc u;
    catch e
        fprintf('FAIL %s: %s\n', sid, e.message);
        fprintf(fid,'%s,%s,ERROR,"%s",0,0,0,0,0\n', sid,coh,e.message);
    end
end
fclose(fid);
fprintf('\nwrote out/tables/breathing_label_scan.csv\n');
end
