function reproc_markDone(p, id)
% Append a completed session ID to a done-marker file (for resume-on-restart).
    fid = fopen(p, 'a');
    if fid > 0, fprintf(fid, '%s\n', id); fclose(fid); end
end
