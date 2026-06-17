function dumpSheet(xlsxPath, outTxt)
% Dump every cell of Sheet1 as tab-separated text so it can be read outside MATLAB.
if nargin < 1 || isempty(xlsxPath)
    xlsxPath = 'C:\Users\Adam\Documents\GitHub\ZelanoLabScripts\dataTracking.xlsx';
end
if nargin < 2 || isempty(outTxt)
    outTxt = 'C:\Users\Adam\Documents\GitHub\ZelanoLabScripts\_dev\sheetDump.txt';
end
C = readcell(xlsxPath, 'Sheet', 'Sheet1');
fid = fopen(outTxt, 'w');
[nr, nc] = size(C);
fprintf(fid, 'ROWS=%d COLS=%d\n', nr, nc);
for r = 1:nr
    parts = cell(1, nc);
    for c = 1:nc
        v = C{r,c};
        if isa(v,'missing') || (isnumeric(v) && isempty(v))
            parts{c} = '<NA>';
        elseif ischar(v)
            parts{c} = v;
        elseif isstring(v)
            parts{c} = char(v);
        elseif isnumeric(v) || islogical(v)
            parts{c} = num2str(v);
        else
            parts{c} = ['<' class(v) '>'];
        end
    end
    fprintf(fid, '%d\t%s\n', r, strjoin(parts, ' | '));
end
fclose(fid);
fprintf('wrote %s\n', outTxt);
end
