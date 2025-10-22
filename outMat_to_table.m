function T = outMat_to_table(outMat, datalabel)
% outMat_to_table  Build a table from nested cell data with special handling:
% If the 5th value (label 'response') is empty, drop the 6th value
% (label 'response_str') and shift subsequent values left.
%
% T = outMat_to_table(outMat, datalabel)
%
% Inputs
%   outMat    : Nx1 cell array; each cell contains a 1x8 or 1x9 cell array.
%   datalabel : 1x9 cell array of labels (expects 'response' at 5, 'response_str' at 6).
%
% Output
%   T         : table with size [numel(outMat) x 9], columns named from datalabel.

    % Basic checks
    if ~iscell(outMat)
        error('outMat must be a cell array.');
    end
    if ~iscell(datalabel) || numel(datalabel) ~= 9
        error('datalabel must be a 1x9 cell array of labels.');
    end

    nRows = numel(outMat);
    nCols = 9;

    % Identify indices of interest (be strict about positions 5 and 6)
    idxResponse   = find(strcmpi(datalabel, 'response'), 1);
    idxResponseStr= find(strcmpi(datalabel, 'response_str'), 1);

    useSpecialRule = ~isempty(idxResponse) && ~isempty(idxResponseStr) ...
                     && idxResponse == 5 && idxResponseStr == 6;

    % Preallocate cell matrix for the table contents
    C = cell(nRows, nCols);

    for r = 1:nRows
        row = outMat{r};

        % Normalize to a 1xK cell row vector
        if isempty(row)
            row = cell(1,0);
        elseif ~iscell(row)
            row = {row};
        end
        row = row(:).';  % force row orientation

        % --- Special rule: if 5th ('response') is empty, skip 6th ('response_str') ---
        if useSpecialRule
            respVal = [];
            if numel(row) >= 5
                respVal = row{5};
            end
            if isEmptyValue(respVal)
                endVals = row(6:8); 
                row{5} = 999; 
                row{6} = []; 
                row(7:9) = endVals; 
               
            end
        end

        % Fill into output (sequentially), leave remaining cells blank
        k = min(numel(row), nCols);
        if k > 0
            C(r, 1:k) = row(1:k);
        end
    end

    % Make valid, unique variable names from datalabel
    varNames = matlab.lang.makeValidName(string(datalabel));
    varNames = matlab.lang.makeUniqueStrings(varNames, {}, namelengthmax);
    varNames = cellstr(varNames);

    % Convert to table
    T = cell2table(C, 'VariableNames', varNames);
end

function tf = isEmptyValue(v)
% Consider empty if:
% - empty cell/content
% - empty char ''
% - empty string ""
% - NaN numeric scalar
% - NaT datetime
    tf = isempty(v);
    if ~tf && ischar(v)
        tf = isempty(v); % '' counts as empty
    elseif ~tf && (isstring(v))
        tf = all(strlength(v) == 0);
    elseif ~tf && isnumeric(v) && isscalar(v)
        tf = isnan(v);
    elseif ~tf && isdatetime(v) && isscalar(v)
        tf = isnat(v);
    end
end
