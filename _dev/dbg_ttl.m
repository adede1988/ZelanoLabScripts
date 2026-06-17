function dbg_ttl()
    repo = 'C:\Users\Adam\Documents\GitHub\ZelanoLabScripts';
    addpath(repo); xlsx = fullfile(repo,'dataTracking.xlsx');
    id = '250623_Dupi_NMH_KS_1';
    root = 'R:\Neurology\Zelano_Lab\Lab_Common\Dupi\';
    S = struct('id', id, 'root', root, 'figPath', tempdir);
    if ~isfolder(fullfile(tempdir,id)), mkdir(fullfile(tempdir,id)); end

    [rawL, PL] = getSessionParams_O15(S);
    PN = applyParams('O15', id, xlsx);

    fprintf('--- P field diffs (PL vs PN) ---\n');
    f = union(fieldnames(PL), fieldnames(PN));
    for i=1:numel(f)
        fn=f{i}; a=[]; b=[];
        if isfield(PL,fn), a=PL.(fn); end
        if isfield(PN,fn), b=PN.(fn); end
        if ~isequaln(a,b)
            fprintf('  P.%s differs\n', fn);
            if isstruct(a)&&isstruct(b)
                sf=union(fieldnames(a),fieldnames(b));
                for j=1:numel(sf)
                    va=[];vb=[];
                    if isfield(a,sf{j}),va=a.(sf{j});end
                    if isfield(b,sf{j}),vb=b.(sf{j});end
                    if ~isequaln(va,vb), fprintf('     .%s : %s vs %s\n', sf{j}, mat2str(va), mat2str(vb)); end
                end
            else
                disp(a); disp(b);
            end
        end
    end

    % run detect on fresh raw copies
    rawA = rawL;  % already has data
    [TTLa, ~] = detect_ttls_O15(rawA, PL);
    rawB = rawL;
    [TTLb, ~] = detect_ttls_O15(rawB, PN);
    fprintf('TTL equal: %d\n', isequal(TTLa,TTLb));
    if ~isequal(TTLa,TTLb)
        Da = table2array(TTLa); Db = table2array(TTLb);
        d = find(~(Da==Db | (isnan(Da)&isnan(Db))));
        fprintf('num differing cells: %d of %d\n', numel(d), numel(Da));
        [r,c]=ind2sub(size(Da), d(1:min(5,end)));
        for k=1:numel(r)
            fprintf('  (%d,%d): %g vs %g\n', r(k),c(k), Da(r(k),c(k)), Db(r(k),c(k)));
        end
    end
end
