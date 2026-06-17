# taskList.md — execution plan

Work the tasks in order. Validate after each. Keep legacy
`getSessionParams_*.m` / `assemble_outDat_*.m` as oracles until the end.

---

## Task 0 — Orient

- [ ] Read `CLAUDE.md`.
- [ ] Read the 4 `getSessionParams_*.m`, the 3/4 `assemble_outDat_*.m`,
  `detect_ttls_O15.m`, the 3 `*_makeOutDat.m`, and the 4 `*_main.m`.
- [ ] Confirm `dataTracking.xlsx` has the parameter columns (schema in CLAUDE.md).
- [ ] Note the loader paths each `getSessionParams_*` uses to load `raw` (you will move
  these into `assemble_outDat_all`).

---

## Task 1 — `detectBeats.m`

Replaces the 25 bespoke `getBeats_*` functions. Each legacy detector had the shape:

```matlab
test = find(arrayfun(@(x,y,z) <conds>, ECGz(c1,a1:end-b1), ECGz(c2,a2:end-b2), ...));
test = test(diff(test) > beatSep);
heartBeats = test;
```

The channel slices are a **relative-lag alignment**. Encode each conjunct as
`chan,lag,op,thr` where `lag = a-1` (start sample minus 1), `op ∈ {gt,lt}`. Conjuncts are
joined by `&`. Example (legacy `…JH_1`): `1,2,gt,5 & 2,0,gt,4 & 3,0,lt,-0.5`.

Algorithm (reproduces the legacy slicing by truncating every channel to the common
length after applying its lag):

```matlab
function heartBeats = detectBeats(ECGz, beatSep, spec)
    conj = parseSpec(spec);              % fields: chan, lag, op, thr
    if isempty(conj), heartBeats = []; return; end
    N      = size(ECGz,2);
    maxLag = max([conj.lag]);
    L      = N - maxLag;
    if L <= 1, heartBeats = []; return; end
    mask = true(1,L);
    for i = 1:numel(conj)
        v = ECGz(conj(i).chan, 1+conj(i).lag : conj(i).lag+L);   % length L
        if strcmpi(conj(i).op,'gt'), mask = mask & (v >  conj(i).thr);
        else,                        mask = mask & (v <  conj(i).thr); end
    end
    test = find(mask);
    test = test(diff(test) > beatSep);   % NB: drops the last element by construction
    heartBeats = test;
end
```

- [ ] Implement `detectBeats.m` with a local `parseSpec` (`strsplit` on `&` then `,`).
- [ ] **Validate** against every distinct legacy `getBeats_*` body: build a synthetic
  `ECGz` (e.g. `randn(3,200000)`), pick a `beatSep`, and assert
  `isequal(detectBeats(ECGz,beatSep,spec), getBeats_<name>(ECGz,beatSep))` for all 25.
  The specs are already in the `beatSpec` column of `dataTracking.xlsx`; the legacy bodies
  are in `getSessionParams_breathingTask.m`. Fix `detectBeats` (not the data) on any mismatch.

---

## Task 2 — `applyParams.m`

Single source of truth. Default path
`R:\Neurology\Zelano_Lab\Lab_Common\Admin\dataTracking.xlsx`; optional 3rd arg overrides.

**Dispatch:** if `sel` is `'makeOutDat'`/`'main'` → Mode A (cfg); otherwise treat `sel` as
a session ID → Mode B (P).

**Sheet read (cache by path+mtime):** `readcell(xlsxPath,'Sheet','Sheet1')`; header row = 2;
data from row 3. Trim header names. Drop fully empty rows. Treat
`missing`/`NaN`/`''`/`[]` as blank.

**Helpers:**
- `canonTask(t)`: lower+strip spaces → `breathingtasks|wavebreathing`=`breathing`,
  `odorcuetask`=`cue`, `o15`=`O15`, `threshold`=`thresh`, else `''`.
- `taskKey('breathingTask'|'cueTask'|'threshTask'|'O15')` → internal canon.
- `studyOf(Type)`: contains `dupi`→`dupi`, `eeg`→`eeg`, else (`obecontrol`/`obe…`)→`obe`.
- `typeStr`: `dupi`→`Dupi`, `obe`→`OBE`, `eeg`→`EEG` (this is `P.type`).
- Coercion: `num_or` (numeric or `str2double`), `bool_or` (`true/1/yes`→true),
  `list_or("a,b,c")`→`str2num(['[' s ']'])`, blank→default.

**Row selection (both modes):** rows where `canonTask(Task)==taskKey(task)` AND raw
extracted (`Raw Data Extracted` non-blank and not `INCOMPLETE`). Dedupe by `sessID`
(case-insensitive, keep first). Stable order = sheet order.

**Mode A — `cfg`:**
- `cfg.sessionIDs` (n×1 cell, sheet IDs).
- `cfg.root` (n×1 cell): the row's `datPre` cell, else Type→root map.
- `cfg.datPre` (row cell): `{Dupi, OBEControl, EEGbreathing}` **always**, then any extra
  distinct roots appended. (Keeps 1/2/3 stable for the `datPrei==` branches.)
- `cfg.datPrei` (1×n): index of each row's root in `cfg.datPre`.
- `cfg.rspIDX`, `cfg.rspFlip` (1×n).
- `cfg.isNewStd` (1×n logical) from the `isNewStd` column.
- `cfg.newIDs` = `cfg.sessionIDs(cfg.isNewStd)`.
- `cfg.paramSource` (1×n cell).

**Mode B — `P`** (find the single matching row by sessID, case-insensitive; error if none):
- Common: `task` (the caller key, e.g. `'breathingTask'`), `type`, `fs_target=500`,
  `debug=false`, `computeResp=true`, `rspIDX`, `rspFlip`, `hasEEG`, `spikeClean`,
  `spikeThresh`, `spikeWin`, `macroRemove`, `paramSource`.
- `breathingTask`: `hasMacros`; `beatSpec` (default `'1,0,gt,3'`);
  `getBeats = @(ECGz,beatSep) detectBeats(ECGz,beatSep,beatSpec)`.
- `cueTask`: `respThresh`,`cuedBackBuff`,`adjWin`;
  `ttlMap = struct('cue',{'cue','Cue','cueOnset'},'target',{'targ','target','TargetOnset'},'resp',{'resp','response','button'})`.
- `threshTask`: `respThresh`,`cuedBackBuff`,`adjWin`; `ttlMap = struct('sniff',{'sniff'})`.
- `O15`: `respThresh`,`cuedBackBuff`,`adjWin`;
  `pd = struct('zthresh',-2,'minPulseSamp',200,'maxPulseSamp',1200,'trialSplitSamp',850)`;
  `ttl = struct('expectedTrialCount',30,'removeTrialMarksIdx',list_or(ttlRemoveIdx,[]),'note',string(ttlNote))`.

**Defaults for blank cells:** `rspIDX=1`, `rspFlip=1`, `spikeThresh=20`, `spikeWin=11`,
`macroRemove=[]`; `hasEEG=true` except O15 (`false`); `spikeClean=true` except O15
(`false`); `respThresh=500`, `cuedBackBuff=150`, `adjWin=500`; `hasMacros=true`.

- [ ] Implement `applyParams.m` with the locals above.
- [ ] **Validate P parity:** for each session in an old switch, compare
  `applyParams('<task>', sessID)` to `getSessionParams_<task>(S)` field-by-field
  (skip `raw`; for breathing compare `func2str(P.getBeats)` resolves to the same spec).

---

## Task 3 — `assemble_outDat_all.m`

`[outDat, raw, TTL] = assemble_outDat_all(S, P)` — combines the raw-loading half of
`getSessionParams_*` with the `assemble_outDat_*` assemblers, branching on `P.task`.

**Figure dir:** `figDir = S.fig` if set, else `fullfile(S.figPath, S.id)`. Set
`raw.paths.fig = figDir`; mkdir as needed.

**Raw load by task** (copy exactly from the legacy loaders):
- `breathingTask`: `<root>/<id>/preProc/<id>_breathingPreProc.mat` → `od` (try
  `outDat`→`out`→`chanDat`). Set `raw.sessID/fs_raw/data/labels/beh`. If `od` has no
  `TTL`, build the 5-min fallback `TTL=0:600000:size(data,2)` with `TTL(1)=1`, drop last;
  then `raw.TTL = round(od.TTL./4)`.
- `cueTask`: `<root>/<id>/preProc/<id>_cueTaskPreProc.mat` → `od.outDat`; `raw.*`;
  `raw.TTL = od.TTL` if present.
- `threshTask`: `<root>/<id>/preProc/<id>_PEA_threshold_preproc.mat` → same as cue.
- `O15`: `<root>/<id>/raw/raw_O15/raw_O15.mat` → `curDat`;
  `raw.fs_raw=curDat.rawData.fsample`, `raw.data=curDat.rawData.trial{1}`,
  `raw.labels=curDat.outLabs`, `raw.ncslabels` if present; behavior CSV
  `<root>/<id>/Behavioral_data/O15/O15_responses_<id>.csv` via `readtable`. Then
  `[TTL, raw] = detect_ttls_O15(raw, P);`.
- All: `raw.type = P.type`; `raw.paths.root = S.root`.

**Build `outDat`:** `behDat, labels, fs, data, sessID=S.id, task=P.task, type=P.type`,
`figs=fullfile(figDir,char(P.task))` (mkdir), `rspIDX=P.rspIDX`, `rspFlip=P.rspFlip`. For
breathing/cue/thresh: `if isfield(raw,'TTL'), outDat.TTL=raw.TTL; end`. For O15 add
`CSClist=raw.ncslabels`, `OGdataDir=fullfile(S.root,S.id)`, `loadFile` (the single
`*.m` in the session dir containing `'LoadData'`; error if not unique),
`preProcScript='O15PreProc.m'`. **Do not** set `outDat.TTL` for O15 (the main script does).

**Outputs:** `TTL` = O15 table from `detect_ttls_O15`, else `[]`.

- [ ] Implement `assemble_outDat_all.m`.
- [ ] **Validate:** for one session per task, compare `outDat` fields to the legacy
  `assemble_outDat_*`; for O15 compare `TTL` to `detect_ttls_O15(raw,P)` and confirm
  `OGdataDir`/`loadFile` resolve.

---

## Task 4 — Rewrite the 4 `*_main.m`

Surgical edits only; everything from the first `outDat = downsample_data(...)` onward is
unchanged.

For each of `breathingTaskPreProc_main.m`, `cueTaskPreProc_main.m`,
`threshPreProc_main.m`, `O15PreProc_main.m`:

- [ ] Delete the array definitions `datPre = {...}`, `datPrei = [...]`, `sessionIDs = {...}`.
- [ ] Immediately before the `for s = ...` loop insert:
  ```matlab
  cfg        = applyParams('<TASK>','main');
  sessionIDs = cfg.sessionIDs;
  ```
- [ ] Change the loop bound to `for s = 1:numel(sessionIDs)`.
- [ ] Replace `S.root = datPre{datPrei(s)};` with `S.root = cfg.root{s};`.
- [ ] Keep the skip-check (the `outDat = load(...)`/`isfield`/`continue` block, or O15's
  `if ~exist(... _O15preproc.mat)` wrapper) **verbatim**.
- [ ] Replace `[raw, P] = getSessionParams_<task>(S);` with
  `P = applyParams('<TASK>', S.id);`.
- [ ] Replace the assemble call:
  - breathing/cue/thresh: `outDat = assemble_outDat_breathing_cue_Task(raw, S, P);`
    (thresh has a trailing comment) → `[outDat, raw, TTL] = assemble_outDat_all(S, P);`.
  - O15: **delete** the line `[TTL, raw] = detect_ttls_O15(raw, P);`, and replace
    `outDat = assemble_outDat_O15(raw, S, P);` with
    `[outDat, raw, TTL] = assemble_outDat_all(S, P);`. Keep the following
    `outDat.TTL = TTL;` line.
- [ ] Leave `codePre`, addpaths, `figPath`, `EEGLOC`, `targTraceDir` (breathing), the
  `success` vector (breathing) and the `try/catch` (breathing) unchanged.

`<TASK>` = `breathingTask` / `cueTask` / `threshTask` / `O15`.

---

## Task 5 — Rewrite the 3 `*_makeOutDat.m`

Surgical edits only; the entire ingestion loop (per-session branches, `newList`/`newSet`
"standard" branch, `datPrei==1/2` type-setting, the save) is unchanged.

For each of `breathingTask_makeOutDat.m`, `cueTask_makeOutDat.m`,
`threshPreProc_makeOutDat.m`:

- [ ] Delete the array definitions: `datPre`, `datPrei`, `sessionIDs`,
  `newList` (breathing) / `newSet` (cue, thresh), `rspIDX`, `rspFlip`.
- [ ] Keep `codePre`, `behDatPath` (cue), `behDatPath_newSet` (thresh), all addpaths,
  `set(0,...)`, `ft_defaults` (thresh) unchanged.
- [ ] Immediately before the loop (`parfor sessi=...` for breathing, `for sessi=...`
  otherwise) insert:
  ```matlab
  cfg        = applyParams('<TASK>','makeOutDat');
  sessionIDs = cfg.sessionIDs;
  datPre     = cfg.datPre;
  datPrei    = cfg.datPrei;
  <newList|newSet> = cfg.newIDs;   % newList for breathing; newSet for cue/thresh
  rspIDX     = cfg.rspIDX;
  rspFlip    = cfg.rspFlip;
  ```
- [ ] Everything else unchanged.

Note: the breathing bespoke branch keyed on `'250904_OBE_NWU_TI'` will not match the
sheet ID `…_TI_1` and will fall through to the standard reader. This is the known TI/TI_1
issue — leave it; do not patch.

---

## Task 6 — Full validation & cleanup

- [ ] Run all parity checks from Tasks 1–3.
- [ ] End-to-end: for one already-processed session per task, run the rewritten
  `_makeOutDat` (where applicable) + `_main` and confirm it completes and matches a
  reference output where one exists.
- [ ] Adding a new participant should now be: fill the sheet row + parameter columns
  (`paramSource='guess'` rows pre-filled by carry-forward, to be reviewed) — no code edits.
- [ ] Only after everything passes, retire `getSessionParams_*.m` and
  `assemble_outDat_*.m` (or leave them parked). The 25 `getBeats_*` bodies can then be
  deleted along with `getSessionParams_breathingTask.m`.

## Out of scope (flagged earlier; do NOT change now)

- Filename case mismatches (`_breathingPreproc` vs `_breathingPreProc`), the `redoAGAIN`
  prefix in breathing `makeOutDat`, and the `250904_OBE_NWU_TI` vs `…_TI_1` naming.
