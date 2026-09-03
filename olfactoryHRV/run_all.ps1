# olfactoryHRV - full pipeline, extraction through report.
# Run from anywhere:  powershell -ExecutionPolicy Bypass -File run_all.ps1
# Paths come from ohrv_config.m (MATLAB) and the script location (Python), so
# nothing here is machine-specific except the two interpreter paths below.

$code = Split-Path -Parent $MyInvocation.MyCommand.Path
$work = Join-Path $code "work"
$ml   = "C:\Program Files\MATLAB\R2026a\bin\matlab.exe"
$py   = "C:\Program Files\PsychoPy\python.exe"     # any python with h5py + numpy

if (-not (Test-Path $ml)) { throw "MATLAB not found at $ml - edit this script" }
if (-not (Test-Path $py)) { throw "Python not found at $py - edit this script" }
New-Item -ItemType Directory -Force -Path $work | Out-Null

function Step($name, $block) {
  "=============================================================="
  "STEP $name   $(Get-Date -Format 'HH:mm:ss')"
  "=============================================================="
  & $block
  "-- $name exit=$LASTEXITCODE  $(Get-Date -Format 'HH:mm:ss')"
}

# Clear cached intermediates. The preprocessed .mat finals get rewritten by the
# preprocessing pipeline; stale extracts silently mix data versions, so always
# start clean rather than reusing what is on disk.
if ($args -notcontains "-keep") {
  Remove-Item "$work\*_slim.mat","$work\*_beats.npz","$work\*.csv","$work\panel_state.mat" `
              -Force -ErrorAction SilentlyContinue
  "cleared cache in $work"
}

Step "1 extract breathing (behDat + RRint)" { & $ml -batch "cd('$code'); rsa_extract" }
Step "2 extract heart beats"                { & $py "$code\get_beats.py" }
Step "3 score olfaction"                    { & $ml -batch "cd('$code'); rsa_olf" }
Step "4 per-session slopes + nulls"         { & $ml -batch "cd('$code'); rsa_analyze" }
Step "5 vagal metric panel"                 { & $ml -batch "cd('$code'); rsa_panel" }
Step "6 change scores + inference"          { & $ml -batch "cd('$code'); rsa_report" }
Step "7 threshold sensitivity"              { & $ml -batch "cd('$code'); thresh_analysis" }
Step "8 figures"                            { & $ml -batch "cd('$code'); rsa_figs; fig_grant; fig_grant_thresh; fig_duration; fig_raw" }
Step "9 build reports"                      { & $py "$code\mkreport.py"; & $py "$code\mkgrant.py" }

"=============================================================="
"ALL DONE $(Get-Date -Format 'HH:mm:ss')"
Get-ChildItem $work -Filter "*.csv" -ErrorAction SilentlyContinue |
  ForEach-Object { "  {0,-28} {1,7:N0} bytes" -f $_.Name, $_.Length }
