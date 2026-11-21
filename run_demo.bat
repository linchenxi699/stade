@echo off
chcp 65001 >nul
setlocal
set "ROOT=%~dp0"
set "PY=D:\coding enviroment\python\python.exe"
cd /d "%ROOT%"

if not exist "%PY%" (
  echo [ERROR] Python not found at %PY%
  echo Edit this file and set PY to your python.exe path.
  pause
  exit /b 1
)

echo ============================================
echo  StaDE demo: simulate -^> DE -^> stability -^> enrichment
echo ============================================
echo.

echo [1/3] generating the demo dataset ...
"%PY%" scripts\99_make_demo_dataset.py
if errorlevel 1 goto :fail

echo.
echo [2/3] preprocessing + differential expression + 200 bootstraps ...
"%PY%" -m stade.cli prep --counts data\demo\demo_counts.tsv.gz --meta data\demo\demo_sample_meta.tsv --meta-index sample_id --group-col group --case-col Tumor --tmm --outdir demo_work
if errorlevel 1 goto :fail
"%PY%" -m stade.cli de --bundle demo_work\bundle.npz --case-col Tumor --n-boot 200 --outdir demo_out
if errorlevel 1 goto :fail

echo.
echo [3/3] KEGG enrichment ...
"%PY%" -m stade.cli enrich --de demo_out\de_stability.tsv --genesets data\raw\genesets\KEGG_2021_Human.txt --library KEGG_demo --robust-only --n-perm 500 --outdir demo_out\enrichment
if errorlevel 1 goto :fail

echo.
echo ============================================
echo  DONE. Results are in demo_out:
echo    demo_out\de_stability.tsv        DE + stability table
echo    demo_out\figs\                   volcano, PCA, stability plots
echo    demo_out\enrichment\             enrichment tables and dot plot
echo ============================================
start "" explorer "%ROOT%demo_out"
pause
exit /b 0

:fail
echo.
echo [FAILED] copy the error above and send it to the assistant.
pause
exit /b 1
