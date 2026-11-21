@echo off
chcp 65001 >nul
setlocal
set "PY=D:\coding enviroment\python\python.exe"

if not exist "%PY%" (
  echo [ERROR] Python not found at %PY%
  echo Edit this file and set PY to your python.exe path.
  pause
  exit /b 1
)

if "%~1"=="" (
  "%PY%" -m stade.cli --help
  echo.
  echo Example:
  echo   stade.bat prep --counts data.tsv --meta meta.tsv --group-col group --case-col Tumor --outdir work
  pause
  exit /b 0
)

"%PY%" -m stade.cli %*
if errorlevel 1 pause
endlocal
