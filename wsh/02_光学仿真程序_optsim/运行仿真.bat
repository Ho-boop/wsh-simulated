@echo off
chcp 65001 >nul
cd /d "%~dp0"
REM Double-click: run the example config.  Drag a .yaml file onto this .bat: run that config.
python optsim.py %*
pause
