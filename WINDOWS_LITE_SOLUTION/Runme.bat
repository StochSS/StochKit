@devenv /upgrade "%~dp0Stochkit.sln"
SET MSBuildUseNoSolutionCache=1
@msbuild "%~dp0Stochkit.sln" /clp:NoSummary /p:configuration=release;DebugSymbols=false;DebugType=None /v:m /nologo /m
@msbuild "%~dp0tools\SBMLconverter\SBMLconverter.vcproj" /clp:NoSummary /p:configuration=release;DebugSymbols=false;DebugType=None /v:m /nologo /m
@xcopy "%~dp0bin" "%~dp0WINDOWS_LITE\bin\" /q >nul
@xcopy "%~dp0*.exe" "%~dp0WINDOWS_LITE\" /q >nul
@xcopy "%~dp0tools\SBMLconverter\bin" "%~dp0WINDOWS_LITE\tools\SBMLconverter\bin\" /q >nul
@copy "%~dp0test.bat" "%~dp0WINDOWS_LITE\" >nul
@copy "%~dp0tools\SBMLconverter\documentation.doc" "%~dp0WINDOWS_LITE\tools\SBMLconverter\" >nul
@copy "%~dp0tools\SBMLconverter\documentation.pdf" "%~dp0WINDOWS_LITE\tools\SBMLconverter\" >nul
@copy "%~dp0tools\SBMLconverter\hsr.xml" "%~dp0WINDOWS_LITE\tools\SBMLconverter\" >nul
@echo.
@echo Done. The WINDOWS_LITE folder is what you want
