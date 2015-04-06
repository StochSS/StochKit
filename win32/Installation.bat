@echo off
echo.
echo Update the Solutions
echo.
cd /d %~dp0
devenv /upgrade Stochkit.sln 2>nul
devenv /upgrade Mixed_Compiled_Solution\Mixed_Compiled_Solution.sln 2>nul
devenv /upgrade custom_drivers\custom_drivers.sln 2>nul
echo.
echo Build Boost library
echo.
cd libs\boost_1_53_0
call bootstrap.bat
b2 variant=release --with-filesystem --with-program_options --with-thread --with-date_time --with-serialization --with-system --with-chrono --stagedir=.\
cd ..\..

echo.
echo Build Stochkit
echo.
SET MSBuildUseNoSolutionCache=1

msbuild Stochkit.sln /clp:NoSummary /p:configuration=release;DebugSymbols=false;DebugType=None /v:m /nologo /m

echo.

echo Finished
@echo on