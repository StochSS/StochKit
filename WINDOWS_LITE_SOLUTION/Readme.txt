This solution is used to generate the executables for windows lite version stochkit. I assume you have Visual Studio Professional.
To build the solution
1. Copy the WINDOWS_LITE_SOLUTION folder to a windows machine
2. Put a copy of source files into this folder
   -copy the StochKit "src" folder into the WINDOWS_LITE_SOLUTION folder
   -copy output.h and readSBML.cpp from the StochKit tools/SBMLconverter directory into the WINDOWS_LITE_SOLUTION/tools/SBMLconverter folder
3. Open AbsoluteDirectories.vsprops and edit the boost paths
(3.5. If you only have visual studio express, open StochKit.sln with visual studio and convert the solution to your version. If you use non-express visual studio, skip this step)
4. Open a visual studio command prompt.
5. Run the Runme.bat in the command prompt (you can do this by dragging the batch file into the command window and clicking enter, or by navigating to the folder where the batch file is and typing "Runme").
6. Runme.bat will create a folder named WINDOWS_LITE. Upload this folder into the svn repository.
7. The package_StochKit2 script will copy the Matlab tools and model files before release.
