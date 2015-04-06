.\ssa -m .\models\examples\dimer_decay.xml -t1 -r10 -f --keep-histograms --out-dir ".\models\examples\ssa output"
"%~dp0tau_leaping" -m "%~dp0models\examples\dimer_decay.xml" -t1 -r10 -f --keep-histograms --out-dir "%~dp0models\examples\tau_leaping_output"
.\tools\SBMLconverter\bin\sbml2stochkit .\tools\SBMLconverter\hsr.xml
