.\ssa -m .\models\examples\dimer_decay.xml -t1 -r10 -f --keep-histograms --out-dir .\models\examples\dimer_decay_test

"%~dp0tau_leaping" -m .\models\examples\schlogl.xml -t1 -r10 -f --keep-histograms --out-dir ".\models\examples\schlogl output"

"%~dp0ssa" -m "%~dp0models\examples\events1.xml" -t10 -r10 -f --keep-histograms --out-dir "%~dp0.\models\examples\events1_test"

.\tools\SBMLconverter\bin\sbml2stochkit .\tools\SBMLconverter\hsr.xml
