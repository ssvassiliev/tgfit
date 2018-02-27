# Global and target analysis of time-resolved spectroscopic data #

May 2013. 
Added Euler propagator.
Fixed bug: non-diagonal transfer matrix was not properly assembled 

Tgfit package now includes two programs:

	tgfit 
	tgfit-auto

tgfit  requires several input files:

1. Run control parameters [~/.tgfit/.tgfitrc],
2. Fitting model describing what and how to fit (answer file). 
    Example aswer files [sasX.ans] are in [~/.tgfit]
3. A set of datafiles
4. Steady-state fluorescence spectrum [./spec.dat]. This file 
   contains one column of numbers representing fuorescence intensity at all wavelength at which kinetics were measured.
tgfit performs least squares fit of the fluorescence decay kinetics.It takes starting parameter values from the answer file.

tgfit-auto is provided to facilitate search of the proper starting 
   parameters. In addition to the described above input files this 
   program takes a set of starting parameters from [~/.tgfit/auto-X.rc] 
   and performs fitting using all combinations of starting parameters.    

To install tgfit:

1. Install plplot:

http://sourceforge.net/projects/plplot/

2. Compile tgfit and tgfit-auto:

make tgfit
make auto

3. Make sure you have directory /<your-home>/bin
   if you are not sure, create it again:
      mkdir ~/bin

4. Perform installation:
   
make install
