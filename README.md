# Rho
Application to quickly process alternating voltage/emf values to calculate resistivity and conductivity

The .exe file (in releases) should be downloaded by the user. 
If the user wishes to alter the source code, then download the .mlapp file which opens in MatLab App Designer. 

The purpose of this app is to facilitate electrical resistvity calculations by allowing the user to input 
the parameters that vary from one experiment to the other, such as the sample length and diameter, current 
and thermocouple type. Users also have the hand to remove irrelevant data from the raw data analysis, and to 
adjust the quality of the fit by playing with the treshold number. 

For feedback or questions, please email me at mberrada@uwo.ca

1. Prerequisites for Deployment 

Verify that version 9.11 (R2021b) of the MATLAB Runtime is installed.   
If not, you can run the MATLAB Runtime installer.
To find its location, enter
  
    >>mcrinstaller
      
at the MATLAB prompt.
NOTE: You will need administrator rights to run the MATLAB Runtime installer. 

Alternatively, download and install the Windows version of the MATLAB Runtime for R2021b 
from the following link on the MathWorks website:

    https://www.mathworks.com/products/compiler/mcr/index.html
   
For more information about the MATLAB Runtime and the MATLAB Runtime installer, see 
"Distribute Applications" in the MATLAB Compiler documentation  
in the MathWorks Documentation Center.

2. Files to Deploy and Package

Files to Package for Standalone 
================================
-Rho.exe
-MCRInstaller.exe 
    Note: if end users are unable to download the MATLAB Runtime using the
    instructions in the previous section, include it when building your 
    component by clicking the "Runtime included in package" link in the
    Deployment Tool.
-This readme file 

3. Definitions

For information on deployment terminology, go to
https://www.mathworks.com/help and select MATLAB Compiler >
Getting Started > About Application Deployment >
Deployment Product Terms in the MathWorks Documentation
Center.
