---
title: 'Rho: Application to analyse electrical resistivity'
tags:
  - Electrical Resistivity
  - Thermal Conductivity
  - Data analysis
  - Voltage drop
  
authors:
  - name: Meryem Berrada
    orcid: 0000-0002-4934-0917
    equal-contrib: true
    affiliation: 1 
  - name: Richard A. Secco
    orchid: 0000-0001-5029-659X
    equal-contrib: true 
    affiliation: 1

affiliations:
 - name: Department of Earth Sciences, University of Western Ontario, London, ON, Canada 
   index: 1

date: 25 July 2022
bibliography: paper.bib

# Summary

Electrical resistivity measurements of metals often use the four-wire method, where data acquisition alternates between temperature and voltage. The temperature and voltage values are similar and manual selection of each set of measurements is time consuming. `Rho` is intended to make electrical resistivity data analysis faster and more reliable by automatically selecting the temperature and voltage measurements using signal processing. The time required to process 52 data points (which consists of approximately 6000 measurements) using `Rho` is less than one minute compared to up to 4 hours for manual processing. `Rho` is stored in FigShare and can be used by any research group conducting alternate temperature and voltage measurements to obtain electrical resistivity at either ambient or extreme conditions.

# Introduction 
Electrical resistivity of materials is a basic property that is measured and used in a large variety of studies. In one type of study, electrical resistivity data are used to calculate thermal conductivity and interior heat flow of terrestrial-type bodies. The electrical resistivity of a sample is often obtained using the four-wire method as it reduces the contribution of the electrodes that is present in the two-wire method. The four-wire method is used in the cooling of industrial devices [@Hapenciuc:2019], soldering processes [@Tillmann:2017], electronic systems [@Krishnamurthy:2019], temperature determination at depth [@Erkan:2017], and simulations of physical properties of terrestrial-type bodies at high pressures and high temperatures. High pressure and high temperature measurements of electrical resistivity (ρ) in high pressure multi-anvil presses are used to calculate thermal conductivity (k) and the heat flow in the interiors of terrestrial-type bodies. Ezenwa and Secco [@Ezenwa:2017a] developed a four-wire method to solve the voltage bias, caused by passing a test current in a single direction, which previously resulted in erroneous analysis of ρ measurements due to parasitic voltages caused by the Seebeck effect or other sources. This polarity switching method allowed successful measurements of ρ of [@Ezenwa:2017a], Nb [@Ezenwa:2017b], Ni [@Silber:2017], Co [@Ezenwa:2017c], Cu [@Ezenwaetal:2017], Fe [@Silber:2018][@Yong:2019], Ag [@Littleton:2018], W [@Littleton:2019], Re [@Littleton:2019], Au [@Berrada:2018], Fe-Si alloys [@Silber:2019][@Berrada:2020][@Berrada:2021], and Fe-S alloys [@Littleton:2021]. Three other papers from different research groups report using similar data collection practices, while several other papers use the four-wire method without a polarity switch [@BerradaSecco:2021]. In these measurements, the sample is located in the middle of a high-pressure cell, between two metal discs which are each contacted by a thermocouple. The metal discs, ensuring contact between the thermocouples and the sample, are typically composed of W. The thermocouples are typically Type-C thermocouples (one leg is made of 95%W/5%Re and the other leg is made of 74%W/26%Re) where voltage (i.e., thermoelectric emf) is measured and then converted to temperature after the experiment. The voltage drop across the sample for resistivity measurement is made between a pair of wires on opposite sides of the sample (one wire leg from each thermocouple), and therefore includes the voltage contribution of the discs. Thus, the ρ of the sample is obtained by subtracting the disc contribution from the total measurement.

The four-wire method incorporates a polarity switch and a mode switch, and results in a pattern of alternating values of temperature (emf) and voltage in a single column when acquired by a single voltmeter. A complete measurement produces a pattern similar to the following: temperature emf before voltage measurement (T~b~), positive voltage drop with polarity in one direction (V^+^), negative voltage drop with polarity in opposite direction (V^-^), temperature emf after voltage measurement (T~a~). The order of V^+^ and V^-^ may be alternated depending on the position of the polarity switch. At the University of Western Ontario, in the High-Pressure High-Temperature Laboratory, a programmable Keysight B2961 power supply is used to provide a constant direct test current of 0.2 A while data are acquired by a programmable Keysight 34470A meter operating at 20 Hz and 1 µV resolution. The application `Rho` calculates the average of 10 data points for each temperature section and all available data points for each voltage section. It is expected that temperature remains stable during these 10 data points. At high temperatures, when the frequency of measurements is increased and the quantity of data points in each section is low (due to the rapid switching between each section), `Rho` selects all the available data for each section. The values in each selection are averaged. Then the V^+^ and V^-^ selections are averaged to give the voltage drop (V), while T~b~ and T~a~ selections are averaged to give the temperature (T). The voltage drop and current (I) are then used to calculate the electrical resistance (R ) of the sample using Ohm’s law: 
$$R=\frac{V}{I}$$	(1)
The electrical resistivity is then calculated using Pouillet’s law:
$$\rho=\frac{RA}{l}$$	(2)
where A is the sample’s cross-sectional area and l is the sample’s length which are both measured on the post-experiment recovered sample. Values of ρ are then used to calculate the electronic component of thermal conductivity (k) via the empirical Wiedemann-Franz law:
$$k=\frac{L_0T}{\rho_{total}}$$	(3)
where T is temperature and Lo is the theoretical Sommerfeld value (L~o~ = 2.44∙10^-8^ WΩK^-2^) of the Lorenz number. The error in temperature measurement corresponds to the standard deviation of the data points of each temperature section. The error in ρ is obtained by error propagation using the uncertainties in sample geometry and standard deviation of the voltage measurements. Similarly, the error in k corresponds to the propagation of the error in ρ. The selection of temperature and voltage values is typically done manually, and, to our knowledge, no software available in the literature is capable of automatically processing the temperature and voltage signals. 

# Implementation and architecture
The main page of `Rho` is displayed in \autoref{Fig.1}(a). The application requires data to be imported in .xlsx, .csv, or .xls format, where the T~b~, V^+^, V^-^, and T~a~ data alternate within the first column of the file (see \autoref{Fig.1}(b)). The first figure (see \autoref{Fig.1}(a)) that is generated when the file is loaded is of the raw data. Noise, or abrupt fluctuation in data, that are not representative of the measurements can be removed by selecting the starting and ending indices (x-axis) of the section. The application attributes NaN values to each data point of this section, which removes them from the future analysis. When ‘Plot Data’ is selected, `Rho` identifies the temperature and voltage sections based on deviations from the main increasing trend, which corresponds to the temperature increase with time. 

![(a) Main screen. The input parameters options are displayed on the left panel. A visualization of the Data selection, resulting resistivity and output results are displayed on the right. An example raw data set of Fe at 8 GPa is displayed in Figure 1 [@Berrada:2020]. (b) Example of raw data Excel file for Fe at 8 GPa [@Berrada:2020] where all measurements are within the first column of the file. The identification of parameter sequence is shown in the right column for the first two sets of data for clarity here only (this not included in the original data file).\label{Fig.1}](fig1.png)

First, the intervals of negative values are labelled as voltage measurements (V~total~). This temporary selection is displayed in Figure 1 of the application (see \autoref{Fig.1} below). In the first scenario, the pattern T~b~ V^+^ V^-^ T~a~ is considered and T~b~ and T~a~ are temporarily identified. To satisfy this scenario, the difference between T~b~ and V~total~ must be smaller than the difference between T~a~ and V~total~ considering that the part of V~total~ closer to T~b~ would be V^+^ and that closer to T~a~ would be V^-^ (and V^-^ < V^+^). Rho uses a loop to evaluate this scenario at each deviation from the main increasing trend. Whenever the scenario is observed, V^+^ is temporarily defined as the first trend in V~total~. In the second scenario, the pattern T~b~ V^-^ V^+^ T~a~ is considered, and a similar loop is applied to identify the intervals where the difference between T~b~ and V~total~ is larger than that of T~a~ and V~total~. Here, the part of V~total~ closer to T~b~ would be V^-^ and that closer to T~a~ would be V^+^. Another analysis of the first scenario identifies the stable intervals on the main increasing trend right before V^+^ as T~b~, and that right after V^-^ as T~a~. In other words, fluctuations, or noise, between temperature and voltage measurements are automatically ignored. A similar analysis of the second scenario identifies the stable intervals on the main increasing trend right before V^-^ as T~b~ and that right after V^+^ as T~a~. The corresponding indexes (x-axis) and data (y-axis) for each selection are combined into variables for T~b~, V^-^, V^+^, and T~a~. The ‘Outliers Degrees of Freedom’ parameter evaluates the differences between the trends of V^-^ and V^+^ and identifies data points outside of these trends as outliers. A larger degree of freedom will include more data points in the final variables for V^-^ and V^+^, while a smaller degree of freedom will ignore more data points. The value 20 seems to be the best fit for this parameter, although it can be adjusted by the user. In other words, fluctuations of 20 times and more are defined as outliers. The final selection is displayed in Figure 2 of the application (see \autoref{Fig.2}(a) and (b)). The contribution of the discs to the measured ρ is subtracted by fitting the ρ of the selected material. The ‘Discs’ dropdown menu contains options of W from 2-5 GPa [@Littleton:2019], Re from 2-5 GPa [@Littleton:2019] and Pt at 1 atm, 10 and 20 GPa [@Gomi:2019]. Supplementary options may be added to future updates of the application as data become available in the literature for different metals used for the discs. Temperature is automatically converted from emf to (K) using the selected thermocouple calibration equation. The ‘Thermocouple’ dropdown menu contains options for Type-C [@Omega:2019a], Type-S at 1 atm [@Omega:2019b] and Type-S from 1-5 GPa thermocouples [@Getting:1970]. The input parameters are then used in equation (1), (2) and (3) to plot ρ in Figure 3 of the application and output values of T (K), ρ (µΩ·cm), error in ρ (µΩ·cm), k (Wm^-1^K^-1^) and error in k (Wm^-1^K^-1^) (see \autoref{Fig.2}(b) and (c)). The values are displayed in the ‘Output’ tab and can be copied and pasted on another platform for further analysis.

![(a) Temporary and final selections of the temperature and voltage sections. (b) Calculated resistivity compared with that of Fe at 1 atm [@Chu:1981]. (c) Output values of T (K), ρ (µΩ·cm), error in ρ (µΩ·cm), k (Wm^-1^K^-1^) and error in k (Wm^-1^K^-1^). (d) Example of the final selections of one measurement for each parameter.\label{Fig.2}](fig2.png)

Examples of running the software are displayed in \autoref{Fig.2}. After loading the raw data file, Figure 1 of the application should automatically display the raw data. If the raw data file is in the incorrect format, an error sound will be set, and Figure 1 will not display the raw data. The users should confirm the temperature and voltage selections since similarity in these values might result in the incorrect identification of some data points. The computed ρ can be compared with that of Fe at 1 atm [@Chu:1981].

# Availability
Users may choose to download Rho.exe file from FigShare or Rho.m from GitHub. 
Rho.exe – Any system capable of running a .exe file. 
Rho.m – Any system capable of running Matlab R2012a or higher. 

Tested in Matlab R2012a, likely to work in earlier versions also. 


Name: FigShare Archive
Persistent identifier: 10.6084/m9.figshare.19175432  
Licence: CC by 4.0
Publisher: Meryem Berrada
Version published: 1.0
Date published: 15/02/22

Name: GitHub
Identifier: https://github.com/meryemberradauwo/Rho.git 
Licence: GNU General Public License v3.0
Date published: 02/02/22

# Statement of need
The application Rho minimizes the data analysis time and uses signal processing to identify the temperature (emf) and voltage drop measurements from a single column of alternating measurements resulting from the four-wire method. The confirmation of the selection by the user remains required considering the limitations of the attributed thresholds for outliers. This software (Rho.exe) can be used to analyse electrical resistivity for applications to the cooling of industrial devices, soldering processes, electronic systems, temperature determination at depth, and simulations of physical properties of terrestrial-type bodies at high pressures and high temperatures. Rho can be used by any research group conducting alternate temperature (emf) and voltage measurements to obtain electrical resistivity at either ambient or extreme conditions. Modifications to the software, such as changing the identification thresholds, adding TC types or disc materials, can be done through Rho.m via the MATLAB platform. The authors are continuously working to improve the software and implement more features. 

We encourage the reader to download Rho via MATLAB [@Berrada:2022] or through this article. Inquiries can be addressed to the corresponding author of this article, or through the contact information displayed in the software.

# Acknowledgements
The application is designed on MATLAB App Designer. MATLAB is a registered trademark of The MathWorks, Inc. This work was supported by funds to R.A.S. from the Natural Sciences and Engineering Research Council of Canada [grant number 2018-05021] and the Canada Foundation for Innovation [project number 11860].

# References 




