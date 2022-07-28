Rho: Application to analyse electrical resistivity

Meryem Berrada1 [mberrada@uwo.ca](mailto:mberrada@uwo.ca) https://orcid.org/0000-0002-4934-0917

Richard A. Secco1 [secco@uwo.ca](mailto:secco@uwo.ca) https://orcid.org/0000-0001-5029-659X

1Department of Earth Sciences, University of Western Ontario, London, ON, Canada

Abstract

Electrical resistivity measurements of metals often use the four-wire method, where data acquisition alternates between temperature and voltage. The temperature and voltage values are similar and manual selection of each set of measurements is time consuming. *Rho* is intended to make electrical resistivity data analysis faster and more reliable by automatically selecting the temperature and voltage measurements using signal processing. The time required to process 52 data points (which consists of approximately 6000 measurements) using *Rho* is less than one minute compared to up to 4 hours for manual processing. *Rho* is stored in *FigShare* and can be used by any research group conducting alternate temperature and voltage measurements to obtain electrical resistivity at either ambient or extreme conditions.

Keywords

Electrical Resistivity; Thermal Conductivity; Data analysis; Voltage drop

Introduction

Electrical resistivity of materials is a basic property that is measured and used in a large variety of studies. In one type of study, electrical resistivity data are used to calculate thermal conductivity and interior heat flow of terrestrial-type bodies. The electrical resistivity of a sample is often obtained using the four-wire method as it reduces the contribution of the electrodes that is present in the two-wire method. The four-wire method is used in the cooling of industrial devices [1], soldering processes [2], electronic systems [3], temperature determination at depth [4], and simulations of physical properties of terrestrial-type bodies at high pressures and high temperatures. High pressure and high temperature measurements of electrical resistivity (ρ) in high pressure multi-anvil presses are used to calculate thermal conductivity (k) and the heat flow in the interiors of terrestrial-type bodies. Ezenwa and Secco [5] developed a four-wire method to solve the voltage bias, caused by passing a test current in a single direction, which previously resulted in erroneous analysis of ρ measurements due to parasitic voltages caused by the Seebeck effect or other sources. This polarity switching method allowed successful measurements of ρ of Zn [5], Nb [6], Ni [7], Co [8], Cu [9], Fe [10,11], Ag [12], W [13], Re [13], Au [14], Fe-Si alloys [15-17], and Fe-S alloys [18]. Three other papers from different research groups report using similar data collection practices, while several other papers use the four-wire method without a polarity switch [19]. In these measurements, the sample is located in the middle of a high-pressure cell, between two metal discs which are each contacted by a thermocouple. The metal discs, ensuring contact between the thermocouples and the sample, are typically composed of W. The thermocouples are typically Type-C thermocouples (one leg is made of 95%W/5%Re and the other leg is made of 74%W/26%Re) where voltage (i.e., thermoelectric emf) is measured and then converted to temperature after the experiment. The voltage drop across the sample for resistivity measurement is made between a pair of wires on opposite sides of the sample (one wire leg from each thermocouple), and therefore includes the voltage contribution of the discs. Thus, the ρ of the sample is obtained by subtracting the disc contribution from the total measurement.

The four-wire method incorporates a polarity switch and a mode switch, and results in a pattern of alternating values of temperature (emf) and voltage in a single column when acquired by a single voltmeter. A complete measurement produces a pattern similar to the following: temperature emf before voltage measurement (Tb), positive voltage drop with polarity in one direction (V+), negative voltage drop with polarity in opposite direction (V), temperature emf after voltage measurement (Ta). The order of V+ and V- may be alternated depending on the position of the polarity switch. At the University of Western Ontario, in the High-Pressure High-Temperature Laboratory, a programmable Keysight B2961 power supply is used to provide a constant direct test current of 0.2 A while data are acquired by a programmable Keysight 34470A meter operating at 20 Hz and 1 µV resolution. The application *Rho* calculates the average of 10 data points for each temperature section and all available data points for each voltage section. It is expected that temperature remains stable during these 10 data points. At high temperatures, when the frequency of measurements is increased and the quantity of data points in each section is low (due to the rapid switching between each section), *Rho* selects all the available data for each section. The values in each selection are averaged. Then the V+ and V- selections are averaged to give the voltage drop (V), while Tb and Ta selections are averaged to give the temperature (T). The voltage drop and current (I) are then used to calculate the electrical resistance (R) of the sample using Ohm’s law:

|   | (1) |
|---|-----|

The electrical resistivity is then calculated using Pouillet’s law:

|   | (2) |
|---|-----|

where *A* is the sample’s cross-sectional area and *l* is the sample’s length which are both measured on the post-experiment recovered sample. Values of ρ are then used to calculate the electronic component of thermal conductivity (k) via the empirical Wiedemann-Franz law:

|   | (3) |
|---|-----|

where T is temperature and Lo is the theoretical Sommerfeld value (Lo = 2.44∙10-8 W∙Ω∙K-2) of the Lorenz number. The error in temperature measurement corresponds to the standard deviation of the data points of each temperature section. The error in ρ is obtained by error propagation using the uncertainties in sample geometry and standard deviation of the voltage measurements. Similarly, the error in k corresponds to the propagation of the error in ρ. The selection of temperature and voltage values is typically done manually, and, to our knowledge, no software available in the literature is capable of automatically processing the temperature and voltage signals.

**Implementation and architecture**

The main page of *Rho* is displayed in **Fig.1a**. The application requires data to be imported in .xlsx, .csv, or .xls format, where the Tb, V+, V-, and Ta data alternate within the first column of the file (see **Fig.1b**). The first figure (see **Fig. 1a**) that is generated when the file is loaded is of the raw data. Noise, or abrupt fluctuation in data, that are not representative of the measurements can be removed by selecting the starting and ending indices (x-axis) of the section. The application attributes NaN values to each data point of this section, which removes them from the future analysis. When ‘Plot Data’ is selected, *Rho* identifies the temperature and voltage sections based on deviations from the main increasing trend, which corresponds to the temperature increase with time.

![](media/f99bac8b1aa700b0bd42097d502a47b2.png)

**Fig.1:** **(a)** Main screen. The input parameters options are displayed on the left panel. A visualization of the Data selection, resulting resistivity and output results are displayed on the right. An example raw data set of Fe at 8 GPa is displayed in Figure 1 [16]. **(b)** Example of raw data Excel file for Fe at 8 GPa [16] where all measurements are within the first column of the file. The identification of parameter sequence is shown in the right column for the first two sets of data for clarity here only (this not included in the original data file).

First, the intervals of negative values are labelled as voltage measurements (Vtotal). This temporary selection is displayed in Figure 1 of the application (see **Fig. 2a** below). In the first scenario, the pattern Tb V+ V- Ta is considered and Tb and Ta are temporarily identified. To satisfy this scenario, the difference between Tb and Vtotal must be smaller than the difference between Ta and Vtotal considering that the part of Vtotal closer to Tb would be V+ and that closer to Ta would be V- (and V- \< V+). *Rho* uses a loop to evaluate this scenario at each deviation from the main increasing trend. Whenever the scenario is observed, V+ is temporarily defined as the first trend in Vtotal. In the second scenario, the pattern Tb V- V+ Ta is considered, and a similar loop is applied to identify the intervals where the difference between Tb and Vtotal is larger than that of Ta and Vtotal. Here, the part of Vtotal closer to Tb would be V- and that closer to Ta would be V+. Another analysis of the first scenario identifies the stable intervals on the main increasing trend right before V+ as Tb, and that right after V- as Ta. In other words, fluctuations, or noise, between temperature and voltage measurements are automatically ignored. A similar analysis of the second scenario identifies the stable intervals on the main increasing trend right before V- as Tb and that right after V+ as Ta. The corresponding indexes (x-axis) and data (y-axis) for each selection are combined into variables for Tb, V-, V+, and Ta. The ‘Outliers Degrees of Freedom’ parameter evaluates the differences between the trends of V- and V+ and identifies data points outside of these trends as outliers. A larger degree of freedom will include more data points in the final variables for V- and V+, while a smaller degree of freedom will ignore more data points. The value 20 seems to be the best fit for this parameter, although it can be adjusted by the user. In other words, fluctuations of 20 times and more are defined as outliers. The final selection is displayed in Figure 2 of the application (see **Fig. 2a** and **2d**). The contribution of the discs to the measured ρ is subtracted by fitting the ρ of the selected material. The ‘Discs’ dropdown menu contains options of W from 2-5 GPa [13], Re from 2-5 GPa [13] and Pt at 1 atm, 10 and 20 GPa [20]. Supplementary options may be added to future updates of the application as data become available in the literature for different metals used for the discs. Temperature is automatically converted from emf to (K) using the selected thermocouple calibration equation. The ‘Thermocouple’ dropdown menu contains options for Type-C [21], Type-S at 1 atm [22] and Type-S from 1-5 GPa thermocouples [23]. The input parameters are then used in equation (1), (2) and (3) to plot ρ in Figure 3 of the application and output values of T (K), ρ (µΩ·cm), error in ρ (µΩ·cm), k (Wm-1K-1) and error in k (Wm-1K-1) (see **Fig. 2b** and **2c**). The values are displayed in the ‘Output’ tab and can be copied and pasted on another platform for further analysis.

![](media/3537476049f5f06fc84d7315f0cad712.png)

**Fig.2:** **(a)** Temporary and final selections of the temperature and voltage sections. **(b)** Calculated resistivity compared with that of Fe at 1 atm [24]. **(c)** Output values of T (K), ρ (µΩ·cm), error in ρ (µΩ·cm), k (Wm-1K-1) and error in k (Wm-1K-1). **(d)** Example of the final selections of one measurement for each parameter.

Examples of running the software are displayed in **Fig.2**. After loading the raw data file, Figure 1 of the application should automatically display the raw data. If the raw data file is in the incorrect format, an error sound will be set, and Figure 1 will not display the raw data. The users should confirm the temperature and voltage selections since similarity in these values might result in the incorrect identification of some data points. The computed ρ can be compared with that of Fe at 1 atm [24].

**Availability**

Users may choose to download *Rho.exe* file from *FigShare* or *Rho.m* from *GitHub*.

*Rho.exe* – Any system capable of running a .exe file.

*Rho.m* – Any system capable of running Matlab R2012a or higher.

Tested in Matlab R2012a, likely to work in earlier versions also.

**Name:** FigShare Archive

**Persistent identifier:** 10.6084/m9.figshare.19175432

**Licence:** CC by 4.0

**Publisher:** Meryem Berrada

**Version published:** 1.0

**Date published:** 15/02/22

**Code repository** (e.g. SourceForge, GitHub etc.) (required)

**Name:** GitHub

**Identifier:** <https://github.com/meryemberradauwo/Rho.git>

**Licence:** GNU General Public License v3.0

**Date published:** 02/02/22

**Statement of need**

The application *Rho* minimizes the data analysis time and uses signal processing to identify the temperature (emf) and voltage drop measurements from a single column of alternating measurements resulting from the four-wire method. The confirmation of the selection by the user remains required considering the limitations of the attributed thresholds for outliers. This software (*Rho.exe*) can be used to analyse electrical resistivity for applications to the cooling of industrial devices, soldering processes, electronic systems, temperature determination at depth, and simulations of physical properties of terrestrial-type bodies at high pressures and high temperatures. *Rho* can be used by any research group conducting alternate temperature (emf) and voltage measurements to obtain electrical resistivity at either ambient or extreme conditions. Modifications to the software, such as changing the identification thresholds, adding TC types or disc materials, can be done through *Rho.m* via the MATLAB platform. The authors are continuously working to improve the software and implement more features.

We encourage the reader to download Rho via MATLAB [25] or through this article. Inquiries can be addressed to the corresponding author of this article, or through the contact information displayed in the software.

**Acknowledgements**

The application is designed on MATLAB App Designer. MATLAB is a registered trademark of The MathWorks, Inc.

**Funding statement**

This work was supported by funds to R.A.S. from the Natural Sciences and Engineering Research Council of Canada [grant number 2018-05021] and the Canada Foundation for Innovation [project number 11860].

**Competing interests**

The authors declare that they have no competing interests.

**References**

1.  Hapenciuc, C. L., Negut, I., Borca-Tasciuc, T., & Mihailescu, I. N. (2019). A steady-state hot-wire method for thermal conductivity measurements of fluids. International Journal of Heat and Mass Transfer, 134, 993-1002. <https://doi.org/10.1016/j.ijheatmasstransfer.2019.01.098>
2.  Tillmann, W., Sievers, N., Henning, T., & Jakimenko, D. (2017). FEM study of analyzing the electrical resistance of brazed joint by the 4-wire technique for quality assurance. Measurement, 104, 43-49. <https://doi.org/10.1016/j.measurement.2017.03.015>
3.  Krishnamurthy, R., Bharatiraja, C., Adedayo, Y., Tariq, M., & Azeem, A. (2019). Locating Wire Fault in Controller Area Network Based on Kelvin (Four-Wire) Resistance Approach. Paper presented at the Applications of Computing, Automation and Wireless Systems in Electrical Engineering, Singapore. <https://doi.org/10.1007/978-981-13-6772-4_98>
4.  Erkan, K., Akkoyunlu, B., Balkan, E., & Tayanç, M. (2017). A portable borehole temperature logging system using the four-wire resistance method. Journal of Geophysics and Engineering, 14(6), 1413-1419. <https://doi.org/10.1088/1742-2140/aa7ffe>
5.  Ezenwa, I. C., & Secco, R. A. (2017a). Constant electrical resistivity of Zn along the melting boundary up to 5 GPa. High Pressure Research, 37(3), 319-333. <https://doi.org/10.1080/08957959.2017.1340473>
6.  Ezenwa, I. C., & Secco, R. A. (2017b). Electronic transition in solid Nb at high pressure and temperature. Journal of Applied Physics, 121(22), 225903. <https://doi.org/10.1063/1.4985548>
7.  Silber, R. E., Secco, R. A., & Yong, W. (2017). Constant electrical resistivity of Ni along the melting boundary up to 9 GPa. Journal of Geophysical Research: Solid Earth, 122(7), 5064-5081. <https://doi.org/10.1002/2017JB014259>
8.  Ezenwa, I. C., & Secco, R. A. (2017c). Invariant electrical resistivity of Co along the melting boundary. Earth and Planetary Science Letters, 474, 120-127. <https://doi.org/10.1016/j.epsl.2017.06.032>
9.  Ezenwa, I. C., Secco, R. A., Yong, W., Pozzo, M., & Alfè, D. (2017). Electrical resistivity of solid and liquid Cu up to 5 GPa: Decrease along the melting boundary. Journal of Physics and Chemistry of Solids, 110, 386-393. <https://doi.org/10.1016/j.jpcs.2017.06.030>
10. Silber, R. E., Secco, R. A., Yong, W., & Littleton, J. A. H. (2018). Electrical resistivity of liquid Fe to 12 GPa: Implications for heat flow in cores of terrestrial bodies. Scientific Reports, 8(1), 10758. <https://doi.org/10.1038/s41598-018-28921-w>
11. Yong, W., Secco, R. A., Littleton, J. A. H., & Silber, R. E. (2019). The Iron Invariance: Implications for Thermal Convection in Earth's Core. Geophysical Research Letters, 46(20), 11065-11070. <https://doi.org/10.1029/2019GL084485>
12. Littleton, J. A. H., Secco, R. A., & Yong, W. (2018). Decreasing electrical resistivity of silver along the melting boundary up to 5 GPa. High Pressure Research, 38(2), 99-106. <https://doi.org/10.1080/08957959.2018.1435786>
13. Littleton, J. A. H., Secco, R. A., Yong, W., & Berrada, M. (2019). Electrical resistivity and thermal conductivity of W and Re up to 5 GPa and 2300 K. Journal of Applied Physics, 125(13), 135901. <https://doi.org/10.1080/08957959.2018.1435786>
14. Berrada, M., Secco, R. A., & Yong, W. (2018). Decreasing electrical resistivity of gold along the melting boundary up to 5 GPa. High Pressure Research, 38(4), 367-376. <https://doi.org/10.1080/08957959.2018.1493476>
15. Silber, R. E., Secco, R. A., Yong, W., & Littleton, J. A. H. (2019). Heat Flow in Earth's Core From Invariant Electrical Resistivity of Fe-Si on the Melting Boundary to 9 GPa: Do Light Elements Matter? Journal of Geophysical Research: Solid Earth, 124(6), 5521-5543. <https://doi.org/10.1029/2019JB017375>
16. Berrada, M., Secco, R. A., & Yong, W. (2020). Electrical resistivity measurements of Fe-Si with implications for the early lunar dynamo. Journal of Geophysical Research: Planets. <https://doi.org/10.1029/2020JE006380>
17. Berrada, M., Secco, R. A., & Yong, W. (2021). Adiabatic Heat Flow in Mercury’s Core from Electrical Resistivity Measurements of Liquid Fe-8.5wt%Si to 24 GPa. Earth and Planetary Science Letters. <https://doi.org/10.1016/j.epsl.2021.117053>
18. Littleton, J. A. H., Secco, R. A., & Yong, W. (2021). Electrical Resistivity of FeS at High Pressures and Temperatures: Implications of Thermal Transport in the Core of Ganymede. Journal of Geophysical Research: Planets, 126(5). <https://doi.org/10.1029/2020JE006793>
19. Berrada, M., & Secco, R. A. (2021). Review of Electrical Resistivity Measurements and Calculations of Fe and Fe-Alloys Relating to Planetary Cores. Frontiers in Earth Science, 9(802). Review. <https://doi.org/10.3389/feart.2021.732289>
20. Gomi, H., & Yoshino, T. (2019). Resistivity, Seebeck coefficient, and thermal conductivity of platinum at high pressure and temperature. Physical Review B, 100(21), 214302. <https://link.aps.org/doi/10.1103/PhysRevB.100.214302>
21. Omega. (2019a). Revised Thermocouple Reference Tables - Type C. In T. T. C. C) (Ed.): Reference Tables N.I.S.T. Monograph 175 Revised to ITS-90.
22. Omega. (2019b). Revised Thermocouple Reference Tables - Type S. In T. T. C. C) (Ed.): Reference Tables N.I.S.T. Monograph 175 Revised to ITS-90.
23. Getting, I. C., & Kennedy, G. C. (1970). Effect of Pressure on the emf of Chromel‐Alumel and Platinum‐Platinum 10% Rhodium Thermocouples. Journal of Applied Physics, 41(11), 4552-4562. <https://doi.org/10.1063/1.1658495>
24. Chu, T. K., & Chi, T. C. (1981). Properties of Selected Ferrous Alloying Elements (Vol. III). Washington: McGraw-Hill.
25. Berrada, M. (2022). Rho (https://www.mathworks.com/matlabcentral/fileexchange/\<...\>), MATLAB Central File Exchange.