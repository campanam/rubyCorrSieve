# rubyCorrSieve Change Log

Michael G. Campana, 2019  
Smithsonian Conservation Biology Institute  
Contact: campanam@si.edu  

Public domain rubyCorrSieve updates (versions >=1.7-0) are US government works  
Unmodified *CorrSieve* source code (versions <= 1.6-5) copyright (c) Michael G. Campana, 2010-2011  

### Version 1.7-0  
Removal of the method 'each_permutation' from the Array class.  
Improvement of the 'pvalue' method.  
Conversion of the main processing section to the method 'process_data'.  
Revisions to the splash screen for updated licensing information.  
Addition of method 'splash_screen' to control revised program.  
Addition of methods 'show_license' and 'gnu_license' to detail licensing information.  

### Version 1.6-5  
Added compatibility with INSTRUCT output.  
Added *F*<sub>ST</sub> statistic optimisation abilities to Ruby version.  
Improved algorithm for *F*<sub>ST</sub> optimisation using Q-matrix correlations.  

### Version 1.6-4  
R version improvements. No Ruby version changes.  

### Version 1.6-3
R version improvements and documentation corrections. No Ruby version changes.  

### Version 1.6-2
R version improvements. No Ruby version changes.  

### Version 1.6-1
Now reads the *K* values and the number of runs per *K* from the raw data.  
Program no longer requires the same number of runs per *K*.  
Program no longer requires the minimum *K* to be 1.  
Program no longer reads correlation data when the correlation option is turned off.  
Fixed bug that did not correctly validate the ln P(D) yes/no entry.  

### Version 1.51  
Minor improvement of yes/no interface algorithm.  
Reworded Ln P(D) and Δ*K* prompt.  
Better explanation of installing and executing Ruby.  

### Version 1.5  
Added some licensing information to beginning screen.  
Added the ability to omit correlation matrices calculation.  
Added the ability to use the average maximum correlation method rather than the rows-and- columns method.  
Changed wording of *F*<sub>ST</sub> prompt.  

### Version 1.4  
Added the ability to calculate Δ*F*<sub>ST</sub>.  
Changed output files from list form to spreadsheet form.  
Corrected a bug that caused *F*<sub>ST</sub> values for *K*'s greater than 9 to be imported incorrectly.  
Corrected a bug that caused cluster assignment data to be imported incorrectly if population information is included.  

### Version 1.32
Added the standard deviation within *F*<sub>ST</sub>s within clusters to mean *F*<sub>ST</sub>s output file.  
Added the overall standard deviation of *F*<sub>ST</sub>s ignoring clusters to mean *F*<sub>ST</sub>s output file.  
Added the mean of the within-cluster standard deviation of *F*<sub>ST</sub>s to mean *F*<sub>ST</sub>s output file.  
Added the standard deviation of the standard deviation within *F*<sub>ST</sub>s within clusters to mean *F*<sub>ST</sub>s output file.  

### Version 1.31  
Generalised function to calculate standard deviation.  
Added the standard deviation of the mean *F*<sub>ST</sub>s to mean *F*<sub>ST</sub>s output file.  

### Version 1.3
Added the ability to calculate mean *F*<sub>ST</sub>s.  

### Version 1.2
Added the ability to calculate Δ*K* (Evanno et al. 2005).  Streamlined the data import algorithm.  

### Version 1.11  
Fixed a bug in the output file which plotted the top row as *K* rather than the number of *K*.  

### Version 1.1  
Added the ability to calculate an estimated *p* in addition to an exact *p*.

### Version 1.02  
Added the ability to change the output files’ save path.  

### Version 1.01  
Streamlined some of the algorithms.  

### Version 1.0  
First functioning version (*Q*-matrix correlations only).  
