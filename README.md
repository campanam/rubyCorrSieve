# rubyCorrSieve  
Michael G. Campana  
Smithsonian Conservation Biology Institute  
Contact: campanam@si.edu  

Ruby implementation of *CorrSieve*  


## Licensing  
Original Ruby source code (*CorrSieve* versions <= 1.6-5) copyright (c) Michael G. Campana, 2010-2011 is licensed under the GNU General Public License (version 3 or later). See included LICENSE file for details.  

Public domain updates by Michael G. Campana (2019) to the original Ruby source code (*CorrSieve* versions >= 1.7-0) are United States government works. These modifications are annotated in the modified source code.  

## Introduction
*CorrSieve* is a Ruby and R package that filters *Q* value output from the programs STRUCTURE (Pritchard et al. 2000) and INSTRUCT (Gao et al. 2007) by correlation values. It outputs matrices showing significant correlations between individual runs for each *K*. It can also calculate Δ*K* (Evanno et al. 2005), mean *F*<sub>ST</sub>s and Δ*F*<sub>ST</sub>. These measures help identify meaningful values of *K*.  

## Installation and Compatibility  
rubyCorrSieve is compatible with Windows, Linux, and UNIX (including macOS) operating systems. rubyCorrSieve requires the Ruby interpreter. Installation files are available at www.ruby-lang.org/en/downloads. Install the appropriate interpreter for your operating system.  

Clone this repository to your system. Using a Linux/UNIX command line this can be performed using `git`.:  
`git clone https://github.com/campanam/rubyCorrSieve`  

You may need to make the CorrSieve-1.7-0.rb file executable. On Linux/UNIX, enter:
`chmod +x rubyCorrSieve/*.rb`  

Move the CorrSieve-1.7-0.rb executable and LICENSE file to your chosen execution directory.  

## Bugs and Contributing
Please report all bugs (and any suggestions for improvements) to Michael G. Campana (campanam@si.edu).  

## CorrSieve Citation  
Campana, M.G. et al. 2011. *CorrSieve*: software for summarizing and evaluating Structure output. Mol. Ecol. Resour. 11:349-352. https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1755-0998.2010.02917.x  

## References
Evanno et al. 2005. Detecting the number of clusters of individuals using the software STRUCTURE: a simulation study. Mol. Ecol. 14: 2611-2620. https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1365-294X.2005.02553.x.  

Gao et al. 2007. A Markov Chain Monte Carlo approach for joint inference of population structure and inbreeding rates from multilocus genotype data. Genetics. 176: 1635-1651. https://www.genetics.org/content/176/3/1635.  
