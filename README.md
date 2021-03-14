# containerized_metabolomics

A virtual machine equipped with (mainly R-based) tools and packages for metabolomics analysis.
  
Main focus is on an initial descriptive statistics module yielding p-values of differentail expression as well as Fold Change values; followed by pathway analyses.
  
Heavily builds on metaboanalyst and two example notebooks simply recreate demos at https://www.metaboanalyst.ca website.
  
The additional advantage of this container (besides easy deployment) is that it hosts R scripts wrapping some metaboanalyst pipelines together to provide a shortcut bot working from a notebook or via command line.