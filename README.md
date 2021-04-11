# containerized_metabolomics

A virtual machine equipped with (mainly R-based) tools and packages for metabolomics analysis.

## Background and rationale

Main focus is on an initial descriptive statistics module yielding p-values of differentail expression as well as Fold Change values; followed by pathway analyses.
  
Heavily builds on metaboanalyst and two example notebooks simply recreate demos at https://www.metaboanalyst.ca website.
  
The additional advantage of this container (besides easy deployment) is that it hosts R scripts wrapping some metaboanalyst pipelines together to provide a shortcut both working from a notebook or via command line.

## Usage

To run the example notebook, simply pull the Docker container with

```
docker pull  szabogtamas/containerized_metabolomics:latest
```

Then start the RStudio server:

```
docker run --rm -p 127.0.0.1:7787:8787 -v $PWD:/home/rstudio/local_files -e USERID=$UID -e DISABLE_AUTH=true szabogtamas/containerized_metabolomics
```

Your work will only get saved if it is in the `local_files` directory of the container (attached to the folder from where you ran the above command).  
The recomended usage is to open a template notebook in `repo_files` (e.g.: `repo_files/notebooks/script_capability_demo.Rmd`) and save it in the `local_files` folder before starting to work on actual data.