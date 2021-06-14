# containerized_metabolomics

A virtual machine equipped with (mainly R-based) tools and packages for metabolomics analysis.

## Background and rationale

Main focus is on an initial descriptive statistics module yielding p-values of differentail expression as well as Fold Change values; followed by pathway analyses.
  
Heavily builds on metaboanalyst and two example notebooks simply recreate demos at https://www.metaboanalyst.ca website.
  
The additional advantage of this container (besides easy deployment) is that it hosts R scripts wrapping some metaboanalyst pipelines together to provide a shortcut both working from a notebook or via command line.

## Interactive notebook

To run the example notebook, simply pull the Docker container with

```
docker pull  szabogtamas/containerized_metabolomics:latest
```

Then start the RStudio server with something like:

```
docker run --rm -p 127.0.0.1:8787:8787 -v $PWD:/home/rstudio/local_files -e USERID=$UID -e PASSWORD=<yourpassword> szabogtamas/containerized_metabolomics
```

Your work will only get saved if it is in the `local_files` directory of the container (attached to the folder from where you ran the above command).  
The recomended usage is to open a template notebook in `repo_files` (e.g.: `repo_files/notebooks/script_capability_demo.Rmd`) and save it in the `local_files` folder before starting to work on actual data.

## Nextflow pipeline

To compile a standardized report from a batch of samples, a Nextflow pipeline is also available on the image. If the measured data is in the `input_data` folder of the current directory (alternatively, using test data in `data_examples/pipeline_input_data`), a basic report can be generated issuing the following command:

```
docker run -it -v $PWD:/home/rstudio/local_files \
  szabogtamas/containerized_metabolomics \
  nextflow run /home/rstudio/repo_files/nextflow/main.nf \
  --input_folder /home/rstudio/local_files/input_data \
  --report_filename "my_test_metabolomics_report.pdf" \
  --report_title "A test report" \
  --report_author "Me"
```

Additional tables will be in `tables` folder and individual figures in `figures` folder.