## schnapps

schnapps needs to be run at the patient level (ie isabl_patient_id), the required inputs are:

* cn binned data from HMMcopy (result_type == 'reads')
* metrics files from SCDNA-ANNOTATION
* COUNT-HAPS files

In this example I've pulled info (for sample OV2295) from isabl (get_paths_isabl.py) and then manually pulled out the relevent paths in isablpaths.csv and added them to the script.

The Dockerfile includes all relevent software. It can be pulled using
``` 
singularity pull docker://marcjwilliams1/schnapps
```

There should be 4 outputs from this app

* Rdata file with schnapps object
* *csv.gz file with allele specific copy number calls
* heatmap.png copy number heatmap
* qc.png figure with various QC plots

`run-schnapps.R` is an R script that runs schnapps, `jobscript.sh` is an example LSF job submission script using singularity.