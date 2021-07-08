import sys
sys.path.append("/home/william1/bin/shahlabdata")
import shahlabdata.dlp
import pandas as pd
import os
os.environ['ISABL_API_URL']='https://isabl.shahlab.mskcc.org/api/v1/'

individual = "OV2295"
analysis_results = shahlabdata.isabl.get_individual_results(individual)

analysis_results.to_csv("isablpaths.csv")
