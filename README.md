# Ovarian-stromal-damage-ST-analysis

See the spec-list.txt and r-spec-list.txt for packages installable with conda.

To recreate the environment on another system, use the following command:
conda create --name <envname> --file spec-list.txt

After creating the environment with specified R and python versions, you can choose to download the specific R packages specified in r_packages.csv.

The preprocessing files are outside of individual folders for result generation for each figure, and should be run before running analysis. The folder other_utils contain handy functions, such as cropping scanpy objects focusing on a specified cluster and so on.