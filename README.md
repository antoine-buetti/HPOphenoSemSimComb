# HPOphenoSemSimComb
Calculates score combination of genomics scores (e.g., numerical ACMG) with HPO semantic similarity (e.g., Resnik distance between patient's phenotype and Human Phenotype Ontology (HPO) terms).

The procedure can be run directly in GitHub through Actions functionalities: just modify the input files in "Input/" and "output-$date/" will be created to be accessed in GitHub. 

Alternatively one can execute locally the different steps as executed by the GitHub Actions workflow, as:

- setup.sh : installs some dependencies
- download_HPO_data.sh :  downloads structured data from HPO ("Human Phenotype Ontology")
- build.sh : compiles and build the executables
- test.sh : performs unit tests
- run.sh : runs the program using the user-provided input found in "Input/", which results in an annotated table found in the directory "output-$date/"

twincacca

