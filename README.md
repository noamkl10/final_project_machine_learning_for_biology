# final_project_machine_learning_for_biology

for running the code, first download the data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106487 and open it.
put the folder in a folder name data, and then run CreatingAndSortingMatrix.py with the line: python CreatingAndSortingMatrix.py

after that go to R, and install celcall:
library(devtools)
devtools::install_github("ShellyCoder/cellcall")

watch out that there is a lot of depndencies to install so the R code will run.

after that run AnalyzeTheData.R and then you will get in the all the pictures. 

both files, espically the second take long time to run.
