Exploratory Data Analysis

This folder contains some unorganized scripts to explore the data in the `DepLabData` and find proteins that change across conditions.

Description of each file:

# eda.R

Given the profiles of all the samples, it computes the number of peaks and the location of the highest peak for each sample, as well as the average location and the average number of peaks for each condition. It also produces the heatmaps of these quantities (that should be turned into interactive heatmaps in the app).

# changing_proteins.Rmd

This is where we define the method to rank the proteins in the order of how much they change across condition. Briefly, we look for proteins that have profiles that are highly correlated between replicates in the same condition and poorly correlated across conditions.

# complexes.Rmd

This was an attempt at identifying new protein complexes. Currently not in a publishing form.

The files `interacting_complexes.txt` and `non-interacting_proteins.txt` are useful for this task.


