## Repository for the multiEditR package 

multiEditR detects and quantifies base edits from Sanger sequencing. Algorithms were developed by Mitch Kluesner, and implemented by Jeremy Chacón. 

multiEditR inputs a sample Sanger sequence, and either a control Sanger or fasta file, a motif within which to look for edits, the base being edited and expected edited base, and return whether edits were found.

Furthermore, if the control sequence is rev-com from the sample sequence, this package detects that and rev-coms the control sequence prior to testing.



To install, first you must have devtools:

```
install.packages("devtools")
```

Then, use devtools to install the multiEditR package:

```
devtools::install_github("jeremymchacon/multiEditR.pckg")

# Note: above is a temporary home. The permanent home will be:
# devtools::install_github("MoriarityLab/multiEditR.pckg")

```

## Basic functionality:

Here, we load in two example sanger sequences, then detect whether "A" bases were edited within the motif "AGTAGCTGGGATTACAGATG" using detect_edits(), the main function of the package. 
```
library(multiEditR)
sample_file = system.file("extdata", "RP272_cdna_wt.ab1", package="multiEditR")
ctrl_file = system.file("extdata", "RP272_cdna_ko.ab1", package="multiEditR")
motif = "AGTAGCTGGGATTACAGATG"
wt = "A"
edit = "G"

fit = detect_edits(
  sample_file = sample_file,
  ctrl_file = ctrl_file,
  p_value = 0.0001, 
  phred_cutoff = 0.0001,
  motif = motif, 
  wt = wt, 
  edit = edit 
)
```

After fitting is complete, multiple graphs and tables can be made:

```

# obtain full data table 
tbl = results(fit)
# optionally save it
writexl::write_xlsx(tbl, "my_results.xlsx") 


# chromatogram of predicted edits in sample, for all motifs found
plot_sample_chromatogram(fit)

# chromatogram of predicts "edits" in control, for all motifs found
plot_control_chromatogram(fit)

# percent signal of basecalls in sample
plot_raw_sample(fit)

# height of significant and non-significant edits
plot_editing_barplot(fit)

# list of all internals
str(fit)
```

## Batch Mode

The package can also take in a parameters spreadsheet and run the edit-detection algorithm on many samples at once, then generate a single html report from all of the samples. 

The package comes with an example which can be loaded like so:

```
params = load_example_params()
head(params)
```

Note that you can save a skeleton parameters spreadsheet to modify with your own data:

```
save_batch_skeleton("./my_parameters.xlsx")
```

We use "detect_edits_batch" to fit all samples simulateously:

```
fits = detect_edits_batch(params)
```

We can access summary tables for all samples:

```
data.tbl = get_batch_results_table(fits)
stats.tbl = get_batch_stats_table(fits)
```

Or access a single fit model

```
fit1 = fits[[1]] 
results(fit1)
plot_sample_chromatogram(fit1)
```

Or make a report containing results from all of the samples:
```
# First set your directory or interest:
setwd("your/directory/of/interest/to/write/files")

# Run the report, which may take a few minutes
create_multiEditR_report(fits, params, "my_html_report.html")
```

Please add any other feature requests or issues you find to the issues page. Thank you!
