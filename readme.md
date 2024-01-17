## Repository for the multieditR package 

multiEditR helps identify edits from Sanger Sequences. The algorithms were developed by Mitch Kleusner, and put into this R package by Jeremy Chacón. 

To install, use:

```
devtools::install_github("jeremymchacon/multieditR")
```

## Basic functionality:

Here, we load in two example sanger sequences, then detect whether "A" bases were edited within the motif "AGTAGCTGGGATTACAGATG" using detect_edits(), the main function of the package. 
```
library(multieditR)
sample_file = system.file("extdata", "RP272_cdna_wt.ab1", package="multieditR")
ctrl_file = system.file("extdata", "RP272_cdna_ko.ab1", package="multieditR")
motif = "AGTAGCTGGGATTACAGATG"
wt = "A"
edit = "G"

fit = detect_edits(
    # Set parameters
  sample_file = sample_file,
  ctrl_file = ctrl_file,
  p_value = 0.0001, 
  phred_cutoff = 0.0001,
  motif = motif, 
  motif_fwd = TRUE,
  use_ctrl_seq = FALSE,
  wt = wt, 
  edit = edit 
)
```

After fitting is complete, multiple graphs and tables can be made:

```

# chromatogram of predicted edits in sample
geom_chromatogram(fit$sample_sanger,fit$sample_locs[1], fit$sample_locs[2])

# chromatogram of predicts "edits" in control
geom_chromatogram(fit$ctrl_sanger,fit$ctrl_locs[1], fit$ctrl_locs[2])

# percent signal of basecalls in sample
plot_raw_sample(fit)

# percent noise of basecalls 
plot_trimmed_sample(fit)

# height of significant and non-significant edits
plot_editing_barplot(fit)

# main table
result = fit$sample_data 
# secondary table
stats = fit$statistical_parameters 
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
print(fit1$sample_name)
geom_chromatogram(fit1$sample_sanger,fit1$sample_locs[1], fit1$sample_locs[2])
```

Or make a report containing results from all of the samples:

```
create_multieditR_report(fits, params, "my_html_report.html")
```