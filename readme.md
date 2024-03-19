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
devtools::install_github("MoriarityLab/multiEditR.pckg")
```

## Basic functionality:

Here, we load in two example sanger sequences, then detect whether "A" bases were edited within the motif "AGTAGCTGGGATTACAGATG" using detect_edits(), the main function of the package. 
```
library(multiEditR.pckg)
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
create_multiEditR_report(fits, params, "my_html_report.html")
```

Known issues:

- sometimes motifs get trimmed due to quality, resulting in chromatograms which seem shorter than expected

- when the control sequence is revcom, printing the chromatogram might not be working correctly

- sometimes the relative position of the edited bases are off by a fixed amount. I think this occurs when indels are removed but not accounted for. 

Please add any other feature requests or issues you find to the issues page. Thank you!
