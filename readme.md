## Repository for the multieditR package 

multiEditR helps identify edits from Sanger Sequences

To install, use:

```
devtools::install_github("jeremymchacon/multieditR")
```

Basic functionality:

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

