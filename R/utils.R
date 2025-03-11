# utils.R

###########################################################################################
# Written by Jeremy Chacon and Mitchell Kluesner
#  
# This file is part of multiEditR (Multiple Edit Deconvolution by Inference of Traces in R)
# 
# Please only copy and/or distribute this script with proper citation of 
# multiEditR publication
###########################################################################################

bases = c("A", "C", "G", "T")
ACGT = bases

phred_scores = data.frame(stringsAsFactors=FALSE,
                          phred = c("!", "“", "$", "%", "&", "‘", "(", ")", "*", "+", ",", "–",
                                    ".", "/", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                                    ":", ";", "<", "=", ">", "?", "@", "A", "B", "C", "D", "E", "F",
                                    "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S",
                                    "T", "U", "V", "W", "X", "Y", "Z", "[", "\\", "]", "^", "_",
                                    "`", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l",
                                    "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y",
                                    "z", "{", "|", "}", "~"),
                          prob = c(1, 0.794328235, 0.501187234, 0.398107171, 0.316227766,
                                   0.251188643, 0.199526232, 0.158489319, 0.125892541, 0.1,
                                   0.079432824, 0.063095734, 0.050118723, 0.039810717, 0.031622777,
                                   0.025118864, 0.019952623, 0.015848932, 0.012589254, 0.01,
                                   0.007943282, 0.006309573, 0.005011872, 0.003981072, 0.003162278,
                                   0.002511886, 0.001995262, 0.001584893, 0.001258925, 0.001, 0.000794328,
                                   0.000630957, 0.000501187, 0.000398107, 0.000316228, 0.000251189,
                                   0.000199526, 0.000158489, 0.000125893, 1e-04, 7.94328e-05,
                                   6.30957e-05, 5.01187e-05, 3.98107e-05, 3.16228e-05, 2.51189e-05,
                                   1.99526e-05, 1.58489e-05, 1.25893e-05, 1e-05, 7.9433e-06,
                                   6.3096e-06, 5.0119e-06, 3.9811e-06, 3.1623e-06, 2.5119e-06, 1.9953e-06,
                                   1.5849e-06, 1.2589e-06, 1e-06, 7.943e-07, 6.31e-07, 5.012e-07,
                                   3.981e-07, 3.162e-07, 2.512e-07, 1.995e-07, 1.585e-07, 1.259e-07,
                                   1e-07, 7.94e-08, 6.31e-08, 5.01e-08, 3.98e-08, 3.16e-08,
                                   2.51e-08, 2e-08, 1.58e-08, 1.26e-08, 1e-08, 7.9e-09, 6.3e-09, 5e-09,
                                   4e-09, 3.2e-09, 2.5e-09, 2e-09, 1.6e-09, 1.3e-09, 1e-09, 8e-10,
                                   6e-10, 5e-10)
)

######
### Analysis Functions
######


# Check that a string is IUPAC notation
checkIUPAC = function(x){all(strsplit(x, "")[[1]] %in% c("A", "C", "G", "T", "U", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"))}

# convert nucleotide to a factor
nucleotide_factor = function(x){factor(x, levels = c("A", "C", "G", "T"))}

# get_trim_points returns a vector of length 2, the first point
# being the first position to be included, and the second being the second 
# position to be included. It is modified from abif_to_fastq.
# If the trimming would leave less than min_seq_len, it returns -1
# similarly, if the sequence is already shorter than min_seq_len, it returns -1
#' @importFrom sangerseqR read.abif
get_trim_points = function (path, cutoff = 0.001,
                            min_seq_len = 20, offset = 33)
{
  abif <- sangerseqR::read.abif(path)
  if (is.null(abif@data$PCON.2)) {
    message(sprintf("failed on %s", path))
    return()
  }
  nucseq <- substring(abif@data$PBAS.2, 1, length(abif@data$PLOC.2))
  if (!typeof(abif@data$PCON.2) == "integer") {
    num_quals <- utf8ToInt(abif@data$PCON.2)[1:length(abif@data$PLOC.2)]
  }else {
    num_quals <- abif@data$PCON.2[1:length(abif@data$PLOC.2)]
  }

  trim_msg <- "Sequence %s can not be trimmed because it is shorter than the trim\n               segment size"
  if (nchar(nucseq) <= min_seq_len) {
    warning(sprintf(trim_msg, path))
    return(-1)
  }
  scores = cutoff - 10^(num_quals/-10)
  running_sum <- rep(0, length(scores) + 1)
  for (i in 1:length(scores)) {
    num <- scores[i] + running_sum[i]
    running_sum[i + 1] <- ifelse(num < 0, 0, num)
  }
  trim_start <- min(which(running_sum > 0)) - 1
  trim_finish <- which.max(running_sum) - 2
  if (trim_finish - trim_start < min_seq_len - 1) {
    warning(sprintf(trim_msg, path))
    return(-1)
  }
  return(c(trim_start, trim_finish))
}

subchar <- function(string, pos, char) {
  for(i in pos) {
    substr(string, i, i) = char
  }
  string
}

#' @importFrom gamlss gamlss
#' @importFrom dplyr filter
#' @importFrom gamlss.dist qZAGA
#' @importFrom magrittr `%>%`
make_ZAGA_df = function(sanger_df, p_adjust){
  nvals <- list()
  nvals$A = dplyr::filter(sanger_df, max_base != "A")$A_area
  nvals$C = dplyr::filter(sanger_df, max_base != "C")$C_area
  nvals$G = dplyr::filter(sanger_df, max_base != "G")$G_area
  nvals$T = dplyr::filter(sanger_df, max_base != "T")$T_area

  n_models =   n_models <-lapply(nvals, FUN = function(x){
    set.seed(1)
    if((unique(x)[1] == 0 & length(unique(x)) == 1) |
       (unique(x)[1] == 0 & length(unique(x)) == 2 & table(x)[2] == 1))
    {x = c(rep(0, 989), 0.00998720389310502, 0.00998813447664401,0.009992887520785,
           0.00999585366068316, 0.00999623914632598, 0.00999799013526835, 0.010001499423723,
           0.0100030237039207, 0.0100045782875701, 0.0100048452355807, 0.0100049548867042)
    message("Replacement vector used for low noise.")
    } # add noise if all 0s, or all 0s and one other value.
    tryCatch(gamlss::gamlss((x)~1, family = gamlss.dist::ZAGA), error=function(e) # Progressively step up the mu.start if it fails
      tryCatch(gamlss::gamlss((x)~1, family = gamlss.dist::ZAGA, mu.start = 1), error=function(e)
        tryCatch(gamlss::gamlss((x)~1, family = gamlss.dist::ZAGA, mu.start = 2), error=function(e)
          tryCatch(gamlss::gamlss((x)~1, family = gamlss.dist::ZAGA, mu.start = 3), error=function(e) # additional step added.
            tryCatch(gamlss::gamlss((x)~1, family = gamlss.dist::ZAGA, mu.start = 4), error=function(e) # additional step added.
              gamlss::gamlss((x)~1, family = gamlss.dist::ZAGA, mu.start = mean(x))
            )
          )
        )
      )
    )
    # throws errors when a completely 0 vector
  })

  null_m_params = lapply(n_models, FUN = function(x){
    mu <- exp(x$mu.coefficients[[1]])
    sigma <- exp(x$sigma.coefficients[[1]])
    nu.logit <- x$nu.coefficients[[1]]
    nu <- exp(nu.logit)/(1+exp(nu.logit))
    fillibens <-cor(as.data.frame(qqnorm(x$residuals, plot = FALSE)))[1,2]
    crit = gamlss.dist::qZAGA(p = 1-p_adjust, mu = mu, nu = nu, sigma = sigma)

    return(data.frame(mu= mu, sigma = sigma, nu = nu, crit = crit, fillibens = fillibens))
  })

  null_m_params %>%
    plyr::ldply(., "data.frame") %>%
    dplyr::rename(base = `.id`) %>%
    return()
}

neg_to_zero = function(x){ifelse(x < 0, 0, x)}

### Adjust the height and percent values based on GBM and the significance
#' @importFrom dplyr filter select mutate 
#' @importFrom magrittr `%>%`
pvalue_adjust = function(sanger_df, wt, boi, motif, sample_file, critical_values, zaga_parameters, p_value){
  sanger_df %>%
    # only bases of interest
    dplyr::filter(grepl(wt, ctrl_max_base)) %>% # To only pull out the bases of interest from the motif
    dplyr::select(index, ctrl_index, ctrl_max_base, max_base, A_area:T_area, A_perc:T_perc) %>%
    # If the channel is not found in the potential edits, then it is not applicable or NA, otherwise
    # If the height of the channel is greater than the critical value for that channel,
    # return TRUE, otherwise, return FALSE
    dplyr::mutate(A_pvalue = mapply(FUN = gamlss.dist::dZAGA, x = A_area,
                             mu = zaga_parameters[1, "mu"],
                             sigma = zaga_parameters[1, "sigma"],
                             nu = zaga_parameters[1, "nu"]),
           C_pvalue = mapply(FUN = gamlss.dist::dZAGA, x = C_area,
                             mu = zaga_parameters[2, "mu"],
                             sigma = zaga_parameters[2, "sigma"],
                             nu = zaga_parameters[2, "nu"]),
           G_pvalue = mapply(FUN = gamlss.dist::dZAGA, x = G_area,
                             mu = zaga_parameters[3, "mu"],
                             sigma = zaga_parameters[3, "sigma"],
                             nu = zaga_parameters[3, "nu"]),
           T_pvalue = mapply(FUN = gamlss.dist::dZAGA, x = T_area,
                             mu = zaga_parameters[4, "mu"],
                             sigma = zaga_parameters[4, "sigma"],
                             nu = zaga_parameters[4, "nu"])) %>%
    dplyr::mutate(A_p_adjust = p.adjust(A_pvalue, "BH"),
           C_p_adjust = p.adjust(C_pvalue, "BH"),
           G_p_adjust = p.adjust(G_pvalue, "BH"),
           T_p_adjust = p.adjust(T_pvalue, "BH")) %>%
    dplyr::mutate(A_sig = if(!grepl("A", boi)) {FALSE} else {ifelse(A_pvalue <= p_value, TRUE, FALSE)},
           C_sig = if(!grepl("C", boi)) {FALSE} else {ifelse(C_pvalue <= p_value, TRUE, FALSE)},
           G_sig = if(!grepl("G", boi)) {FALSE} else {ifelse(G_pvalue <= p_value, TRUE, FALSE)},
           T_sig = if(!grepl("T", boi)) {FALSE} else {ifelse(T_pvalue <= p_value, TRUE, FALSE)}) %>%
    dplyr::mutate(A_sig_adjust = if(!grepl("A", boi)) {FALSE} else {ifelse(A_p_adjust<= p_value, TRUE, FALSE)},
           C_sig_adjust = if(!grepl("C", boi)) {FALSE} else {ifelse(C_p_adjust <= p_value, TRUE, FALSE)},
           G_sig_adjust = if(!grepl("G", boi)) {FALSE} else {ifelse(G_p_adjust <= p_value, TRUE, FALSE)},
           T_sig_adjust = if(!grepl("T", boi)) {FALSE} else {ifelse(T_p_adjust <= p_value, TRUE, FALSE)}) %>%
    ### also selecting for the sample index added 07.09.2019
    dplyr::select(index, ctrl_index, ctrl_max_base, max_base, A_area:T_area, A_perc:T_perc, A_sig:T_sig, A_sig_adjust:T_sig_adjust, A_pvalue:T_pvalue, A_p_adjust:T_p_adjust) %>%
    dplyr::mutate(motif = motif, sample_file = sample_file)
}

### Reverse complements a DNA character string

#' @importFrom Biostrings reverseComplement DNAString 
revcom = function(x){as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))}



#' @importFrom dplyr filter left_join mutate
#' @importFrom gamlss.dist dZAGA
#' @importFrom magrittr `%>%`
calculate_edit_pvalue = function(motif_part_of_sample, zaga_parameters, wt, edit, p_value){
  # this reproduces the functionality of "pvalue_adjust" from Mitch's code
  zaga_params_edit_base_only = zaga_parameters %>%
    dplyr::filter(base == edit)
  
  # get the potential edited rows
  potential_edits = motif_part_of_sample %>%
    dplyr::filter(expected_motif == wt)
  suppressMessages(
    motif_part_of_sample %>%
      dplyr::left_join(
        potential_edits %>%
          dplyr::mutate(edit_pvalue = mapply(FUN = gamlss.dist::dZAGA, x = .[[paste0(edit, "_area")]],
                                      mu = zaga_params_edit_base_only[1, "mu"],
                                      sigma = zaga_params_edit_base_only[1, "sigma"],
                                      nu = zaga_params_edit_base_only[1, "nu"])) %>% 
          dplyr::mutate(edit_padjust = p.adjust(edit_pvalue, "BH")) %>% 
          dplyr::mutate(edit_sig = edit_padjust < p_value)
      )
  )
}

#' @importFrom dplyr as_tibble
#' @importFrom sangerseqR makeBaseCalls
#' @importFrom magrittr `%>%`
make_sample_df = function(sample_sanger){
  # this creates the basic dataframe we use; it contains the peaks per base
  # for all positions in the sanger--we add trimming information later
  # the gist of getting peak data is as follows:
  # for each position, find the trace location where the peak is.
  # if a base did not have a peak, then figure out which base had the highest peak
  # and where. grab the value from that trace location for the NA bases. 
  # peakPosMatrix tells us where peaks were, if there was one. 
  
  peak_locs = sangerseqR::makeBaseCalls(sample_sanger)@peakPosMatrix 
  colnames(peak_locs) = bases
  peak_locs = peak_locs %>% 
    dplyr::as_tibble()
  peak_locs$max_base = sangerseqR::makeBaseCalls(sample_sanger)@primarySeq %>%
    as.character() %>% strsplit("") %>% {.[[1]]}
  
  ### Once this runs, it's the same as samp_peakAmpDF from before
  for (i in 1:nrow(peak_locs)){
    row = unlist(peak_locs[i,1:4])
    peak_vals = sapply(1:4, FUN = function(x){
      if (is.na(row[x])){return(NA)}
      else{
        return(sample_sanger@traceMatrix[row[x], x])
      }})
    peak_locs$A[i] = ifelse(is.na(peak_locs$A[i]), row[peak_locs$max_base[i]], peak_locs$A[i])
    peak_locs$C[i] = ifelse(is.na(peak_locs$C[i]), row[peak_locs$max_base[i]], peak_locs$C[i])
    peak_locs$G[i] = ifelse(is.na(peak_locs$G[i]), row[peak_locs$max_base[i]], peak_locs$G[i])
    peak_locs$T[i] = ifelse(is.na(peak_locs$T[i]), row[peak_locs$max_base[i]], peak_locs$T[i])
    
  }
  
  peak_locs$A_area = sample_sanger@traceMatrix[peak_locs$A, 1]
  peak_locs$C_area = sample_sanger@traceMatrix[peak_locs$C, 2]
  peak_locs$G_area = sample_sanger@traceMatrix[peak_locs$G, 3]
  peak_locs$T_area = sample_sanger@traceMatrix[peak_locs$T, 4]
  peak_locs$A_perc = peak_locs$A_area/ rowSums(peak_locs[c("A_area","C_area","G_area","T_area")])
  peak_locs$C_perc = peak_locs$C_area/ rowSums(peak_locs[c("A_area","C_area","G_area","T_area")])
  peak_locs$G_perc = peak_locs$G_area/ rowSums(peak_locs[c("A_area","C_area","G_area","T_area")])
  peak_locs$T_perc = peak_locs$T_area/ rowSums(peak_locs[c("A_area","C_area","G_area","T_area")])
  
  peak_locs$position = 1:nrow(peak_locs)
  peak_locs
}

#' @importFrom Biostrings pairwiseAlignment
is_revcom_ctrl_better = function(init_sample_seq,
                                 init_ctrl_seq){
  init_sample_seq = as.character(init_sample_seq)
  init_ctrl_seq = as.character(init_ctrl_seq)
  fwd_score = Biostrings::score(Biostrings::pairwiseAlignment(init_sample_seq, init_ctrl_seq))
  rev_score = Biostrings::score(Biostrings::pairwiseAlignment(init_sample_seq, 
                                               revcom(init_ctrl_seq)))
  return(rev_score > fwd_score)
}

#' @importFrom sangerseqR readsangerseq
is_file_ab1 = function(filepath){
  # checks if the file is .ab1, looks like a fasta, or neither
  if (!file.exists(filepath)){
    return(FALSE)
  }
  result = tryCatch({
    sangerseqR::readsangerseq(filepath)
    return(TRUE)
  },error =
    function(e){
      FALSE
    })
  return(result)
}

###
# Data retrieval and Plotting Functions
###


#' Get complete fit results
#'
#' @param fit the result of detect_edits
#' @return A dataframe with the following columns:
#' raw_sample_position - position in the sample sanger
#' passed_trimming - whether the sample sanger passed trimming (TRUE) or was low quality (FALSE)
#' motif_number - The numbered motif found in the control. -1 means not a match, 1 is the first match, 2 the second, etc.
#' edit_sig - whether the WT was significantly changed
#' control_base - The control sequence. Also serves as the motif sequence when motif_number > 0
#' sample_primary_call - The highest defined peak in the sample
#' sample_secondary_call - The second highest defined peak in the sample
#' A-T_area - The maximum sample height of the base peak
#' A-T_perc - The percent of this base peak's height over all bases
#' edit_pvalue/adjust 
#' motif_seq - the motif_seq as used. Meaning, if motif_fwd == FALSE, the rev-com of the motif you supplied
#' sample_file
#' control_file
#' @export
#' @importFrom dplyr mutate left_join
#' @import ggplot2
#' @importFrom magrittr `%>%`
results = function(fit){
  suppressMessages(
    fit$intermediate_data$raw_sample_df %>%
    dplyr::mutate(motif_seq = fit$sample_data$motif[1],
           sample_file = basename(fit$sample_data$sample_file[1]),
           control_file = basename(fit$intermediate_data$sample_alt$ctrl_file[1]),
           ) %>%
    dplyr::left_join(fit$sample_data %>%
                       dplyr::mutate(raw_sample_position = index) %>%
                       dplyr::rename(position_in_motif = target_base) %>%
                       dplyr::mutate(expected_change = fit$expected_change) %>%
                dplyr::select(raw_sample_position,position_in_motif,expected_change,
                              edit_pvalue, edit_padjust, edit_sig)
                ) %>%
    dplyr::mutate(passed_trimming = !is_trimmed) %>%
    dplyr::rename(control_base = control_primary_call) %>%
    dplyr::rename(motif_number = motif) %>%
    dplyr::select(raw_sample_position, passed_trimming, motif_number, edit_sig, control_base, 
                  expected_change, sample_primary_call, sample_secondary_call, A_area, C_area,
                  G_area, T_area, A_perc, C_perc, G_perc, T_perc, edit_pvalue, 
                  edit_padjust, motif_seq, sample_file, control_file)
  )
}

#' Plot the primary call percentage and trim / motif detection points
#'
#' @param fit the result of detect_edits
#' @return A ggplot object
#' @export
#' @importFrom dplyr mutate
#' @import ggplot2
#' @importFrom magrittr `%>%`
plot_raw_sample = function(fit){
  tmp = fit$intermediate_data$raw_sample_df 
  tmp %>%
    dplyr::mutate(max_base_percent = max_base_height / Tot.Area) %>%
    dplyr::mutate(trimmed = ifelse(is_trimmed, 100, 0)) %>%
    dplyr::mutate(motif = ifelse(motif > -1, 100, 0)) %>%
    ggplot(aes(x = index, y= 100*max_base_percent))+
    geom_bar(aes(y = trimmed), fill = "grey", stat = "identity", width = 2)+
    geom_bar(aes(y = motif),fill = "#53BCC2", stat = "identity", width = 2)+
    geom_line()+
    xlab("Position in sample Sanger sequence") +
    ylab("Percent signal of basecall")+
    ylim(0,100)+
    theme_classic(base_size = 18)+
    annotate("rect", xmin =  max(tmp$index) *  0.05, xmax = max(tmp$index) *0.1,
             ymin = 20, ymax = 30, fill = "grey")+
    annotate("text", x = max(tmp$index) *0.12, y = 25, label = "low Phred trimmed bases", hjust = 0)+
    annotate("rect", xmin =  max(tmp$index) *  0.05, xmax = max(tmp$index) *0.1,
             ymin = 8, ymax = 18, fill = "#53BCC2")+
    annotate("text", x = max(tmp$index) *0.12, y = 13, label = "motif bases", hjust = 0)
}



# Function for plotting trimmed sample data

#' Plot the primary call percentage just in the non-trimmed data
#'
#' @param fit the result of detect_edits
#' @return A ggplot object
#' @importFrom dplyr mutate filter
#' @import ggplot2
#' @importFrom magrittr `%>%`
#' @export
plot_trimmed_sample = function(fit){
  tmp = fit$intermediate_data$raw_sample_df 
  tmp  %>%
    dplyr::mutate(max_base_percent = max_base_height / Tot.Area) %>%
    dplyr::filter(!is_trimmed) %>%
    dplyr::mutate(motif = ifelse(motif > -1, 100, 0)) %>%
    ggplot(aes(x = index, y= 100 - 100*max_base_percent))+
    geom_bar(aes(y = motif),fill = "#53BCC2", stat = "identity", width = 2)+
    geom_line()+
    ylim(0,100)+
    xlab("Position in sample Sanger sequence") +
    ylab("Percent noise in WT basecall")+
    theme_classic(base_size = 18)+
    annotate("rect", xmin =  max(tmp$index) *  0.05, xmax = max(tmp$index) *0.1,
             ymin = 8, ymax = 18, fill = "#53BCC2")+
    annotate("text", x = max(tmp$index) *0.12, y = 13, label = "motif bases", hjust = 0)
}

#' Plot statistics in the tested wt bases 
#'
#' @param fit the result of detect_edits
#' @return A ggplot object
#' @export
#' @importFrom dplyr mutate filter
#' @import ggplot2
#' @importFrom magrittr `%>%`
plot_editing_barplot = function(fit){
  editing_data = fit$intermediate_data$output_sample_alt
  editing_data %>%
    ggplot(aes(x = sig, y = perc, color = sig, fill = sig)) +
    scale_y_continuous(limits = c(0,100), breaks = seq(0,100,10)) +
    geom_bar(data = . %>%
               group_by(sig, base) %>%
               summarize(mean_height = mean(perc)),
             aes(y = mean_height),
             stat = "identity",alpha = 0.6, color = "black") +
    geom_point(size = 2, alpha = 1) +
    ylab("Percent height of Sanger trace") +
    xlab("") +
    guides(fill = NULL) +
    theme_classic(base_size = 18) +
    theme(aspect.ratio = 4/1,
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +#sum(grepl(edit, bases))) +
    facet_wrap(.~base, nrow = 1, scales = "free_x")
}


#' @importFrom magrittr `%>%`
#' @importFrom dplyr group_by summarize rename
tableEditingData = function(editing_data) {
  editing_data %>%
    dplyr::group_by(base, sig) %>%
    dplyr::summarize(Mean = mean(perc), Max = max(perc), Min = min(perc), SD = sd(perc), N = length(perc)) %>%
    dplyr::rename(Base = base, Significance = sig)
}

#' @importFrom dplyr as_tibble distinct case_when summarize pull rename mutate arrange filter inner_join ungroup
#' @importFrom tidyr gather pivot_longer
#' @importFrom magrittr `%>%`
#' @import ggplot2
#' @importFrom sangerseqR peakPosMatrix makeBaseCalls traceMatrix
plot_chromatogram_at_motif = function(sanger, raw_sample_df, motif, 
                                      sanger_fwd = TRUE,
                                      motif_fwd = TRUE){
  motif_locations = raw_sample_df$raw_sample_position[raw_sample_df$motif != -1]
  
  # first collect the data which will make up the tiling
  sanger_peaks = raw_sample_df %>%
    dplyr::rename(A = A_area) %>%
    dplyr::rename(C = C_area) %>%
    dplyr::rename(G = G_area) %>%
    dplyr::rename(T = T_area) %>%
    dplyr::rename(base_called = sample_primary_call) %>%
    dplyr::rename(peak = max_base_height) %>%
    dplyr::rename(rowid = raw_sample_position) %>%
    dplyr::select(rowid, A, C, G, T, base_called, peak) %>%
    dplyr::filter(base_called %in% c("A","C","G","T")) %>%
    dplyr::distinct()
  
  
  # we use those locations to rearrange peak_mat
  sanger_peaks_at_motif = sanger_peaks %>%
    dplyr::filter(rowid %in% motif_locations) %>%
    dplyr::mutate(motif_pos = if(sanger_fwd){
      1:nrow(.) - 0.5
    }
    else{
      rev(1:nrow(.)) - 0.5
    }) %>%
    dplyr::select(-rowid, - peak) %>%
    tidyr::pivot_longer(c(A,C,`T`,G),
                 names_to = "base", values_to = "peak") %>%
    dplyr::group_by(motif_pos) %>%
    dplyr::mutate(peak_percent = as.character(round(100*peak/sum(peak)))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(base_pos = case_when(base == "A" ~ -1/8, base == "C" ~ -2/8,
                                base == "G" ~ -3/8, base == "T" ~ -4/8)) 
  
  
  # now figure out the trace_ends for each motif
  motif_ends = raw_sample_df %>%
    dplyr::filter(motif != -1) %>%
    dplyr::group_by(motif) %>%
    dplyr::summarize(motif_start = min(raw_sample_position),
              motif_end = max(raw_sample_position))
  
  # peakPosMatrix tells us where the peaks for the motif start and end
  trace_to_position_df = sangerseqR::peakPosMatrix(sangerseqR::makeBaseCalls(sanger)) %>% 
    as.data.frame() %>%
    dplyr::mutate(base_position = if(sanger_fwd){
      1:nrow(.)
    }else{
      1 + rev(1:nrow(.))
    }) %>%
    tidyr::pivot_longer(-base_position, names_to = "base", values_to = "rowid") %>%
    dplyr::filter(!is.na(rowid)) %>%
    dplyr::group_by(base_position) %>%
    dplyr::slice_head(n =1)
  
  # this holds starts and ends in trace and base positions per motif
  motif_ends$trace_start = trace_to_position_df %>%
    dplyr::filter(base_position %in% motif_ends$motif_start) %>%
    dplyr::pull(rowid)
  motif_ends$trace_end = trace_to_position_df %>%
    dplyr::filter(base_position %in% motif_ends$motif_end) %>%
    dplyr::pull(rowid)
  
  
  ## get the trace df
  full_trace_df = sangerseqR::traceMatrix(sangerseqR::makeBaseCalls(sanger)) %>%
    as.data.frame() %>%
    dplyr::rename(A = V1, C = V2, G = V3, `T` = V4) %>%
    dplyr::mutate(rowid = if(sanger_fwd){
      1:nrow(.)
    }else{
      rev(1:nrow(.))
    }) %>%
    tidyr::pivot_longer(-rowid, names_to = "base", values_to = "peak") %>%
    dplyr::mutate(base = if(sanger_fwd){
      base
    }else{
      dplyr::case_when(base == "A" ~ "T", base == "T" ~ "A", base == "C" ~ "G", base == "G" ~ "C")
    }) %>%
    dplyr::mutate(base_fill_code = case_when(base == "A" ~ "1000", base == "C" ~ "2000",
                                             base == "G" ~ "3000", base == "T" ~ "4000")) 
  
  motif_trace_df = data.frame()
  for (motif_id in unique(motif_ends$motif)){
    motif_trace_df = rbind(motif_trace_df,
                           full_trace_df %>%
                             dplyr::filter(rowid >= motif_ends$trace_start[motif_ends$motif == motif_id] -5,
                                    rowid <= motif_ends$trace_end[motif_ends$motif == motif_id] +5) %>%
                             dplyr::mutate(rowpos = nchar(motif) * (motif_id -1 ) + nchar(motif) * (rowid - min(rowid)) / (max(rowid) - min(rowid)))
                           )
  }
  motif_trace_df = motif_trace_df %>%
  dplyr::mutate(rowpos = if(sanger_fwd){
      rowpos
    }else{
      rev(rowpos)
    }) %>%
    dplyr::group_by(base) %>%
    dplyr::mutate(trace = peak / max(peak)+1/16) %>%
    dplyr::ungroup()

  
  # establish the colors in the fill_key to use a character set for the 
  # tile AND trace colors
  fill_colors = c(rep("#FFFFFF",7),
                  colorRampPalette(colors = c("white", "#e41a1c"))(13),
                  colorRampPalette(colors = c("#e41a1c", "#e41a1c"))(30),
                  colorRampPalette(colors = c("#e41a1c", "#377eb8"))(51),
                  "#32CD32", "#4169E1", "#121212", "#FC4338")
  names(fill_colors) = c(0:6, 7:19, 20:49, 50:100, 1000, 2000, 3000, 4000)
  
  #### We plot the raw trace which has been subsetted to just the motif locations
  
  # we modify the labels depending on how many motifs were found
  motif_order = 1:nchar(motif)
  motif_current = motif
  if (!motif_fwd){
    motif_current = revcom(motif)
    motif_order = rev(motif_order)
  }
  motif_label = strsplit(
    paste0(rep(motif_current, nrow(motif_ends)), collapse = ""), "")[[1]]
  motif_pos_label = rep(motif_order, nrow(motif_ends))
  motif_x = 1:length(motif_pos_label)-0.5
  motif_ends$motif_colors = c("lightgray","white")[1 + motif_ends$motif %% 2]
  
  # put the pieces together
  motif_trace_df %>%
    ggplot(aes(x = rowpos, ymin = 1/16, ymax = trace, 
               color =base, fill = base_fill_code))+
    geom_ribbon(alpha = 0.1)+
    scale_color_manual(values = c("A" = "#4daf4a", "C" = "#377eb8", "G" = "black", "T" = "#e41a1c"))+
    geom_tile(data = sanger_peaks_at_motif,
              inherit.aes = FALSE,
              aes(x = motif_pos, y = base_pos,  fill = peak_percent),
              color = "black")+
    geom_text(data = sanger_peaks_at_motif,
              inherit.aes = FALSE, color = "black",
              aes(x = motif_pos, y = base_pos, label = peak_percent))+
    scale_fill_manual(values = fill_colors)+
    theme_void()+
    theme(legend.position = "none",
          aspect.ratio = 1/1.5)+
    annotate("text", x = -0.5, y = -c(1:4)/8, label = c("A","C","G","T"))+
    annotate("text", x = motif_x, y = -0, 
             label = motif_label)+
    annotate("rect", xmin = (motif_ends$motif-1)*nchar(motif), 
             xmax = (motif_ends$motif)*nchar(motif), 
             ymin = rep(-0.575, nrow(motif_ends)), 
             ymax = rep(-0.675, nrow(motif_ends)),
             fill = motif_ends$motif_colors)+
    annotate("text", x = motif_x, y = -5/8, 
             label = motif_pos_label
         )+
    annotate("text", x = 1, y = 1, 
             label = paste0("motif given: ",motif, "\nmotif revcom below?  ", !motif_fwd), 
             hjust = 0, vjust = 0, size = 2)
}

#' Plot the chromatogram + percentage tiles for all detected motifs in the sample
#' 
#'
#' @param fit the result of detect_edits
#' @return A ggplot object
#' @export
plot_sample_chromatogram = function(fit){
  motif = fit$motif
  motif_fwd = fit$motif_fwd
  plot_chromatogram_at_motif(fit$sample_sanger, raw_sample_df = fit$intermediate_data$raw_sample_df, 
                             motif = motif, motif_fwd = motif_fwd)
}

#' Plot the chromatogram + percentage tiles for all detected motifs in the control
#' 
#'
#' @param fit the result of detect_edits
#' @return A ggplot object
#' @export
#' @import ggplot2
#' @importFrom magrittr `%>%`
#' @importFrom dplyr as_tibble group_by ungroup select left_join distinct case_when summarize pull rename mutate arrange filter inner_join ungroup
plot_control_chromatogram = function(fit){
  if (str_detect(fit$intermediate_data$sample_alt$ctrl_file[1], "\\.fa")){
    message("control was not .ab1, returning empty ggplot")
    return(ggplot()+
             theme_void()+
             theme(aspect.ratio = 1/1.5)+
             annotate("text", x = 0, y = 0, label = "control not .ab1"))
  }
  motif = fit$motif
  motif_fwd = fit$motif_fwd
  control_fwd = !fit$ctrl_is_revcom
  sanger = fit$ctrl_sanger
  ctrl_df = suppressWarnings(
    make_sample_df(sanger)%>% 
      dplyr::rename(raw_control_position = position) %>%
      dplyr::mutate(Tot.Area = A_area+ C_area+G_area+ T_area) %>%
      dplyr::group_by(raw_control_position) %>%
      dplyr::mutate(max_base_height = max(A_area, C_area, G_area, T_area)) %>%
      dplyr::ungroup() %>%
      dplyr::select(c(raw_control_position, max_base, A_area, C_area, G_area, T_area, 
             Tot.Area, max_base_height,
             A_perc, C_perc,G_perc,T_perc)) %>%
      dplyr::left_join(fit$intermediate_data$raw_sample_df %>%
                select(-c(A_area, C_area, G_area, T_area, Tot.Area,
                          A_perc, C_perc, G_perc, T_perc, max_base, max_base_height))) %>%
      dplyr::mutate(index = raw_control_position) %>%
      dplyr::mutate(raw_sample_position = raw_control_position)
  )
  
  
  plot_chromatogram_at_motif(fit$ctrl_sanger, raw_sample_df = ctrl_df,
                             motif = motif,
                             sanger_fwd = control_fwd,
                             motif_fwd = motif_fwd)
}


