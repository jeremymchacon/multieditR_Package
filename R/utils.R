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
#' @importFrom gamlss.dist qZAGA
#' @importFrom magrittr `%>%`
make_ZAGA_df = function(sanger_df, p_adjust){
  nvals <- list()
  nvals$A = filter(sanger_df, max_base != "A")$A_area
  nvals$C = filter(sanger_df, max_base != "C")$C_area
  nvals$G = filter(sanger_df, max_base != "G")$G_area
  nvals$T = filter(sanger_df, max_base != "T")$T_area

  n_models =   n_models <-lapply(nvals, FUN = function(x){
    set.seed(1)
    if((unique(x)[1] == 0 & length(unique(x)) == 1) |
       (unique(x)[1] == 0 & length(unique(x)) == 2 & table(x)[2] == 1))
    {x = c(rep(0, 989), 0.00998720389310502, 0.00998813447664401,0.009992887520785,
           0.00999585366068316, 0.00999623914632598, 0.00999799013526835, 0.010001499423723,
           0.0100030237039207, 0.0100045782875701, 0.0100048452355807, 0.0100049548867042)
    message("Replacement vector used for low noise.")
    } # add noise if all 0s, or all 0s and one other value.
    tryCatch(gamlss::gamlss((x)~1, family = ZAGA), error=function(e) # Progressively step up the mu.start if it fails
      tryCatch(gamlss::gamlss((x)~1, family = ZAGA, mu.start = 1), error=function(e)
        tryCatch(gamlss::gamlss((x)~1, family = ZAGA, mu.start = 2), error=function(e)
          tryCatch(gamlss::gamlss((x)~1, family = ZAGA, mu.start = 3), error=function(e) # additional step added.
            tryCatch(gamlss::gamlss((x)~1, family = ZAGA, mu.start = 4), error=function(e) # additional step added.
              gamlss::gamlss((x)~1, family = ZAGA, mu.start = mean(x))
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

#' @importFrom pwalign pairwiseAlignment
is_revcom_ctrl_better = function(init_sample_seq,
                                 init_ctrl_seq){
  init_sample_seq = as.character(init_sample_seq)
  init_ctrl_seq = as.character(init_ctrl_seq)
  fwd_score = score(pwalign::pairwiseAlignment(init_sample_seq, init_ctrl_seq))
  rev_score = score(pwalign::pairwiseAlignment(init_sample_seq, 
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
# Plotting Functions
###


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

# Define geom_chromatogram function
# input the sanger object, start index, end index and output chromatogram with underlying base percents
#' @importFrom sangerseqR makeBaseCalls
#' @importFrom dplyr as_tibble as_data_frame summarize rename mutate arrange filter inner_join
#' @importFrom tibble rowid_to_column
#' @importFrom tidyr gather 
#' @importFrom magrittr `%>%`
#' @import ggplot2
geom_chromatogram = function(sanger, start_index, end_index){

  # define bases
  bases = c("A", "C", "G", "T")

  # define the basecalls  in the sanger object initially to set the proper phasing to the data
  sanger = sangerseqR::makeBaseCalls(sanger)

  # make rawtraces from the sanger trace matrix
  rawtraces = sanger@traceMatrix %>%
    # the trace matrix is the scatter of points that make up the chromatogram
    # convert the large matrix of traces into a tibble object
    dplyr::as_data_frame() %>%
    # name the columns by base channel
    dplyr::rename(A = V1, C = V2, G = V3, `T` = V4) %>%
    # create a column 1 through n trace measurements taken
    tibble::rowid_to_column("x")

  # create the peakPosMatrix giving the x coordinate of each basecall
  peakPosMatrix = sanger@peakPosMatrix %>% dplyr::as_data_frame()
  # create column names for each channel
  colnames(peakPosMatrix) = c("A", "C", "G", "T")
  # format peakPosMatrix
  peakPosMatrix = peakPosMatrix %>%
    # add the basecall for each position
    dplyr::mutate(basecall = sanger@primarySeq %>% as.character %>% strsplit(split = "") %>% unlist()) %>%
    dplyr::rowid_to_column("index") %>%
    tidyr::gather(base, x, A:`T`) %>%
    dplyr::filter(basecall == base) %>%
    dplyr::arrange(index) %>%
    dplyr::select(index, basecall, x)

  # Make the traces data frame
  traces = peakPosMatrix %>%
    # filter only the trace information at positions that are registered as the position of the peak height
    dplyr::inner_join(., rawtraces) %>%
    # define the previous base index
    dplyr::mutate(x_1 = lag(x, default = 1)) %>%
    # start the chromatogram for each base as the position between the two bases divided by two
    # start the chromatogram as the x coordinate bewteen the first  base of interest and the n-1 base
    dplyr::mutate(start = floor((x + x_1) / 2)) %>%
    # end the chromatogram for each base as the position
    # empirically derived I believe
    dplyr::mutate(end = lead(start, default = (max(x) + 6))) %>%
    # define the total area for each basecall
    dplyr::mutate(Tot.Area = A + C + G + `T`) %>%
    # calculate percent heights for each base from the total area
    dplyr::mutate(A_perc = round(A*100/Tot.Area),
                  C_perc = round(C*100/Tot.Area),
                  G_perc = round(G*100/Tot.Area),
                  T_perc = round(`T`*100/Tot.Area))

  # The position list is for every basecall, these are the corresponding x-coordinates
  # This is a list object, where the primary index of the list represents the basecall index
  # The vector element under each index represents the x-coordinates of the peakAmpMatrix for each basecall
  position_list = mapply(FUN = seq,
                         from = c(1, traces$end),
                         to = c(traces$end, max(traces$end)+10),
                         by = 1) %>%
    head(., -1)

  # make the names of the position_list indentical to the index of the position list
  names(position_list) = c(1:length(position_list))

  indices = ldply(position_list, "data.frame") %>%
    dplyr::rename(index = ".id", x = "data.frame") %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(index = as.numeric(index))

  # This joins the basecall index to the rawtrace x coordinate based data
  rawtraces %<>% dplyr::inner_join(., indices)

  # Enter the plotting based on the indices
  # defines the basecall_indices as the start and end indices of interest for making the chromatogram
  basecall_indices = c(start_index, end_index)

  # determine the start and end x coordinates of the rawtraces to filter on for plotting
  index_filters = rawtraces %>%
    # filter the rawtraces data to just include the indices that are within the region of interest
    # some of the bleedthrough is still observed at this point
    dplyr::filter(index >= basecall_indices[1] & index <= basecall_indices[2]) %>%
    # the x coordinates corresponding to the region of interest
    .$x %>%
    # takes the min and max of the peakAmpMatrix x cooridinates
    quantile(., c(0, 1)) %>%
    # adds an x coordinate cushion
    {. + c(-6, 6)}

  # gather the rawtraces data for ggplot format
  # filter in only traces that fall within the x coordinates rawtraces index filter
  # could merge this line of code with the preceeding index_filters chunk
  plotTraces = rawtraces %>%
    tidyr::gather(base, height, A:`T`) %>%
    dplyr::filter(x >= index_filters[1] & x <= index_filters[2]) %>%
    # This join operation is used to creat a base code that later allows a manual gradient to be employed
    dplyr::inner_join(.,
               data.frame(stringsAsFactors = FALSE,
                          base = c("A", "C", "G", "T"),
                          base_code = c("1000", "2000", "3000", "4000"))
    ) %>%
    # calculate the min and max height for geom_ribbon
    dplyr::mutate(max_height = pmax(height), min_height = 0)

  # calculate the y coordinates for the chromatogram
  y_bases = max(plotTraces$height)*(1/2) %>%
    {.*-c(2/11,4/11,6/11,8/11,10/11)}

  # calculate the number of bases involved
  n_bases = basecall_indices[2]-basecall_indices[1]+1

  # calculate the x coordinates for the chromatogram
  x_bases = (index_filters[2]-index_filters[1]) %>%
    {.*(seq(2, 2*n_bases+1, 2)/(2*n_bases + 2))} %>%
    {.+index_filters[1]}

  tile_plot = traces %>%
    dplyr::filter(index >= basecall_indices[1] & index <= basecall_indices[2]) %>%
    dplyr::select(basecall, A_perc, C_perc, G_perc, T_perc) %>%
    tidyr::gather(col, labels, basecall:T_perc) %>%
    dplyr::mutate(y = rep(y_bases, each = n_bases),
                  x = rep(x_bases, times = 5))

  base_annotation = data.frame(label = bases,
                               x = rep(min(plotTraces$x), 4),
                               y = tail(rev(sort(unique(tile_plot$y))), -1))

  ## Create functions for the fill color palettes
  colors1 = colorRampPalette(colors = c("white", "white"))
  colors2 = colorRampPalette(colors = c("white", "#e41a1c"))
  colors3 = colorRampPalette(colors = c("#e41a1c", "#e41a1c"))
  colors4 = colorRampPalette(colors = c("#e41a1c", "#377eb8"))

  # establish a fill_key vector
  fill_key = c(0:6, 7:19, 20:49, 50:100, 1000, 2000, 3000, 4000)

  # establish the colors in the fill_key
  fill_colors = c(colors1(length(0:6)),
                  colors2(length(7:19)),
                  colors3(length(20:49)),
                  colors4(length(50:100)),
                  "#32CD32", "#4169E1", "#121212", "#FC4338")

  # tie together the fill_key and fill_colors using the names() function
  names(fill_colors) = fill_key

  chromatogram = plotTraces %>%
    dplyr::mutate(max_height = pmax(height), min_height = 0) %>%
    ggplot(data = ., aes(x = x, y = height)) +
    geom_ribbon(aes(ymin = min_height, ymax = max_height, color = base, fill = base_code), alpha = 0.1) +
    scale_color_manual(values = c("A" = "#4daf4a", "C" = "#377eb8", "G" = "black", "T" = "#e41a1c")) +
    geom_tile(data = filter(tile_plot, col != "basecall"),
              aes(y = y, x = x, fill = as.character(as.numeric(labels))), color = "black") +
    scale_fill_manual(values = fill_colors) +
    geom_text(data = tile_plot, aes(x = x, y = y, label = labels), color = "black") +
    geom_text(data = base_annotation, aes(x = x, y = y, label = label), color = "black") +
    theme_void(base_size = 36) +
    theme(legend.position = "none",
          aspect.ratio = 1/1.5)

  return(chromatogram)
}


#' @importFrom dplyr as_tibble distinct case_when summarize pull rename mutate arrange filter inner_join ungroup
#' @importFrom tidyr gather pivot_longer
#' @importFrom magrittr `%>%`
#' @import ggplot2
#' @importFrom sangerseqR peakPosMatrix makeBaseCalls traceMatrix
plot_chromatogram_at_motif = function(sanger, raw_sample_df, motif, 
                                      sanger_fwd = TRUE,
                                      motif_fwd = TRUE){
  motif_positions = c(motif_start = min(raw_sample_df$raw_sample_position[raw_sample_df$motif == 1]),
                      motif_end = max(raw_sample_df$raw_sample_position[raw_sample_df$motif == 1])+1)
  
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
    dplyr::filter(rowid >= motif_positions[["motif_start"]],
           rowid < motif_positions[["motif_end"]]) %>%
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
  
  # peakPosMatrix tells us where the peaks for the motif start and end
  trace_ends = sangerseqR::peakPosMatrix(sangerseqR::makeBaseCalls(sanger)) %>% 
    as.data.frame() %>%
    dplyr::mutate(base_position = if(sanger_fwd){
      1:nrow(.)
    }else{
      1 + rev(1:nrow(.))
    }) %>%
    tidyr::pivot_longer(-base_position, names_to = "base", values_to = "rowid") %>%
    dplyr::filter(!is.na(rowid)) %>%
    dplyr::filter(base_position %in% c(motif_positions[["motif_start"]], 
                                motif_positions[["motif_end"]]-1)) %>%
    dplyr::pull(rowid) %>% sort
  # sometimes there are two rowid's per base position
  trace_ends = c(min(trace_ends), max(trace_ends))
  
  
  
  # establish the colors in the fill_key to use a character set for the 
  # tile AND trace colors
  fill_colors = c(rep("#FFFFFF",7),
                  colorRampPalette(colors = c("white", "#e41a1c"))(13),
                  colorRampPalette(colors = c("#e41a1c", "#e41a1c"))(30),
                  colorRampPalette(colors = c("#e41a1c", "#377eb8"))(51),
                  "#32CD32", "#4169E1", "#121212", "#FC4338")
  names(fill_colors) = c(0:6, 7:19, 20:49, 50:100, 1000, 2000, 3000, 4000)
  
  #### We plot the raw trace which has been subsetted to just the motif locations
  
  sangerseqR::traceMatrix(sangerseqR::makeBaseCalls(sanger)) %>%
    as.data.frame() %>%
    dplyr::rename(A = V1, C = V2, G = V3, `T` = V4) %>%
    dplyr::mutate(rowid = if(sanger_fwd){
      1:nrow(.)
    }else{
      1:nrow(.)
    }) %>%
    tidyr::pivot_longer(-rowid, names_to = "base", values_to = "peak") %>%
    dplyr::mutate(base = if(sanger_fwd){
      base
    }else{
      dplyr::case_when(base == "A" ~ "T", base == "T" ~ "A", base == "C" ~ "G", base == "G" ~ "C")
    }) %>%
    dplyr::mutate(base_fill_code = case_when(base == "A" ~ "1000", base == "C" ~ "2000",
                                      base == "G" ~ "3000", base == "T" ~ "4000")) %>%
    dplyr::filter(rowid >= trace_ends[1] - 5,
           rowid <= trace_ends[2] + 5) %>%
    dplyr::mutate(rowpos = nchar(motif) * (rowid - min(rowid)) / (max(rowid) - min(rowid))) %>%
    dplyr::mutate(rowpos = if(sanger_fwd){
      rowpos
    }else{
      rev(rowpos)
    }) %>%
    dplyr::group_by(base) %>%
    dplyr::mutate(trace = peak / max(peak)+1/16) %>%
    dplyr::ungroup() %>%
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
    annotate("text", x = 1:nchar(motif)-0.5, y = -0, 
             label = if(motif_fwd){
               str_split(motif, "")[[1]]
               }else{
                 str_split(revcom(motif), "")[[1]]
               })+
    annotate("text", x = 1:nchar(motif)-0.5, y = -5/8, 
             label = if(motif_fwd){
               as.character(1:nchar(motif))
             }else{
               as.character(rev(1:nchar(motif)))
             }
    )+
    annotate("text", x = 1, y = 1, label = paste0("motif given: ",motif, "\nmotif revcom below?  ", !motif_fwd), 
             hjust = 0, vjust = 0, size = 2)
}

#' @export
plot_sample_chromatogram = function(fit){
  motif = fit$motif
  motif_fwd = fit$motif_fwd
  plot_chromatogram_at_motif(fit$sample_sanger, raw_sample_df = fit$intermediate_data$raw_sample_df, 
                             motif = motif, motif_fwd = motif_fwd)
}

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


