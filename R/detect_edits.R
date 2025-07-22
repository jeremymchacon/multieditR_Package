# detect_edits.R

###########################################################################################
# Written by Jeremy Chacon and Mitchell Kluesner
#  
# This file is part of multiEditR (Multiple Edit Deconvolution by Inference of Traces in R)
# 
# Please only copy and/or distribute this script with proper citation of 
# multiEditR publication
###########################################################################################


# edits for this version:
# 1. put control sequence positions into sample
# 2. find the motif in the control and propagate positions to sample
# 3. allow there to be multiple motif positions


## Algorithm in brief:
# 1. load sequences, rev-com control and motif if necessary
# 2. pairwise align sample and control to determine relative positions
# 3. use Mott's algorithm to determine what should be trimmed in sample and flag
# 4. Find all instances of motif in control sequence, allowing zero mismatches (because otherwise small motifs would have too many false positives)
# 5. extract the trimmed, non-motif parts of the sample and use it to build the null distribution with ZAGA
# 6. calculate the p-value for all wt -> edit transitions using the ZAGA statistics and BH adjustments
#' the primary function to run edit detection. 
#' 
#'
#' @param sample_file path to the sample .ab1 file
#' @param ctrl_file path to the control .ab1 or .fa(sta) file
#' @param motif string holding the motif of interest
#' @param motif_fwd boolean, true if the motif should be in the control as-is. False if the motif requires rev-com.
#' @param wt the base to be tested for editing in motif locations
#' @param edit the predicted edit for the base
#' @param phred_cutoff the cutoff for trimming. lower is more stringent. default = 0.001
#' @param p_value the cutoff for BH-adjusted significance. 
#' @return A multiEditR object (list of lists)
#' @export
#' @importFrom magrittr `%>%`
#' @importFrom dplyr group_by ungroup select left_join distinct case_when summarize pull n rename mutate arrange filter inner_join ungroup
#' @importFrom sangerseqR readsangerseq makeBaseCalls
#' @importFrom Biostrings countPattern matchPattern pairwiseAlignment alignedSubject alignedPattern
#' @examples
#' 
#' sample_file = system.file("extdata", "RP272_cdna_wt.ab1", package="multiEditR")
#' ctrl_file = system.file("extdata", "RP272_cdna_ko.ab1", package="multiEditR")
#' motif = "AGTAGCTGGGATTACAGATG"
#' fit = detect_edits(sample_file, ctrl_file,
#'             motif = motif, motif_fwd = TRUE,
#'             wt = "A", edit = "G",
#'             phred_cutoff = 0.001, p_value = 0.05)
#'             
#' tbl = results(fit)
#' writexl::write_xlsx(tbl, "my_results.xlsx") 
#' 
#' plot_sample_chromatogram(fit)
#' plot_raw_sample(fit)       
#'             
#' 
detect_edits = function(sample_file, ctrl_file, motif, motif_fwd, wt, edit,
                            phred_cutoff = 0.00001, p_value = 0.01){
  
  # get the sequences
  sample_sanger = sangerseqR::readsangerseq(sample_file)
  sample_seq = sangerseqR::makeBaseCalls(sample_sanger) %>% 
    {.@primarySeq} %>%
    as.character()
  secondary_seq = sangerseqR::makeBaseCalls(sample_sanger) %>% 
    {.@secondarySeq} %>%
    as.character()
  # get the control sequence and put it into "ctrl_seq"
  if (is_file_ab1(ctrl_file)){
    ctrl_sanger = sangerseqR::readsangerseq(ctrl_file)
    base_calls = sangerseqR::makeBaseCalls(ctrl_sanger)
    ctrl_seq = base_calls@primarySeq %>% as.character()
  }else{
    ctrl_sanger = NA
    fasta_lines = read_lines(ctrl_file)
    ctrl_seq = paste0(fasta_lines[2:length(fasta_lines)], collapse = "")
  }
  
  # figure out if the control should be rev-com
  ctrl_is_revcom = FALSE
  if (is_revcom_ctrl_better(sample_seq, ctrl_seq)){
    message("control sequence aligns better to sample sequence when rev-com. Applying revcom to control sequence.")
    ctrl_seq = revcom(ctrl_seq)
    ctrl_is_revcom = TRUE
  }
  
  # revcom motif if necessary
  motif_orig = motif
  if (!motif_fwd){
    motif = revcom(motif)
  }
  
  # make sure at least one motif is findable in the control
  n_motif_alignments_to_ctrl = Biostrings::countPattern(pattern = Biostrings::DNAString(motif), 
                                            subject =Biostrings::DNAString(ctrl_seq),
                                            max.mismatch = 1)
  if (n_motif_alignments_to_ctrl == 0){
    stop("motif not found in control sequence while allowing 1 mismatch. are you sure you have motif_fwd correct?")
  }
  
  
  #### this will hold our main sequence table
  sample_df = make_sample_df(sample_sanger) %>%
    dplyr::rename(raw_sample_position = position)
  
  
  # apply the control sequence to the sample sequence. 
  control_alignment = Biostrings::pairwiseAlignment(pattern = DNAString(ctrl_seq), 
                                        subject = DNAString(sample_seq))
  
  aligned_sample = as.character(Biostrings::alignedSubject(control_alignment))
  aligned_sample = strsplit(aligned_sample, "")[[1]]
  aligned_control = as.character(Biostrings::alignedPattern(control_alignment))
  aligned_control = strsplit(aligned_control, "")[[1]]
  
  # figure out the raw sample position
  raw_sample_position = cumsum(aligned_sample != "-")
  raw_control_position = cumsum(aligned_control != "-")
  
  # Create a dataframe with the alignment and relative positions
  alignment_df = data.frame(raw_sample_position,
                            sample_primary_call = aligned_sample, 
                            raw_control_position,
                            control_primary_call = aligned_control, 
                            stringsAsFactors = FALSE)
  
  # bind the sample_df, which holds quality scores and percentages
  alignment_df = suppressMessages(
    dplyr::left_join(alignment_df, sample_df)
  )
  # also bind the secondary call in the sample
  secondary_call = data.frame(sample_secondary_call = strsplit(secondary_seq, "")[[1]],
                              raw_sample_position = 1:nchar(secondary_seq))
  alignment_df = suppressMessages(
    dplyr::left_join(alignment_df, secondary_call)
  )
  
  
  # find out what should be trimmed using Mott's algo, which Mitch implemented
  trim_points = get_trim_points(sample_file,cutoff = phred_cutoff)
  
  if (trim_points[1] == -1){
    warning("low quality scores. Trimming would remove most of sequence. Skipping trimming.")
    start_pos = -1
    end_pos = Inf
  }else{
    start_pos = trim_points[1]
    end_pos = trim_points[2]    
  }
  
  alignment_df$trimmed = dplyr::case_when(alignment_df$raw_sample_position < start_pos ~ TRUE,
                                   alignment_df$raw_sample_position >= end_pos ~ TRUE,
                                   .default = FALSE)
  
  # find all alignments of the motif to the control sequence, now we allow zero
  # mismatches because if the motif is tiny, we don't want clutter
  
  motif_alignments_to_ctrl = Biostrings::matchPattern(pattern = DNAString(motif), 
                                          subject = DNAString(ctrl_seq), 
                                          max.mismatch = 0)
  motif_alignments = motif_alignments_to_ctrl@ranges %>% as.data.frame()
  message(paste0(nrow(motif_alignments), " alignments of motif to control sequence found."))
  
  # add in the motif positions. -1 for non-motif, otherwise counting from 1 up for each detection
  alignment_df$motif_found = -1
  for (i in 1:nrow(motif_alignments)){
    start = motif_alignments[i, "start"]
    end = motif_alignments[i, "end"]
    alignment_df$motif[alignment_df$raw_control_position %in% seq(from = start, to = end, by = 1)] = i
  }  
  
  
  
  motif_part_of_sample = alignment_df %>%
    filter(motif != -1)
  
  ## lets check if any of the motif is being suggested for trimming, and warn if so
  if (any(motif_part_of_sample$trimmed)){
    warning("Part of the sample sanger sequence overlapping the motif is low quality and flagged for trimming.
            It is not being trimmed during edit detection, but this may affect results")
  }
  
  # this is all just renames to be compatible with make_ZAGA_df
  motif_part_of_sample$expected_motif = motif_part_of_sample$control_primary_call 
  motif_part_of_sample$ctrl_max_base = motif_part_of_sample$expected_motif 
  motif_part_of_sample$index = motif_part_of_sample$raw_sample_position
  motif_part_of_sample$ctrl_index = motif_part_of_sample$raw_control_position
  # this should be used for generating the NULL; only keep the trimmed part
  nonmotif_part_of_sample = alignment_df %>%
    dplyr::filter(motif == -1) %>%
    dplyr::filter(!trimmed)
  
  
  # pass it to the ZAGA function
  zaga_parameters = make_ZAGA_df(nonmotif_part_of_sample, p_adjust = p_value) %>%
    dplyr::mutate(sample_file = sample_file) %>%
    dplyr::select(everything())
  
  output_stats = calculate_edit_pvalue(motif_part_of_sample,
                                       zaga_parameters, wt, edit, p_value)
  
  # rearrange results to match previous version
  sample_data = suppressMessages(
    output_stats %>%
      dplyr::left_join(alignment_df %>% 
                select(-motif))  %>%
      dplyr::select(-motif) %>%
      dplyr::mutate(motif = motif) %>%
      dplyr::mutate(target_base = row.names(.)) %>% 
      dplyr::filter(!is.na(edit_sig)) %>%
      dplyr::mutate(ctrl_max_base = expected_motif) %>%
      dplyr::mutate(sample_file = sample_file) %>%
      dplyr::mutate(expected_base = expected_motif) %>%
      dplyr::mutate(sample_file = sample_file) %>%
      dplyr::mutate(passed_trimming = !trimmed) %>%
      dplyr::select(passed_trimming, target_base, motif, ctrl_max_base, expected_base, max_base, sample_secondary_call,
           A_perc, C_perc, G_perc, T_perc, 
           edit_pvalue, edit_padjust, edit_sig, index, sample_file)
  )
  
  if (!motif_fwd){
    sample_data$target_base = max(sample_data$target_base) - (sample_data$target_base - 1)
  }
  
  raw_sample_df = alignment_df %>%
    dplyr::mutate(Tot.Area = A_area + C_area + G_area + T_area) %>%
    dplyr::group_by(raw_sample_position) %>%
    dplyr::mutate(max_base_height = max(A_area, C_area, G_area, T_area)) %>% 
    dplyr::mutate(index = raw_sample_position) %>%
    dplyr::mutate(is_trimmed = trimmed) %>%
    ungroup() %>%
    select(raw_sample_position, raw_control_position, sample_primary_call,
           control_primary_call, sample_secondary_call, motif, is_trimmed, 
           A_area, C_area, G_area, T_area, Tot.Area, 
           A_perc, C_perc, G_perc, T_perc, 
           index, max_base, max_base_height) %>%
    ungroup()
  
  output_sample_alt = suppressMessages(
    output_stats %>%
      dplyr::select(-motif) %>%
      dplyr::left_join(alignment_df) %>%
      dplyr::mutate(perc = 100 * .[[paste0(edit, "_perc")]]) %>%
      dplyr::mutate(index = raw_sample_position) %>%
      dplyr::mutate(base = edit) %>%
      dplyr::mutate(sig = ifelse(edit_sig, "Significant", "Non-significant")) %>%
      dplyr::group_by(sig) %>%
      dplyr::mutate(tally = dplyr::n()) %>%
      dplyr::select(index, base, perc, sig, tally) %>%
      dplyr::filter(!is.na(sig))
  )
  
  
  
  motif_positions = motif_part_of_sample %>%
    dplyr::mutate(ctrl_post_aligned_index = index) %>%
    dplyr::mutate(motif_id = 1) %>% # relic from when > 1 motif could be found
    dplyr::select(ctrl_post_aligned_index, motif_id)
  
  output = list(
    "sample_data" = sample_data,
    "statistical_parameters" = zaga_parameters,
    "sample_sanger" = sample_sanger,
    "sample_fastq" = sample_seq,
    "ctrl_sanger" = ctrl_sanger,
    "ctrl_fastq" = ctrl_seq,
    "ctrl_is_revcom" = ctrl_is_revcom, 
    "motif" = motif_orig,
    "motif_fwd" = motif_fwd,
    "expected_change" = edit,
    "intermediate_data" = list("raw_sample_df"=raw_sample_df,
                               "sample_alt"=motif_part_of_sample %>%
                                 filter(expected_motif == wt) %>%
                                 mutate(ctrl_file = ctrl_file) %>%
                                 mutate(sample_file = sample_file),
                               "output_sample_alt" = output_sample_alt,
                               "motif_positions" = motif_part_of_sample$position)
  )
  class(output) = "multieditR"
  
  output
  
}
