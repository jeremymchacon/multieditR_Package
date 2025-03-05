# detect_edits.R

###########################################################################################
# Written by Jeremy Chacon and Mitchell Kluesner
#  
# This file is part of multiEditR (Multiple Edit Deconvolution by Inference of Traces in R)
# 
# Please only copy and/or distribute this script with proper citation of 
# multiEditR publication
###########################################################################################


detect_edits = function(sample_file, ctrl_file, motif, motif_fwd, wt, edit,
                            phred_cutoff = 0.00001, p_value = 0.01){
  sample_sanger = sangerseqR::readsangerseq(sample_file)
  sample_seq = makeBaseCalls(sample_sanger) %>% 
    {.@primarySeq} %>%
    as.character()
  secondary_seq = makeBaseCalls(sample_sanger) %>% 
    {.@secondarySeq} %>%
    as.character()
  # get the control sequence and put it into "ctrl_seq"
  if (is_file_ab1(ctrl_file)){
    ctrl_sanger = sangerseqR::readsangerseq(ctrl_file)
    base_calls = makeBaseCalls(ctrl_sanger)
    ctrl_seq = base_calls@primarySeq %>% as.character()
  }else{
    ctrl_sanger = NA
    fasta_lines = read_lines(ctrl_file)
    ctrl_seq = paste0(fasta_lines[2:length(fasta_lines)], collapse = "")
  }
  ctrl_is_revcom = FALSE
  if (is_revcom_ctrl_better(sample_seq, ctrl_seq)){
    message("control sequence aligns better to sample sequence when rev-com. Applying revcom to control sequence.")
    ctrl_seq = revcom(ctrl_seq)
    ctrl_is_revcom = TRUE
  }
  
  # by now, the sample and control should be in the same order. 
  motif_orig = motif
  if (!motif_fwd){
    motif = revcom(motif)
  }
  
  # lets find the motif in both of them; first check.
  n_motif_alignments_to_ctrl = countPattern(pattern = Biostrings::DNAString(motif), 
                                            subject =Biostrings::DNAString(ctrl_seq),
                                            max.mismatch = 1)
  n_motif_alignments_sample = countPattern(Biostrings::DNAString(motif), 
                                           Biostrings::DNAString(sample_seq), 
                                           max.mismatch = 4)
  if (n_motif_alignments_to_ctrl == 0){
    stop("motif not found in control sequence. are you sure you have motif_fwd correct?")
  }
  if (n_motif_alignments_sample == 0){
    stop("motif not found in sample sequence with up to 4 allowed mismatches. Are you sure you have motif_fwd correct?")
  }
  
  
  # lets align by grabbing just the motif section
  motif_alignment_to_ctrl = matchPattern(pattern = DNAString(motif), 
                                         subject = DNAString(ctrl_seq),
                                         max.mismatch = 1)
  motif_alignment_to_sample = matchPattern(pattern = DNAString(motif), 
                                           subject = DNAString(sample_seq), 
                                           max.mismatch = 4)
  
  #### this will hold our main sequence table
  sample_df = make_sample_df(sample_sanger)
  
  ### lets add in where the motif is
  sample_df$motif = case_when(sample_df$position < motif_alignment_to_sample@ranges@start ~ FALSE,
                              sample_df$position >= motif_alignment_to_sample@ranges@start +
                                motif_alignment_to_sample@ranges@width ~ FALSE,
                              .default = TRUE)
  ## and the secondary call
  sample_df$secondary_base_call = strsplit(secondary_seq, "")[[1]]
  
  
  # sample trimming should happen here (because we need the trimmed bit)
  # to get the null distribution. We do this by trimming with abif_to_fastq,
  # then finding the alignment, and consider the parts outside, "trimmed"
  # Note: a more stringent trimming threshold is better, because a weak one
  # can allow too many N's, leading to failed alignment of sample to trimmed sample
  abif = abif_to_fastq("sample", sample_file,cutoff = phred_cutoff)
  
  alignment <- pairwiseAlignment(pattern = DNAString(abif$seq), 
                                 subject = DNAString(sample_seq), type = "local")
  
  # Extract start and end positions
  start_pos <- start(subject(alignment))
  end_pos <- end(subject(alignment))
  
  sample_df$trimmed = case_when(sample_df$position < start(subject(alignment)) ~ TRUE,
                                sample_df$position >= end(subject(alignment)) ~ TRUE,
                                .default = FALSE)
  
  
  motif_part_of_sample = sample_df %>%
    filter(motif)
  
  ## lets check if any of the motif is being suggested for trimming, and warn if so
  if (any(motif_part_of_sample$trimmed)){
    warning("Part of the sample sanger sequence overlapping the motif is low quality and flagged for trimming.
            It is not being trimmed during edit detection, but this may affect results")
  }
  
  motif_part_of_sample$expected_motif = strsplit(motif, "")[[1]]
  motif_part_of_sample$ctrl_max_base = motif_part_of_sample$expected_motif
  motif_part_of_sample$index = motif_part_of_sample$position
  motif_part_of_sample$ctrl_index = motif_part_of_sample$position
  # this should be used for generating the NULL; only keep the trimmed part
  nonmotif_part_of_sample = sample_df %>%
    filter(!motif) %>%
    filter(!trimmed)
  
  
  # pass it to the ZAGA function
  zaga_parameters = make_ZAGA_df(nonmotif_part_of_sample, p_adjust = p_value) %>%
    dplyr::mutate(sample_file = sample_file) %>%
    dplyr::select(everything())
  
  output_stats = calculate_edit_pvalue(motif_part_of_sample,
                                       zaga_parameters, wt, edit, p_value)
  
  # rearrange results to match previous version
  sample_data = output_stats %>%
    left_join(sample_df %>% 
                select(-motif))  %>%
    select(-motif) %>%
    mutate(motif = motif) %>%
    mutate(target_base = row.names(.)) %>% 
    filter(!is.na(edit_sig)) %>%
    mutate(ctrl_max_base = expected_motif) %>%
    mutate(sample_file = sample_file) %>%
    mutate(expected_base = expected_motif) %>%
    mutate(sample_file = sample_file) %>%
    select(target_base, motif, ctrl_max_base, expected_base, max_base, 
           A_perc, C_perc, G_perc, T_perc, 
           edit_pvalue, edit_padjust, edit_sig, index, sample_file)
  
  raw_sample_df = sample_df %>%
    mutate(Tot.Area = A_area + C_area + G_area + T_area) %>%
    mutate(primary_base_call = max_base) %>%
    group_by(position) %>%
    mutate(max_base_height = max(A_area, C_area, G_area, T_area)) %>% 
    mutate(index = position) %>%
    ungroup() %>%
    select(A_area, C_area, G_area, T_area, primary_base_call, secondary_base_call,
           max_base, Tot.Area, A_perc, C_perc, G_perc, T_perc, index, max_base_height) %>%
    ungroup()
  
  output_sample_alt = output_stats %>%
    select(-motif) %>%
    left_join(sample_df) %>%
    mutate(perc = 100 * .[[paste0(edit, "_perc")]]) %>%
    mutate(index = position) %>%
    mutate(base = edit) %>%
    mutate(sig = ifelse(edit_sig, "Significant", "Non-significant")) %>%
    mutate(tally = row.names(.)) %>%
    select(index, base, perc, sig, tally)
  
  motif_positions = motif_part_of_sample %>%
    mutate(ctrl_post_aligned_index = index) %>%
    mutate(motif_id = 1) %>% # relic from when > 1 motif could be found
    select(ctrl_post_aligned_index, motif_id)
  
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
    "intermediate_data" = list("raw_sample_df"=raw_sample_df,
                               "sample_alt"=motif_part_of_sample %>%
                                 filter(expected_motif == wt) %>%
                                 mutate(ctrl_file = ctrl_file),
                               "output_sample_alt" = output_sample_alt,
                               "motif_positions" = motif_part_of_sample$position)
  )
  class(output) = "multieditR"
  
  output
  
}

load_ctrl_seq = function(ctrl_file,
                         ctrl_is_fasta,
                         phred_cutoff){
  # Make sangerseq objects
  # Need to flesh out the TRUE statement branch
  if(ctrl_is_fasta){ 
    fasta_lines = read_lines(ctrl_file)
    input_seq = paste0(fasta_lines[2:length(fasta_lines)], collapse = "")
    init_ctrl_seq = input_seq
    ctrl_fastq = list()
    ctrl_fastq$seq = input_seq
    ctrl_df = data.frame(max_base = init_ctrl_seq %>% base::strsplit(., split = "") %>% unlist(),
                         base_call = init_ctrl_seq %>% base::strsplit(., split = "") %>% unlist()) %>%
      mutate(index = 1:NROW(max_base))
    ctrl_sanger = NULL
  } else{
    ctrl_sanger = readsangerseq(ctrl_file)
    ctrl_df = make_ctrl_sanger_df(ctrl_sanger)
    init_ctrl_seq = ctrl_df$base_call %>% paste0(., collapse = "")
    ctrl_fastq = abif_to_fastq(path = ctrl_file, cutoff = phred_cutoff)
  }
  return(list(
    "init_ctrl_seq" = init_ctrl_seq,
    "ctrl_fastq" = ctrl_fastq,
    "ctrl_df" = ctrl_df,
    "ctrl_sanger" = ctrl_sanger
  ))
}

#
# reverse_ctrl_df rev-coms as many of the objects in the ctrl_df,
# generated by load_ctrl_df, as possible, then returns the df.
# note that it doesn't reverse the sanger seq as this doesn't seem possible
# rather, in use, the sanger is often set to NULL after running this,
# so that a confusing chromatogram is not made with this.
#
reverse_ctrl_df = function(ctrl_df){
  # Note: this doesn't reverse the actual sanger seq, as I don't 
  # think that is easily possible
  ctrl_df$init_ctrl_seq = revcom(ctrl_df$init_ctrl_seq)
  ctrl_df$ctrl_fastq$seq = revcom(ctrl_df$ctrl_fastq$seq)
  ctrl_df$ctrl_df$orig_index = ctrl_df$ctrl_df$index
  ctrl_df$ctrl_df$index = (1 + max(ctrl_df$ctrl_df$index)) - ctrl_df$ctrl_df$index
  ctrl_df$ctrl_df$max_base = unname(sapply(ctrl_df$ctrl_df$max_base, FUN = revcom))
  ctrl_df$ctrl_df$base_call = unname(sapply(ctrl_df$ctrl_df$base_call, FUN = revcom))
  
  ctrl_df$ctrl_df = arrange(ctrl_df$ctrl_df, index)
  ctrl_df
}


#
# This function takes in two character strings, or DNAString objects,
# and then tries to align them as-is and after rev-coming the second string.
# it returns TRUE is the alignment score is higher after rev-coming the second
# string, otherwise returns false
is_revcom_ctrl_better = function(init_sample_seq,
                                 init_ctrl_seq){
  init_sample_seq = as.character(init_sample_seq)
  init_ctrl_seq = as.character(init_ctrl_seq)
  fwd_score = score(pairwiseAlignment(init_sample_seq, init_ctrl_seq))
  rev_score = score(pairwiseAlignment(init_sample_seq, 
                                      revcom(init_ctrl_seq)))
  return(rev_score > fwd_score)
}

##
# This function takes in a presumed file path and only returns TRUE
# if it can be loaded as a sanger sequence. it return FALSE if the file
# couldn't be found or if it couldn't be loaded as a sanger sequence.
#
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

