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
    mutate(is_motif = motif) %>%
    mutate(is_trimmed = trimmed) %>%
    ungroup() %>%
    select(A_area, C_area, G_area, T_area, primary_base_call, secondary_base_call,
           max_base, Tot.Area, A_perc, C_perc, G_perc, T_perc, index, max_base_height,
           is_motif, is_trimmed) %>%
    ungroup()
  
  output_sample_alt = output_stats %>%
    select(-motif) %>%
    left_join(sample_df) %>%
    mutate(perc = 100 * .[[paste0(edit, "_perc")]]) %>%
    mutate(index = position) %>%
    mutate(base = edit) %>%
    mutate(sig = ifelse(edit_sig, "Significant", "Non-significant")) %>%
    mutate(tally = row.names(.)) %>%
    select(index, base, perc, sig, tally) %>%
    filter(!is.na(sig))
  
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


calculate_edit_pvalue = function(motif_part_of_sample, zaga_parameters, wt, edit, p_value){
  # this reproduces the functionality of "pvalue_adjust" from Mitch's code
  zaga_params_edit_base_only = zaga_parameters %>%
    filter(base == edit)
  
  # get the potential edited rows
  potential_edits = motif_part_of_sample %>%
    filter(expected_motif == wt)
  
  motif_part_of_sample %>%
    left_join(
      potential_edits %>%
        mutate(edit_pvalue = mapply(FUN = gamlss.dist::dZAGA, x = .[[paste0(edit, "_area")]],
                                    mu = zaga_params_edit_base_only[1, "mu"],
                                    sigma = zaga_params_edit_base_only[1, "sigma"],
                                    nu = zaga_params_edit_base_only[1, "nu"])) %>% 
        mutate(edit_padjust = p.adjust(edit_pvalue, "BH")) %>% 
        mutate(edit_sig = edit_padjust < p_value)
    )
}


make_sample_df = function(sample_sanger){
  # this creates the basic dataframe we use; it contains the peaks per base
  # for all positions in the sanger--we add trimming information later
  # the gist of getting peak data is as follows:
  # for each position, find the trace location where the peak is.
  # if a base did not have a peak, then figure out which base had the highest peak
  # and where. grab the value from that trace location for the NA bases. 
  # peakPosMatrix tells us where peaks were, if there was one. 
  
  peak_locs = sangerseqR::makeBaseCalls(sample_sanger)@peakPosMatrix %>% 
    as_tibble()
  names(peak_locs) = bases
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

revcom = function(x){as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))}

is_revcom_ctrl_better = function(init_sample_seq,
                                 init_ctrl_seq){
  init_sample_seq = as.character(init_sample_seq)
  init_ctrl_seq = as.character(init_ctrl_seq)
  fwd_score = score(pairwiseAlignment(init_sample_seq, init_ctrl_seq))
  rev_score = score(pairwiseAlignment(init_sample_seq, 
                                      revcom(init_ctrl_seq)))
  return(rev_score > fwd_score)
}

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
