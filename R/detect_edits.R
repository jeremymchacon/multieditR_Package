# detect_edits.R

###########################################################################################
# Written by Jeremy Chacon and Mitchell Kluesner
#  
# This file is part of multiEditR (Multiple Edit Deconvolution by Inference of Traces in R)
# 
# Please only copy and/or distribute this script with proper citation of 
# multiEditR publication
###########################################################################################


detect_edits = function(
    # Set parameters
  sample_file,
  ctrl_file,

  motif = "A", # Use IUPAC notation
  motif_fwd = TRUE,

  wt = "A", # Enter wt bases of interest with | separation
  edit = "G", # Enter edit of interest with | separation
  p_value = 0.01,
  phred_cutoff = 0.001){
    message("initialized.")
    start = Sys.time()
    suppressWarnings({

      # previously these were parameters but seem better to keep static
      if (motif_fwd){
        boi = paste0(wt, "|", edit)
      }else{
        boi = paste0(revcom(wt), "|", revcom(edit))
      }
      bases = c("A", "C", "G", "T")
      trim = TRUE

      # check files
      if (!is_file_ab1(sample_file)){
        stop("sample_file could not be loaded as a sanger sequence. Are you sure the file exists and is an ab1-like file?")
      }
      ctrl_is_fasta = FALSE
      if (!is_file_ab1(ctrl_file)){
        ctrl_is_fasta = TRUE
      }

      # Make sangerseq objects
      ctrl_info = load_ctrl_seq(ctrl_file,
                                ctrl_is_fasta,
                                phred_cutoff)
      # temporarily read the sample sanger to check if we need to revcom
      sample_sanger = readsangerseq(sample_file)
      init_sample_seq = primarySeq(sample_sanger)
      # IF the sample needs to be reversed, then reverse complement the ctrl_seq
      ctrl_is_revcom = FALSE
      if (is_revcom_ctrl_better(init_sample_seq, ctrl_info$init_ctrl_seq)){
        ctrl_info = reverse_ctrl_df(ctrl_info)
        message("control aligns better to sample when rev-com. Applying revcom to control. Setting control sanger to NULL because it cannot be revcom")
        ctrl_is_revcom = TRUE
      }
      init_ctrl_seq = ctrl_info[["init_ctrl_seq"]]
      ctrl_df = ctrl_info[["ctrl_df"]]
      ctrl_fastq = ctrl_info[["ctrl_fastq"]]
      ctrl_sanger = ctrl_info[["ctrl_sanger"]]

      # Make sangerseq object
      # Generate samp sanger data frame
      # Generate samp primary basecalls
      sample_sanger = readsangerseq(sample_file)
      sample_df = make_samp_sanger_df(sample_sanger, init_ctrl_seq)
      init_sample_seq = sample_df$primary_base_call %>% paste0(., collapse = "")
      # Genereate phred scores for ctrl and samp, trimming is built in using mott's algorithm
      sample_fastq = abif_to_fastq(path = sample_file, cutoff = phred_cutoff)

      # Align the both the ctrl and samp to their fastq filtered sequences
      # reasonable to assume phred scored sequence will always be smaller than the primary seq
      # use high gap penalty to force SNP alignment
      sample_alignment = pairwiseAlignment(pattern = sample_fastq$seq, subject = init_sample_seq)
      ctrl_alignment = pairwiseAlignment(pattern = ctrl_fastq$seq, subject = init_ctrl_seq)

      # Save unfiltered dataframes
      raw_sample_df = sample_df
      raw_ctrl_df = ctrl_df

      message("samples loaded.")

      # Filter dfs on high phred sequence
      sample_df %<>%
        filter(index >= sample_alignment@subject@range@start) %<>%
        filter(index <= sample_alignment@subject@range@start + sample_alignment@subject@range@width - 1) %<>%
        mutate(post_filter_index = 1:NROW(index))
      ctrl_df %<>%
        filter(index >= ctrl_alignment@subject@range@start) %<>%
        filter(index <= ctrl_alignment@subject@range@start + ctrl_alignment@subject@range@width - 1) %<>%
        mutate(post_filter_index = 1:NROW(index))

      # Regenerate primary basecalls
      ctrl_seq = ctrl_df$base_call %>% paste0(., collapse = "")
      sample_seq = sample_df$max_base %>% paste0(., collapse = "")

      # Save the pre-cross alignment dataframes
      pre_cross_align_sample_df = sample_df
      pre_cross_align_ctrl_df = ctrl_df

      message("samples filtered.")
      ### Bring samples together ###
      ### 01.07.19, if a ctrl sequence is used instead of a ctrl sequence, it would enter here.
      ### To use the context correction, would you still need to have the ctrl sequence to apply the GBM to ?
      # Align sample_seq to ctrl_seq
      trimmed_alignment = align_and_trim(sample_seq, ctrl_seq, min_continuity = 15)
      
      
      
      samp_alignment_seq = trimmed_alignment$alignment@pattern %>% as.character()
      ctrl_alignment_seq = trimmed_alignment$alignment@subject %>% as.character()

      # Align the trimmed sequences to the sequences from the data frame
      sample_trimmed_alignment = pairwiseAlignment(pattern = trimmed_alignment$pattern, subject = sample_seq, gapOpening = 1000, gapExtension = 1000, type = "local")
      ctrl_trimmed_alignment = pairwiseAlignment(pattern = trimmed_alignment$subject, subject = ctrl_seq, gapOpening = 1000, gapExtension = 1000, type = "local")

      # Filter dfs to aligned sequences
      sample_df %<>%
        filter(post_filter_index >= sample_trimmed_alignment@subject@range@start) %<>%
        filter(post_filter_index <= sample_trimmed_alignment@subject@range@start + sample_trimmed_alignment@subject@range@width - 1) %>%
        mutate(post_aligned_index = 1:NROW(index))
      ctrl_df %<>%
        filter(post_filter_index >= ctrl_trimmed_alignment@subject@range@start) %<>%
        filter(post_filter_index <= ctrl_trimmed_alignment@subject@range@start + ctrl_trimmed_alignment@subject@range@width - 1) %>%
        mutate(post_aligned_index = 1:NROW(index))

      message("sample aligned.")


      ### Generate post-filter, post-aligned sequence
      ctrl_seq = ctrl_df$base_call %>% paste0(., collapse = "")
      samp_seq = sample_df$max_base %>% paste0(., collapse = "")

      ### Save a df for NGS analysis
      pre_aligned_sample_df = sample_df

      ### Filter out any base positons that have an indel
      samp_indel = samp_alignment_seq %>% gregexpr("-", .) %>% unlist
      ctrl_indel = ctrl_alignment_seq %>% gregexpr("-", .) %>% unlist

      ### Code added 10.27.18, as when there is no indels in one sample it returned a -1 instead of 0, which introduced an error
      if(length(samp_indel) == 1 && samp_indel == -1){samp_indel = 0}
      if(length(ctrl_indel) == 1 && ctrl_indel == -1){ctrl_indel = 0}

      ctrl_df = ctrl_df %>%
        mutate(.,
               indel_filter = ctrl_alignment_seq %>%
                 subchar(., samp_indel, "_") %>%
                 gsub("-", "", .) %>%
                 strsplit(x = ., split = "") %>%
                 .[[1]]
        ) %>%
        filter(indel_filter != "_") %>%
        dplyr::select(-indel_filter) %>%
        mutate(filtered_index = 1:NROW(max_base))

      sample_df = sample_df %>%
        mutate(.,
               indel_filter = samp_alignment_seq %>%
                 subchar(., ctrl_indel, "_") %>%
                 gsub("-", "", .) %>%
                 strsplit(x = ., split = "") %>%
                 .[[1]]
        ) %>%
        filter(indel_filter != "_") %>%
        dplyr::select(-indel_filter) %>%
        mutate(filtered_index = 1:NROW(max_base)) %>%
        mutate(ctrl_post_aligned_index = ctrl_df$post_aligned_index)

      message("Indels removed.")
      ### join the initial ctrl index to the sample to give a reference
      sample_df = ctrl_df %>%
        dplyr::select(post_aligned_index, index) %>%
        dplyr::rename(ctrl_post_aligned_index = post_aligned_index, ctrl_index = index) %>%
        inner_join(., sample_df) #%>%
      #dplyr::select(everything(), ctrl_index, ctrl_post_aligned_index, pre_trinucleotide, post_trinucleotide)

      ### Assign the sample and ctrl file names
      sample_df %<>% mutate(sample_file = sample_file, ctrl_file = ctrl_file)
      ctrl_df %<>% mutate(sample_file = sample_file, ctrl_file = ctrl_file)

      ### reassign the sequences post filtering
      samp_df_seq = sample_df$max_base %>% paste0(., collapse = "")
      ctrl_df_seq = ctrl_df$max_base %>% paste0(., collapse = "")

      message("Indices adjusted.")
      ### ENTER MOTIF ISOLATION ###
      # Will want to be compare this to post alignment seq and index, but before indel removal
      # ctrl_seq will have the same indexing as post_aligned_index, but the indexing is not the same across the two
      # Will need to figure out how to be able to pull the same indexing out of both of them
      # Could take the post_aligned_index from ctrl and apply it to the sample

      # Reverse complement motif if needed
      motif_orig = motif
      wt_orig = wt
      if(motif_fwd){}else{
        motif = revcom(motif)
        wt = revcom(wt)
        edit = revcom(edit)
        }

      # Align the motif of interest to the ctrl_seq and check if the motif can be found
      motif_alignment = matchPattern(pattern = DNAString(motif), 
                                     subject = DNAString(ctrl_df_seq), fixed = FALSE)
      n_alignments = motif_alignment@ranges %>% length()
      if (n_alignments == 0){
        stop("motif not found in control sequence. are you sure you have motif_fwd correct?")
      }
      # repeat for sample
      motif_samp_alignments = countPattern(DNAString(motif), 
                                           DNAString(samp_df_seq), 
                                           max.mismatch = 4)
      if (motif_samp_alignments == 0){
        stop("motif not found in sample sequence with up to 4 allowed mismatches. Are you sure you have motif_fwd correct?")
      }
      
      
      motif_positions = mapply(FUN = seq,
                               from = motif_alignment@ranges@start,
                               to = (motif_alignment@ranges@start + nchar(motif) - 1)) %>% as.vector()
      names(motif_positions) = rep(x = c(1:n_alignments), each = nchar(motif))
      if (!motif_fwd){
        motif_positions = motif_positions + 1
      }

      # Append the sequences from the ctrl df to the sample df
      sample_df %<>% mutate(ctrl_max_base = ctrl_df$max_base, ctrl_base_call = ctrl_df$base_call)

      # Generate null and alternative samples for distribution
      # Perform for both the sample df
      sample_null = sample_df %>% filter(!(ctrl_post_aligned_index %in% motif_positions)) # This dataframe consists of all bases in the sample where the motif is not found in the ctrl sequence
      sample_alt = sample_df %>% filter(ctrl_post_aligned_index %in% motif_positions) # This dataframe consists of all bases in the sample where the motif is found in the ctrl sequence


      # Find all potential events of significant noise
      filtered_sample_alt = sample_alt %>%
        #dplyr::filter(grepl(wt, ctrl_max_base)) %>% # Use the ctrl wt base for determining data
        dplyr::filter(ctrl_max_base == wt) %>%
        dplyr::rename(A = A_area, C = C_area, G = G_area, `T` = T_area) %>%
        tidyr::gather(base, height, A:`T`) %>%
        dplyr::filter(base == edit) # Filter out hypothetical mutations that are not of interest

      message("Motifs of interest mapped and subsetted.")

      # Adjust p-value
      n_comparisons = NROW(filtered_sample_alt)
      # Holm-sidak correction, may need to use the smirnov correction for family-wise error rates
      # Will need to read original paper and cite appropriately
      # if(adjust_p) {p_adjust = p.adjust(p_value, method = "holm", n = n_comparisons)} else {p_adjust = p_value}

      # Generate zG models for each base
      # uses the sample_null to calculate
      zaga_parameters = make_ZAGA_df(sample_null, p_adjust = p_value) %>%
        dplyr::mutate(sample_file = sample_file) %>%
        dplyr::select(everything())

      critical_values = zaga_parameters$crit

      ### Find significant edits and then apply the GBM adjustments
      # Determine which values are significant
      # Keep significant values and replace all n.s. values with NA
      # Use the significant values to calculated an adjusted height using formula 1
      # Find all significant edits
      suppressWarnings({
        output_sample_alt = pvalue_adjust(sample_alt, wt, boi, motif, sample_file, critical_values, zaga_parameters, p_value)
      })

      # Create Editing index for output_sample_alt
      # 02.02.2020
      output_sample_alt = output_sample_alt %>%
        dplyr::group_by(ctrl_max_base) %>% # technically shouldn't change anything as it's the same across samples
        dplyr::mutate(AEI_sanger = (sum(G_perc) / (sum(A_perc) + sum(G_perc))),
                      CEI_sanger = (sum(T_perc) / (sum(C_perc) + sum(T_perc))),
                      GEI_sanger = (sum(A_perc) / (sum(G_perc) + sum(A_perc))),
                      TEI_sanger = (sum(C_perc) / (sum(T_perc) + sum(C_perc)))
        ) %>%
        ungroup() %>%
        # Keep only the EI_index for each reference base
        tidyr::gather(EI_base, EI_sanger, AEI_sanger:TEI_sanger) %>%
        dplyr::mutate(EI_base = gsub("EI_sanger", "", EI_base)) %>%
        dplyr::filter(EI_base == ctrl_max_base) %>%
        dplyr::select(-EI_base)

      # define control chromatogram indices to set base position
      sample_chromatogram_indices = range(sample_alt$index) %>% sort
      ctrl_chromatogram_indices = range(sample_alt$ctrl_index) %>% sort

      output_sample = output_sample_alt %>%
        mutate(target_base = if (motif_fwd){
          which(str_split(motif_orig, "")[[1]] == wt_orig)
          }else{
            rev(which(str_split(motif_orig, "")[[1]] == wt_orig))
          }) %>%
        dplyr::select(target_base, `motif`, ctrl_max_base,max_base, A_perc:T_perc, A_sig:T_sig, 
                     A_pvalue:T_pvalue, index, ctrl_index, sample_file)

      output_sample_null = sample_null %>%
        dplyr::select(ctrl_index, ctrl_max_base, max_base, A_area:T_area, A_perc:T_perc) %>%
        dplyr::mutate(motif = motif, sample_file = sample_file)


      if (!exists("ctrl_sanger")){
        ctrl_sanger = NULL
      }
      pre_cross_align_sample_df <- pre_cross_align_sample_df %>%
        mutate(max_base_perc = 100*max_base_height / Tot.Area)

      output_sample_alt = output_sample_alt %>%
        dplyr::select(index, ctrl_index, A_perc:T_perc) %>%
        gather(base, perc, A_perc:`T_perc`) %>%
        inner_join(., output_sample_alt %>%
                     dplyr::select(index, A_sig:T_sig) %>%
                     gather(base_sig, sig, A_sig:T_sig)
        ) %>%
        mutate(base = gsub("_perc", "", base), base_sig = gsub("_sig", "", base_sig)) %>%
        mutate(perc = perc*100, sig = {ifelse(sig, "Significant", "Non-significant")}) %>%
        filter(base == base_sig) %>%
        filter(grepl(pattern = edit, x = base)) %>%
        mutate(tmp = 1) %>%
        group_by(base) %>%
        mutate(tally = cumsum(tmp)) %>%
        dplyr::select(-tmp, -base_sig)

      motif_positions = data.frame(ctrl_post_aligned_index = motif_positions, motif_id = names(motif_positions))

      ###
      output = list(
        "sample_data" = output_sample,
        "statistical_parameters" = zaga_parameters,
        "sample_sanger" = sample_sanger,
        "sample_fastq" = sample_fastq,
        "ctrl_sanger" = ctrl_sanger,
        "ctrl_fastq" = ctrl_fastq,
        "ctrl_is_revcom" = ctrl_is_revcom, 
        "motif" = motif_orig,
        "motif_fwd" = motif_fwd,
        "intermediate_data" = list("raw_sample_df"=raw_sample_df,
                                 "sample_alt"=sample_alt,
                                 "pre_cross_align_sample_df"=pre_cross_align_sample_df,
                                 "output_sample_alt" = output_sample_alt,
                                 "motif_positions" = motif_positions)
      )
      class(output) = "multieditR"

      # Reset start directory
      message(paste0(round(Sys.time() - start, 2), " seconds elapsed."))
    })
    return(output)

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

