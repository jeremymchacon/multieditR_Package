
load_parameters_file = function(path){
  readxl::read_excel(path,
                     col_names = c("sample_name", "sample_file", "ctrl_file", "motif",
                                   "motif_fwd","wt", "edit", "ctrl_seq_fwd", "use_ctrl_seq",
                                   "p_value"), 
                     skip = 1)
}

load_example_params = function(){
  params = load_parameters_file(paste0(system.file("extdata", package = "multieditR"), "/parameters.xlsx"))
  ext_path = system.file("extdata", package = "multieditR")
  params$sample_file = paste0(ext_path, "/", params$sample_file)
  params$ctrl_file = paste0(ext_path, "/", params$ctrl_file)
  params
}

save_batch_skeleton = function(path){
  params = load_example_params()
  writexl::write_xlsx(params, path)
}



detect_edits_batch = function(params = NULL){
  if (is.null(params)){
    cat("detect_edits_batch must get a params dataframe, or the path to an xlsx file which can be coerced to the correctly formatted data.frame. 
    
    To save an example xlsx spreadsheet, 
         use save_batch_skeleton('path/to/save/spreadsheet')
         
    The dataframe needs the following columns:
         sample_name : 'name of this sample'
         sample_file : /full/path/to/sample/ab1
         ctrl_file : /full/path/to/control/ab1/or/fasta
         motif : sequence/to/look/for/edits/within
         wt : base to be edited
         edit: base the wt should have become
         
    The dataframe can optionally have these columns:
        ctrl_seq_fwd : boolean, whether the ctrl sequence is in forward orientation
        use_ctrl_seq : (deprecated), true if the ctrl is a fasta and not ab1
        p_value : threshold for significant edits
        phred_cutoff : threshold,
        trim : boolean, whether to trim bases below phred (TRUE by default)
        boi : bases of interest (defaults to paste0(wt, '|', edit))
        bases : bases in sequence (defaults to ACGT)
        
   The dataframe should look like this:
   
  sample_name sample_file       ctrl_file           motif                motif_fwd wt    edit  ctrl_seq_fwd use_ctrl_seq p_value
  <chr>       <chr>             <chr>               <chr>                <lgl>     <chr> <chr> <lgl>        <lgl>          <dbl>
1 Test1       RP272_cdna_wt.ab1 RP272_cdna_ko.ab1   AGTAGCTGGGATTACAGATG TRUE      A     G     TRUE         FALSE           0.01
2 Test2       RP272_cdna_wt.ab1 RP272_cdna_ko.fasta AGTAGCTGGGATTACAGATG TRUE      A     G     TRUE         TRUE            0.01
3 Test3       RP272_cdna_wt.ab1 RP272_cdna_ko.ab1   CGTATTTTTGTTAGAGATGG TRUE      A     G     TRUE         FALSE           0.01")
    stop("detect_edits_batch requires a parameters dataframe")
  }
  
  if (length(class(params)) == 1 && class(params) == "character"){
    print("params is a string, assuming it is a path to an xlsx sheet containing the parameters. attempting to load. ")
    params = load_parameters_file(params)
  }
  
  fits = lapply(1:nrow(params),
                FUN = function(i){
                  fit = detect_edits(sample_file = params$sample_file[i],
                                     ctrl_file = params$ctrl_file[i],
                                     motif = params$motif[i],
                                     wt = params$wt[i],
                                     edit = params$edit[i],
                                     motif_fwd = ifelse(is.null(params$motif_fwd[i]), TRUE, params$motif_fwd[i]),
                                     ctrl_seq_fwd = ifelse(is.null(params$ctrl_seq_fwd[i]), TRUE, params$ctrl_seq_fwd[i]),
                                     use_ctrl_seq = ifelse(is.null(params$use_ctrl_seq[i]), TRUE, params$use_ctrl_seq[i]),
                                     p_value = ifelse(is.null(params$p_value[i]), 0.01, params$p_value[i]),
                                     phred_cutoff = ifelse(is.null(params$phred_cutoff[i]), 0.001, params$phred_cutoff[i]),
                                     trim = ifelse(is.null(params$trim[i]), TRUE, params$trim[i]),
                                     boi = ifelse(is.null(params$boi[i]), paste0(params$wt[i], "|", params$edit[i]), params$boi[i]),
                                     bases = ifelse(is.null(params$bases[i]), c("A","C","G","T"), params$bases[i])
                                     )
                  fit$sample_data$sample_name = params$sample_name[i]
                  fit$statistical_parameters$sample_name = params$sample_name[i]
                  fit$sample_name = params$sample_name[i]
                  fit
                 })
  fits
}  

get_batch_results_table = function(fits){
  fits %>% lapply(., "[[", 1) %>%
    plyr::ldply(., "tibble") %>%
    select(sample_name, target_base, motif, ctrl_max_base, A_perc, C_perc, G_perc, T_perc,
           A_sig, C_sig, G_sig, T_sig, A_pvalue, C_pvalue, G_pvalue, T_pvalue, index, ctrl_index,
           sample_file)
}

get_batch_stats_table = function(fits){
  fits %>% lapply(., "[[", 2) %>%
    plyr::ldply(., "tibble") %>%
    mutate(base = paste0(base, " fillibens coef.")) %>%
    dplyr::select(sample_name, base, fillibens) %>%
    pivot_wider(names_from = base, values_from = fillibens)
}

create_multieditR_report = function(batch_results, params, path = "./multieditR_batch_results.html"){
  print("rendering document, this could take a couple of minutes. When it is done, open the html in a web browser.")
  template = paste0(system.file(package = "multieditR"), "/batch_report_template.Rmd")
  rmarkdown::render(template,
                    params = list(params.tbl = params,
                                  results.list = batch_results),
                    output_file = path, envir = new.env())
}