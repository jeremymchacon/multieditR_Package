# detect_edits_batch.R

###########################################################################################
# Written by Jeremy Chacon and Mitchell Kluesner
#  
# This file is part of multiEditR (Multiple Edit Deconvolution by Inference of Traces in R)
# 
# Please only copy and/or distribute this script with proper citation of 
# multiEditR publication
###########################################################################################
#' @export
#' @importFrom readxl read_excel
load_parameters_file = function(path){
  x = readxl::read_excel(path)
  colnames(x)[1:7] =  c("sample_name", "sample_file", "ctrl_file", "motif",
                  "motif_fwd","wt", "edit")
  x
}

#' load an example parameters table
#' 
#' @return A dataframe
#' @export
load_example_params = function(){
  params = load_parameters_file(paste0(system.file("extdata", package = "multiEditR"), "/parameters.xlsx"))
  ext_path = system.file("extdata", package = "multiEditR")
  params$sample_file = paste0(ext_path, "/", params$sample_file)
  params$ctrl_file = paste0(ext_path, "/", params$ctrl_file)
  params
}

#' Saves an example parameters excel spreadsheet at the given location
#' 
#' @param path A writable file path, e.g. params.xlsx
#' @export
#' @importFrom writexl write_xlsx
save_batch_skeleton = function(path){
  params = load_example_params()
  writexl::write_xlsx(params, path)
}

#' run edit detection over multiple samples
#' @param params A data.frame with the following columns: sample_file, control_file,
#' motif, motif_fwd, wt, edit, phred_cutoff, p_value. see multiEditR::detect_edits for details.
#' 
#' @return A list of multiEditR objects
#' @export
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
         
    Optional Column (must be named exactly):
        p_value : threshold for significant edits
        phred_cutoff : threshold

   The dataframe should look like this:
   
  sample_name sample_file       ctrl_file           motif                motif_fwd wt    edit  p_value
  <chr>       <chr>             <chr>               <chr>                <lgl>     <chr> <chr> <dbl>
1 Test1       RP272_cdna_wt.ab1 RP272_cdna_ko.ab1   AGTAGCTGGGATTACAGATG TRUE      A     G     0.01
2 Test2       RP272_cdna_wt.ab1 RP272_cdna_ko.fasta AGTAGCTGGGATTACAGATG TRUE      A     G     0.01
3 Test3       RP272_cdna_wt.ab1 RP272_cdna_ko.ab1   CGTATTTTTGTTAGAGATGG TRUE      A     G     0.01")
    stop("detect_edits_batch requires a parameters dataframe")
  }
  
  if (length(class(params)) == 1 && class(params) == "character"){
    message("params is a string, assuming it is a path to an xlsx sheet containing the parameters. attempting to load. ")
    params = load_parameters_file(params)
  }
  
  fits = lapply(1:nrow(params),
                FUN = function(i){
                  tryCatch({
                    fit = detect_edits(sample_file = params$sample_file[i],
                                       ctrl_file = params$ctrl_file[i],
                                       motif = params$motif[i],
                                       wt = params$wt[i],
                                       edit = params$edit[i],
                                       motif_fwd = ifelse(is.null(params$motif_fwd[i]), TRUE, params$motif_fwd[i]),
                                       p_value = ifelse(is.null(params$p_value[i]), 0.01, params$p_value[i]),
                                       phred_cutoff = ifelse(is.null(params$phred_cutoff[i]), 0.001, params$phred_cutoff[i])
                                       )
                    fit$sample_data$sample_name = params$sample_name[i]
                    fit$statistical_parameters$sample_name = params$sample_name[i]
                    fit$sample_name = params$sample_name[i]
                    fit$completed = TRUE
                    fit
                    },
                    error = function(e){
                      list(sample_name = params$sample_name[i],
                           completed = FALSE,
                           error = e)
                    })
                 })
  fits
}  

#' get a single dataframe containing all test results
#' 
#' @param fits The result of multiEditR::detect_edits_batch 
#' @return A data.frame
#' @export
#' @importFrom magrittr `%>%`
#' @importFrom dplyr mutate select tibble
#' @importFrom plyr ldply
get_batch_results_table = function(fits){
  # toss fits which failed
  fits = fits[sapply(fits, FUN = function(x){x$completed})]
  fits %>% lapply(., "[[", 1) %>%
    plyr::ldply(., "dplyr::tibble") %>%
    dplyr::mutate(target_position_in_motif = target_base) %>%
    dplyr::mutate(position_in_sample_trace = index) %>%
    dplyr::mutate(sample_max_base = max_base) %>%
    dplyr::mutate(sample_secondary_base = sample_secondary_call) %>%
    dplyr::select(sample_name, passed_trimming, target_position_in_motif, motif, 
                  expected_base, ctrl_max_base, sample_max_base, sample_secondary_base,
                  edit_sig, edit_pvalue, edit_padjust, A_perc, C_perc, G_perc, T_perc, sample_file)
}

#' get a single dataframe containing all test critical statistics
#' 
#' @param fits The result of multiEditR::detect_edits_batch 
#' @return A data.frame
#' @export
#' @importFrom magrittr `%>%`
#' @importFrom dplyr mutate select tibble
#' @importFrom tidyr pivot_wider
#' @importFrom plyr ldply
get_batch_stats_table = function(fits){
  # toss fits which failed
  fits = fits[sapply(fits, FUN = function(x){x$completed})]
  fits %>% lapply(., "[[", 2) %>%
    plyr::ldply(., "dplyr::tibble") %>%
    dplyr::mutate(base = paste0(base, " fillibens coef.")) %>%
    dplyr::select(sample_name, base, fillibens) %>%
    tidyr::pivot_wider(names_from = base, values_from = fillibens)
}

#' Creates an html report similar to the shiny app
#' 
#' @param batch_results The result of multiEditR::detect_edits_batch 
#' @param params The parameters dataframe given to multiEditR::detect_edits_batch 
#' @param path The file to write the html
#' @export
#' @importFrom rmarkdown render
create_multiEditR_report = function(batch_results, params, path = "./multiEditR_batch_results.html"){
  message("rendering document, this could take a couple of minutes. When it is done, open the html in a web browser.")
  template = paste0(system.file(package = "multiEditR"), "/batch_report_template.Rmd")
  rmarkdown::render(template,
                    params = list(params.tbl = params,
                                  results.list = batch_results),
                    output_dir = getwd(),
                    output_file = path, envir = new.env())
}
