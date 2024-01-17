\name{detect_edits_batch}
\alias{detect_edits_batch}
\title{detect_edits_batch}
\usage{
detect_edits_batch(parameters_dataframe)

detect_edits_batch must get a params dataframe, or the path to an xlsx file which can be coerced to the correctly formatted data.frame. 
    
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
3 Test3       RP272_cdna_wt.ab1 RP272_cdna_ko.ab1   CGTATTTTTGTTAGAGATGG TRUE      A     G     TRUE         FALSE           0.01
}
\description{
Compares multiple sanger sequences to control sequences to predict edits. Comparisons are enumerated in the params. 
}
\examples{
params = load_example_params()
head(params)
fits = detect_edits_batch(params)
data.tbl = get_batch_results_table(fits)
stats.tbl = get_batch_stats_table(fits)
create_multieditR_report(fits, params, "my_html_report.html")
}
\arguments{
  \item{params}{path to the xlsx parameters sheet or a dataframe loaded from that sheet}

}
