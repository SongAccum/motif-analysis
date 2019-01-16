suppressMessages(library(dplyr))

args<-commandArgs(T)
concatenate_seq <- function(fasta_file) {
    biozhuoer::read_fasta(fasta_file) %>% 
        dplyr::mutate(name = stringr::str_extract(name, 'ENST[\\d\\.]+')) %>% 
        dplyr::group_by(name) %>% dplyr::summarise(seq = paste0(seq, collapse = '')) %>% 
        biozhuoer::write_fasta(fasta_file)
}
concatenate_seq(args[1])