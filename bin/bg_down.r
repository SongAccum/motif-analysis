args<-commandArgs(T)
slide <- function(input_bed, output_bed, n = 1000) {
    col_names <- c('chr', 'start', 'end', 'name', 'score', 'strand');

    original <- readr::read_tsv(input_bed, col_names) %>% 
        dplyr::group_by_at(-2:-3) %>% 
        dplyr::summarise(length = sum(end - start), end = max(end)) %>% 
        dplyr::ungroup()
    
    if (n > 0) {
       slide <- original %>% dplyr::mutate(start = end + n, end = start + length)
    } else {
       slide <- original %>% dplyr::mutate(end = start + n, start = end - length)
    }
    
    slide %>% dplyr::select(chr, start, end, name, score, strand) %>% 
        readr::write_tsv(output_bed, col_names = F)
}
slide(args[1], args[2])