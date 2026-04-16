
alpha_div_get <- function(otu, sample_info){
  alpha_df <- data.frame(
    Sample_ID = colnames(otu),
    Shannon = diversity(t(otu), index = "shannon"),
    Simpson = diversity(t(otu), index = "simpson"),
    Observed = specnumber(t(otu))
  )
  alpha_df <- left_join(alpha_df, sample_info, by = "Sample_ID") %>% 
    select(Sample_ID, Groups, everything())
  return(alpha_df)
}
