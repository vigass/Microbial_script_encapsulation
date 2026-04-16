
library(vegan)
otu_temp <- t(otu_rel)
dis_bray <- vegdist(otu_temp, method = "bray")
dis_jaccard <- vegdist(otu_temp, method = "jaccard")


pcoa_get <- function(dis, sample_info){
  
  pcoa_res <- cmdscale(dis, k = 2, eig = TRUE)
  
  pcoa_df <- as.data.frame(pcoa_res$points)
  colnames(pcoa_df) <- c("PCoA1","PCoA2")
  pcoa_df$Sample_ID <- rownames(pcoa_df)
  
  var_exp <- round(pcoa_res$eig / sum(pcoa_res$eig) * 100, 2)
  
  pcoa_df <- pcoa_df %>% 
    left_join(sample_info, by = "Sample_ID") %>% 
    select(Sample_ID, Groups, everything())
  
  attr(pcoa_df, "variance") <- var_exp[1:2]
  
  return(pcoa_df)
}



nmds_get <- function(otu, sample_info, dist_method = "bray"){
  if (dist_method == "jaccard") {
    otu <- vegan::decostand(otu, method = "pa")
  }
  
  nmds_res <- vegan::metaMDS(
    otu, 
    distance = dist_method, 
    k = 2,
    trymax = 100
  )
  
  nmds_df <- as.data.frame(nmds_res$points)
  colnames(nmds_df) <- c("NMDS1","NMDS2")
  nmds_df$Sample_ID <- rownames(nmds_df)
  
  nmds_df <- nmds_df %>% 
    left_join(sample_info, by = "Sample_ID") %>% 
    select(Sample_ID, Groups, everything())
  
  attr(nmds_df, "stress") <- nmds_res$stress
  
  return(nmds_df)
}
