#rm(list=ls())

#-- provide directory with htseq-count output files
dir = (".")
#-- define the pattern of files to be analysed, the file shold end as _HtSeqCount.txt
pattern="Org_.*_HtSeqCount.txt"

#-- provide metadata file containg file names and their the replicate information as below
metadata_file <- "metadata.txt"

#' Title: deseq_from_htseqcount(dir, pattern,metadata_file, "treatment_control")
#' deseq_from_htseqcount(dir, pattern,metadata_file, "treatment_control")
#' @param dir 
#' @param pattern 
#' @param metadata_file 
#'treatment_set1	treated
# treatment_set2	treated
# control_set1	untreated
# control_set2	untreated
#' @param outfile 
#' @author Pooja Sethiya (yb57662@umac.mo)
#' University of Macau
#' @return
#' @export
#'
#' @examples
deseq_from_htseqcount <- function(dir=dir, pattern, metadata_file, outfile){
          
          #--- Load package
          library(tidyverse)
          library(data.table)
          library(DESeq2)
          library(purrr)
          library(GGally)
          library(xlsx)
          #--- Load DESeq2 files
          dd <- list.files(path = dir,
                           pattern=pattern, 
                           full.names = F)
          print(dd)
          #--- get count matrix          
          df_count_matrix <- data_frame(file_name = dd) %>% 
                    mutate(file_cont = map(file_name,fread,data.table = F))  %>%
                    unnest() %>% 
                    mutate(file_name = gsub(pattern="_HtSeqCount.txt",replacement="",file_name))  %>% 
                    spread(key = file_name , value = V2) %>%
                    dplyr::slice(6:nrow(.)) #remove first sixlines of alignment summary
          
          
          #--- write count matrix in a file
          write_delim(df_count_matrix, paste(outfile, "_count_matrix.tab", sep=""),col_names = TRUE,delim="\t")
          
          #--- filter genes if read count is less than 10 in all columns
          df_count_matrix_filtered <- df_count_matrix %>% filter_at(vars(-V1), any_vars(. >= 10)) %>% as.tibble(rownames="V1")
          
          df_count_matrix_filtered <-  data.frame(df_count_matrix_filtered) 
          row.names(df_count_matrix_filtered) <- df_count_matrix_filtered$V1 
          df_count_matrix_filtered <- subset(df_count_matrix_filtered, select=-c(V1))
          
          #--- add conditions
          condition <- read_delim(metadata_file, delim="\t", col_names = FALSE)
          condition=subset(condition, condition$X1 %in% colnames(df_count_matrix_filtered))
          
          
          colData <- as.data.frame(colnames(df_count_matrix_filtered))
          colData$condition <- condition$X2[match(colData[,1], condition$X1)]
          
          colnames(colData) <- c("colData", "condition")
          print( ncol(df_count_matrix_filtered) == nrow(colData))
          
          colData$condition=factor(colData$condition, levels=c("untreated", "treated"))
          
          #--- make matrix for deseq2
          dds <- DESeqDataSetFromMatrix(countData = df_count_matrix_filtered,
                                        colData = colData,design = ~ condition)
          
          #dds <- estimateSizeFactors(dds) # run this or DESeq() first
          
          dds$condition <- factor(dds$condition)
          
          #--- Compute DESeq2
          dds  <- DESeq(dds) 
          
          #--- get DEG
          res <- results(dds)
          cat(summary(res))
          
          ##significance= rep(NA, nrow(res))
          resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
          
          names(resdata)[1] <- "Gene"
          #resdata = data.frame(resdata,significance) 
          print(head(resdata))
          
          expression_matrix <- resdata[,c(1,8:ncol(resdata))]
          print(head(expression_matrix))
          
          
          
          #--- get correlation matrix
          gg <- ggpairs(expression_matrix, columns = 2:ncol(expression_matrix), upper = list(continuous = wrap("cor", size = 5, color="black"))) +
                    theme_bw() +
                    theme(legend.text = element_text(size = 25),
                          legend.title = element_text(size = 20),         
                          axis.title.x = element_text(size = 15),        
                          axis.title.y = element_text(size = 15),
                          axis.text.x = element_text(size = 10),
                          axis.text.y = element_text(size = 10)) 
          
          png(paste(outfile, "_correlation1_plot.png", sep=""),width = 1200, height = 1200, units = "px", pointsize = 12)
          print(gg )
          dev.off()
          #--- corrplot
          rcorr = cor(as.matrix(expression_matrix[,2:ncol(expression_matrix)]))
          max <- max(rcorr)
          min <- min(rcorr)
          gc <- ggcorr(expression_matrix,label=TRUE,label_size = 8, hjust = 0.75, size = 5,label_color = "white",label_round = 2,layout.exp = 1)
          #limits = c(min,max),
          pdf(paste(outfile, "_correlation2_plot.pdf", sep=""),width = 12, height = 12)
          print(gc) 
          dev.off()
          print(gc)
          #--- write deseq output
          up_deg <- subset(resdata, resdata$log2FoldChange> 0.6 & pvalue <= 0.05)
          down_deg <- subset(resdata, resdata$log2FoldChange< -0.6 & pvalue <= 0.05)
          
          write.xlsx(resdata,sheetName = "deseq_output", file=paste(outfile, "_deseq_output.xlsx", sep=""),row.names = F, append=F)
          write.xlsx(up_deg,sheetName = "up_DEG",file=paste(outfile, "_deseq_output.xlsx", sep=""),row.names = F, append=T)
          write.xlsx(down_deg, sheetName ="down_DEG",file=paste(outfile, "_deseq_output.xlsx", sep=""),row.names = F, append=T)
          
}

