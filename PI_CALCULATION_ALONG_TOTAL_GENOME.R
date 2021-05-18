#####Pp_genome_assembly
library(dplyr)
library(data.table)
library(ggplot2)

T3 <- read.table("AF_Malawi_CHIKOPA_Pp_genome_assembly_columns_passed.txt", sep = "")
print("Finished read AF_TASH_SC1")

T3 <- mutate(T3, V8 = paste(V1, V2))

T3 <- T3[,c(1,3:7)]

print("start joining data")

T3_T3 <- inner_join(T3, T3, by = "V8")
T3_T3 <- T3_T3[,c(1:6,8:11)]
rm(T3)
#T3_T3 <- T3_T3[!duplicated(T3_T3$V8), ]
#T3_T3[duplicated(T3_T3$V8), ]
end_first_join <- Sys.time()
colnames(T3_T3) <- c("CONTIG", "REF-1", "ALT-1", "GT-1","AF-1", "SHARED_POS", "REF-2", "ALT-2", "GT-2", "AF-2")
T3_T3 <- T3_T3[!(T3_T3$`GT-1` == "./." | T3_T3$`GT-2` == "./."),]
T3_T3 <- T3_T3[!(T3_T3$`GT-1` == "1/2" | T3_T3$`GT-2` == "1/2"),]

T3_T3 <- within(T3_T3, `ALT-1`[`GT-1`=="0/0"] <- `REF-1`[`GT-1`=="0/0"])
T3_T3 <- within(T3_T3, `ALT-2`[`GT-2`=="0/0"] <- `REF-2`[`GT-2`=="0/0"])

T3_T3 <- within(T3_T3, `REF-1`[`GT-1`=="1/1"] <- `ALT-1`[`GT-1`=="1/1"])
T3_T3 <- within(T3_T3, `REF-2`[`GT-2`=="1/1"] <- `ALT-2`[`GT-2`=="1/1"])

print("Measurment_GT_based_on_allele_frequency-start")
####Measurment_GT_based_on_allele_frequency

#T3_T3

T3_T3$AF_REF1 <- (1 - T3_T3$`AF-1`)
T3_T3$AF_REF2 <- (1 - T3_T3$`AF-2`)

T3_T3 <- T3_T3[, c("REF-1", "ALT-1", "AF_REF1", "AF-1", "REF-2", "ALT-2","AF_REF2", "AF-2", "SHARED_POS", "CONTIG")]

First_random_T3_T3  <- unname(apply(T3_T3[c(1:4)], 1, FUN = function(x) {return(sample(x[1:2], size = 1, prob = x[3:4]))}))
Second_random_T3_T3 <- unname(apply(T3_T3[c(5:8)], 1, FUN = function(x) {return(sample(x[1:2], size = 1, prob = x[3:4]))}))
Complete_sequnece_T3_T3 <- as.data.frame(First_random_T3_T3, stringsAsFactors = F)

rm(First_random_T3_T3)
Complete_sequnece_T3_T3$Second_population <- as.data.frame(Second_random_T3_T3, stringsAsFactors = F)
rm(Second_random_T3_T3)
Complete_sequnece_T3_T3_1 <- cbind(Complete_sequnece_T3_T3, T3_T3$CONTIG)
Complete_sequnece_T3_T3 <- Complete_sequnece_T3_T3_1 
rm(Complete_sequnece_T3_T3_1)
rm(T3_T3)
colnames(Complete_sequnece_T3_T3) <- c("First_pop", "Second_pop", "Contig")
Complete_sequnece_T3_T3$First_pop <- sapply(Complete_sequnece_T3_T3$First_pop, as.character)
Complete_sequnece_T3_T3$Second_pop <- sapply(Complete_sequnece_T3_T3$Second_pop, as.character)
Complete_sequnece_T3_T3 <- na.omit(Complete_sequnece_T3_T3)
Complete_sequnece_T3_T3$row_names <- c(1:nrow(Complete_sequnece_T3_T3))
q1 <- nrow(Complete_sequnece_T3_T3)

t3_t3 <- apply(Complete_sequnece_T3_T3[,1:4], 1, function(x) {if(nchar(x[1]) == 1 & nchar(x[2]) == 1 & x[1] != x[2]) {return(1)}else{return(0)}}) %>% sum()
average_pi <- t3_t3/q1

write(c("average_pi= ", average_pi, "t3_t3 = ", t3_t3, "q1 = ", q1), file = "PI_DIFF_TOTAL_GENOME_Pp_data.txt", append = T, sep = "\t")

t3_t3=0
w_tash3=0

windows_tash3_tash3 <- list()
for(i in 1:nrow(Complete_sequnece_T3_T3)){
  if(nchar(Complete_sequnece_T3_T3$First_pop[i]) == 1 & nchar(Complete_sequnece_T3_T3$Second_pop[i]) ==1 & Complete_sequnece_T3_T3$First_pop[i] != Complete_sequnece_T3_T3$Second_pop[i]) {
    t3_t3 = t3_t3 +1 
  }
  if(Complete_sequnece_T3_T3$row_names[i] > w_tash3){
    windows_tash3_tash3[[ as.character(w_tash3)]] <- c(t3_t3, as.character(Complete_sequnece_T3_T3$Contig[i]))
    t3_t3=0
    w_tash3=as.integer(w_tash3)+5000
  }
}

tash_tash <- transpose(as.data.frame(windows_tash3_tash3, row.names = TRUE))
tash_tash$row_names <- c(1:nrow(tash_tash))
colnames(tash_tash) <- c("Number_of_gt_diff", "contig", "row_number")
tash_tash$Number_of_gt_diff=as.integer(tash_tash$Number_of_gt_diff)
tash_tash$Number_of_gt_diff <- tash_tash[,1]/5000

write.csv(tash_tash, file = "TASH_TASH_PI_5kbp_along_all_genome.csv",  quote = FALSE)


