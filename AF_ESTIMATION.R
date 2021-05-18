library(data.table)
library(tidyr)
library(dplyr)

AF_TASH3_1MLN <- fread("Malawi_ACTTGA.sort.deduplicated.columns.08.07.19.txt", sep = " ", header = F)
AF_TASH3_1MLN <- separate(AF_TASH3_1MLN, col = V6, sep = ",", into = c("DP4_FOR_REF", "DP4_REV_REF", "DP4_FOR_ALT", "DP4_REV_ALT"))
AF_TASH3_1MLN$DP4_FOR_REF <- as.numeric(as.character(AF_TASH3_1MLN$DP4_FOR_REF))
AF_TASH3_1MLN$DP4_REV_REF <- as.numeric(as.character(AF_TASH3_1MLN$DP4_REV_REF))
AF_TASH3_1MLN$DP4_FOR_ALT <- as.numeric(as.character(AF_TASH3_1MLN$DP4_FOR_ALT))
AF_TASH3_1MLN$DP4_REV_ALT <- as.numeric(as.character(AF_TASH3_1MLN$DP4_REV_ALT))
AF_TASH3_1MLN$AF <- (AF_TASH3_1MLN$DP4_FOR_ALT+AF_TASH3_1MLN$DP4_REV_ALT) / ((AF_TASH3_1MLN$DP4_FOR_REF+AF_TASH3_1MLN$DP4_REV_REF)+(AF_TASH3_1MLN$DP4_FOR_ALT+AF_TASH3_1MLN$DP4_REV_ALT))
AF_TASH3_1MLN <- AF_TASH3_1MLN[ ,c("V1", "V2", "V3", "V4", "V5", "AF")]
AF_TASH3_1MLN <- AF_TASH3_1MLN[!(AF_TASH3_1MLN$AF == "NaN"),]
write.table(x = AF_TASH3_1MLN, file = "Malawi_ACTTGA.sort.deduplicated.columns.08.07.19.af.txt", sep = "\t", quote = F, col.names = F, row.names = F)
