#!/usr/local/bin/Rscript
# CAD_mnus421_0002
# RNASeq for mouse Control vs No T follicular helper cells.
# 11 samples 5 WT and 6 NT
#
# Analysis Performed by Xiaohui Zhao
# School of Medicine, University of Cambridge, UK
# Copyright Xiaohui Zhao (xz289@cam.ac.uk)
#
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------

message("+-------------------------------------------------------------------------------+")
message("+  Generate the full sample and samplesheet for Nextflow input                  +")
message("+-------------------------------------------------------------------------------+")

setwd("/rds/project/rds-OsGdhNr6c2U/CAD-Projects/mnus421/CAD_mnus421_0002/Original_data")

SamTab     <- read.csv("SLX-20245.HH5FVDRXY.s_2.contents.csv", header=T)
SamTab$Sample.name[10] <- "GTC246-26-NT4-10"
SamTab$Sample.name[11] <- "GTC246-27-WT3-11"
Sex        <- rep(c("Male", "Female"), c(4,7))
Cond       <- unlist(lapply(SamTab$Sample.name, function(x) strsplit(x, "-")[[1]][3]))
Replicates <- c(1,2,1,2,3,3,4,4,5,6,5)
Condition  <- gsub("[1234]", "", Cond)
Condition  <- ifelse(Condition=="N", "NT", Condition)
fastq_1    <- list.files(path=".", pattern="*.r_1.fq.gz")
fastq_2    <- list.files(path=".", pattern="*.r_2.fq.gz")
seqID_1    <- unlist(lapply(fastq_1, function(x) strsplit(x, split="[.]")[[1]][2]))
seqID_1    <- gsub("_", "-", seqID_1)
seqID_2    <- unlist(lapply(fastq_2, function(x) strsplit(x, split="[.]")[[1]][2]))
seqID_2    <- gsub("_", "-", seqID_2)
f1s        <- cbind(fastq_1, seqID_1); 
f2s        <- cbind(fastq_2, seqID_2);  ## check f1/f2 in the same order
fastq.files<- data.frame(fastq_1, fastq_2, seqID_1)
colnames(fastq.files) <- c("fastq_1", "fastq_2", "Barcode")
SamTall    <- cbind(SamTab, Sex, Replicates, Condition)


SamTable   <- merge(SamTall, fastq.files, by = "Barcode")


write.csv(SamTable, file = "CAD_mnus421_0002_SampleTable.csv", row.names=F, quote=T)

message("+-----  Generate the sample sheet for Nextflow input          -----+")

sample     <- paste0(SamTall$Condition, "_REP", SamTable$Replicates)
fastq_1    <- SamTable$fastq_1
fastq_2    <- SamTable$fastq_2
strandedness <- "reverse"

nfTable    <- data.frame(sample, fastq_1, fastq_2, strandedness)

write.csv(nfTable, file = "CAD_mnus421_0002_nextflow_SampleTable.csv", row.names=F)

message("+----------------------------Finish -------------------------------+")


