library(seqinr)
library(stringi)
library(tidyverse)
library(here)
library(data.table)
library(parallel)
library(XML)

source("scripts/functions.R")

# reg_pattern <-'[HGITS" "]E{4,8}[" "TS]{3,5}E{4,8}[HGITS" "]'
# reg2 <- '[" "TS]{3,5}'
reg_pattern <-'[HGITS" "]E{4,20}[" "TS]{1,5}E{4,20}[HGITS" "]'
reg2 <- '[" "TS]{1,5}'
# 
# ss_fasta <- seqinr::read.fasta("pdb_data/secondary_structure.fasta",
#                        as.string = T,
#                        forceDNAtolower = F)
# 
# ss_names <- getName(ss_fasta) %>%
#   str_subset("secstr")
# 
# ss_dt <- ss_fasta2dt(ss_fasta,
#                      ss_names)
# 
# keep_pdbs <- xmlParse("pdb_data/representatives_90redundant_7_22_2020.xml") %>%
#   xmlToList() %>%
#   unlist() %>% 
#   unname() %>%
#   str_replace_all("\\.",":")
# 
# ss_dt <- ss_dt[id %in% keep_pdbs]
# 
# ss_dt <-
#   filter_secondary_structure(
#     ss_dt,
#     reg_pattern = reg_pattern)
# 
# # dir_dt <- data.table(files = dir("pdb_data/PDBs/",
# #     pattern = "*.cif.gz"))[,pdb_id := str_remove(files, ".cif.gz")
# #                           ]
# # dirz <- ss_dt[!(pdb_id %in% dir_dt$pdb_id)]
# 
# # keep_track <-
# #   mclapply(
# #     unique(dirz$pdb_id),
# #     dl_pdbs,
# #     "pdb_data/PDBs",
# #     mc.cores = 40)
# 
# ss_dt <- mclapply(1:nrow(ss_dt),
#                sec_str_location,
#                ss_dt,
#                secstr = "sec_str",
#                amino_acid = "aa_seq",
#                reg_pattern = reg_pattern,
#                "full",
#                use_regex = TRUE,
#                trim_end = 1,
#                trim_start = 1,
#                index_col = TRUE, mc.cores = 20) %>%
#   data.table::rbindlist()
# 
# ss_dt <-
#   mclapply(1:nrow(ss_dt),
#            sec_str_location,
#            ss_dt,
#            secstr = "full_ss",
#            amino_acid = "full_aa",
#            reg_pattern = reg2,
#            "subseq",
#            use_regex = TRUE,
#            trim_end = 0,
#            trim_start = 0,
#            index_col = TRUE,
#            turn = TRUE,
#            mc.cores = 24) %>%
#   rbindlist()
# ss_dt <- ss_dt[!duplicated(full_aa)]
# fwrite(ss_dt,"data_mid/seqs_of_interest_bigger.csv")
ss_dt <- fread("data_mid/seqs_of_interest_bigger.csv")
ss_dt <- ss_dt[ ,!c("aa_seq",
                    "sec_str")]
rowsers <- c(1:nrow(ss_dt))
# rowsers <- c(1:3000)

system.time({
  dt_back <- 
  mclapply(rowsers,
         count_contacts,
         ss_dt,
         "subseq_start",
         "subseq_stop",
         aa_seqz = "full_aa",
         "full_ss",
         "pdb_id",
         "subseq_Beta1",
         "subseq_Beta2",
         "chain_id",
         "pdb_data",
         mc.cores = 50) 
  
  # dt_back <-
  #   lapply(c(1000:2000),
  #            count_contacts,
  #            ss_dt,
  #            "subseq_start",
  #            "subseq_stop",
  #            aa_seqz = "full_aa",
  #            "full_ss",
  #            "pdb_id",
  #            "subseq_Beta1",
  #            "subseq_Beta2",
  #            "chain_id",
  #            "pdb_data")
  
  dt_rbind <- dt_back %>%
  rbindlist(fill = TRUE)
})

fwrite(dt_rbind,"data_mid/pair_data_bigger.csv")
