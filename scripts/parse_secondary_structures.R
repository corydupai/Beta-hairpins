library(seqinr)
library(stringi)
library(tidyverse)
library(here)
library(data.table)
library(parallel)
library(XML)

source("scripts/functions.R")

# Set variables.
reg_pattern <-'[HGITS" "]E{4,20}[" "TS]{1,5}E{4,20}[HGITS" "]' # regex for base Beta hairpins (4-20 strands, 1-5 turns)
reg2 <- '[" "TS]{1,5}' # regex for 1-5 turn region
sois <- "data_mid/seqs_of_interest_bigger.csv" # Save intermediate file with sequences of interest
# data_out <- "data_mid/pair_data_bigger.csv" # Save intermediate file with full summary info from PDB
data_out <- "data_mid/pair_data_bigger.csv"
angles_out <- "data_mid/angles_info.csv"

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
ss_dt <- fread(sois)
ss_dt <- ss_dt[ ,!c("aa_seq",
                    "sec_str")]
rowsers <- c(1:nrow(ss_dt))
# rowsers <- c(1:10)
# ttt <- fread(paste0("pdb_data/dssp/12E8.dssp"), header = TRUE,
#       fill = TRUE
# )[, acc := as.numeric(acc)
#   ][, Position := as.character(Position)
#     ][, Chain := as.character(Chain)]

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
         mc.cores = 25)
  
  dt_rbind <- dt_back %>%
  rbindlist(fill = TRUE)
})

angles_dt <- 
  dt_rbind[!duplicated(paste0(id,full_aa)), 
           c("id", "full_aa","subseq_aa","PHI","PSI","ALPHA","KAPPA","TCO", "touching_u")
           ][, ncs := nchar(subseq_aa)
           ][, lenz := length(unlist(PHI)), by = c("id","full_aa", "subseq_aa", "touching_u")
           ][ ncs == lenz
           ][ , list( PHI = as.numeric(unlist( PHI )),
                      PSI = as.numeric(unlist( PSI )),
                      ALPHA = as.numeric(unlist( ALPHA )),
                      KAPPA = as.numeric(unlist( KAPPA )),
                      TCO = as.numeric(unlist( TCO )),
                      Letter = as.character(unlist(str_split(subseq_aa, "")))) , by = c("id","full_aa", "subseq_aa", "touching_u") 
           ][ , t_res_num := 1
           ][, t_res_num := cumsum(t_res_num), by = c("id","full_aa", "subseq_aa", "touching_u")
           ][, turn_length := .N, by = c("id","full_aa", "subseq_aa", "touching_u")]

dt_rbind <- unique(dt_rbind[,!c("PHI","PSI","ALPHA","KAPPA","TCO")])
fwrite(angles_dt, angles_out)
fwrite(dt_rbind,data_out)
