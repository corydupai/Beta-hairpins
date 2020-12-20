count_contacts <- function(row_pos,
                      dt,
                      turn_start,
                      turn_end,
                      aa_seqz,
                      ss_seq,
                      ID_col,
                      B1,
                      B2,
                      chain,
                      pdb_path
                      ){
  three_lets <- c()
  indieseas <- c()
  dt <- dt[row_pos, ]

  dt[,touching_u:="NO"]
    #add variables for important info from dt
    turn1 <- dt[ , get(turn_start)]
    turn2 <- dt[ , get(turn_end)]
    full_seq <- dt[ , get(aa_seqz)]
    sec_str <- unlist(strsplit(dt[ , get(ss_seq)],""))
    B_start <- dt[ , get(B1)]
    B_end <- dt[ , get(B2)]
    chainz <- dt[ , get(chain)]
    pdb_ID <-  dt[ , get(ID_col)]

      #Read in the pdb to be cleaned/organized into a dataframe
      pdb_loc <- paste0(pdb_path,"/PDBs/",pdb_ID,".cif.gz")
      pdb_main <- read_file(pdb_loc)
      
      pdb_main <- unlist(strsplit(unlist(
        strsplit(pdb_main,
                 "_atom_site.pdbx_PDB_model_num \n"))[2],"#"))[1]
      pdb_main <- gsub("\\?[[:space:]]+\\?[[:space:]]+\\?[[:space:]]+\\?[[:space:]]+\\?[[:space:]]+\\?","?", pdb_main, fixed=FALSE)
      pdb_main <- unique(fread(pdb_main)[,!c("V1","V2","V7","V14","V15",
                                           "V16","V18","V20")])
      
      colnames(pdb_main)<-c("TypeA","TypeB","Period",
                            "AA","Chain2", "Position_true",
                            "Sub_pos","X","Y","Z","Position","Chain","Rep")
      pdb_main <- pdb_main[Rep==1 & AA!="HOH"]
      pdb_main[, Sub_pos := gsub("\\?","",Sub_pos)
             ][, Position := paste0(Position, Sub_pos)
             ][, Position_true := as.integer(Position_true)
             ][, Chain := as.character(Chain)]
      
      access <-  
        tryCatch({
          fread(paste0(pdb_path,"/dssp/", pdb_ID,".dssp"), header = TRUE
          )[, acc := as.numeric(acc)
          ][, Position := as.character(Position)
          ][, Chain := as.character(Chain)]
          # }, warning = function(w) {
          #   dt[ , touching_u:="NO DSSP"]
          #   return(dt)
          }, error = function(e) {
            dt[ , touching_u:="NO DSSP"]
            return(dt)
          }, finally = {
          })
      if(is.null(access$acc)){
        return(dt)
      }
 
      pdb_main <- merge(pdb_main,access[,c("Position","acc","Chain","TCO","KAPPA","ALPHA","PHI","PSI")],by=c("Position","Chain"))
      setorder(pdb_main, Position_true)

      #lists to convert 3 letter amino acids to 1 letter
      aa.list <- c("A","C","D","E","F","G","H","I","K","L",
                   "M","N","P","Q","R","S","T","V","W","Y","M","U",
                   "A","R","N","D","C","Q","E","H","I","L","K","M",
                   "F","P","S","T","W","Y","V","O","C","W","E","C",
                   "C","E","D","C","K","K","C","K","Y","Y","F","A",
                   "L","C","H","C","T","A","G","V","C","K","K","F",
                   "D","C","C","d")
      three.list <- c("ALA","CYS","ASP","GLU","PHE","GLY","HIS",
                      "ILE","LYS","LEU","MET","ASN","PRO","GLN",
                      "ARG","SER","THR","VAL","TRP","TYR","MSE","SEC",
                      "DAL","DAR","DSG","DAS","DCY","DGN","DGL","DHI",
                      "DIL","DLE","DLY","MED","DPN","DPR","DSN","DTH",
                      "DTR","DTY","DVA","PYL","OCS","TRO","PCA","SMC",
                      "CSO","CGU","BHD","CSD","APK","MLY","CME","KCX",
                      "IYR","TYI","PHI","MAA","MLE","SAH","AVI","CSS",
                      "TPO","ORN","SAR","MVA","6V1","M3L","MLZ","MEA",
                      "IAS","CAF","YCM","CYS")
      for(pos in 1:length(three.list)){
        pdb_main[AA==three.list[pos],AA:=aa.list[pos]]
      }
      #get rid of alpha carbons and amine nitrogens and oxygens
      pdb_main <- pdb_main[(!(TypeB%chin%c("C","O","N","CA"))|
                            (!(TypeB%chin%c("C","O","N"))&AA=="G"))
                         ]
      #get rid of duplicated atoms, get mean position per residues and
      #keep only Beta carbon as representative atom per residue
      pdb_main <- pdb_main[!duplicated(pdb_main[,c("TypeA","TypeB",
                                                 "Chain","Position")])
                         ][,X:=mean(X),by=c("AA","Chain","Position")
                           ][,Y:=mean(Y),by=c("AA","Chain","Position")
                             ][,Z:=mean(Z),by=c("AA","Chain","Position")
                               ][(TypeB=="CA" & AA=="G")|(TypeB=="CB")]

    pdb_file <- pdb_main[Chain==chainz] 
    # xyz <- as.data.table(t(colMeans(pdb_file[,c("X","Y","Z")])))
    
    #Counter for non-standard residues not in list
    let3 <- pdb_file[nchar(AA)>1,AA]
    #Don't include nonstandard amino acids in single letter code
    protein <- paste(pdb_file[nchar(AA)==1,AA],collapse = "")
    indies <- unlist(stri_locate_all_fixed(protein,full_seq))
    if(is.na(indies[1])){
      indies <- unlist(stri_locate_all_regex(protein,paste0(B_start,".{0,5}",B_end)))
    }
    if(length(indies)>2){
      indies <- c(indies[1],indies[(length(indies)/2+1)])
      dt[ , touching_u:="2 REG MATCHES"]
      return(dt)
    }
    
    if(is.na(indies[1])){
      dt[ , touching_u:="0 REG MATCHES"]
      return(dt)
    } else{
      protein <- substr(protein,indies[1],indies[2]) # Beta hairpin sequence
      pdb_file <- pdb_file[nchar(AA)==1]
      pdb_file <- pdb_file[c(indies[1]:indies[2])] # Subset to beta hairpin residues
      n_B_end <- nchar(B_end) # Length of C-term beta sheet
      B1 <- nchar(B_start) # Length of N-term beta sheet
      l_prot <- nchar(protein) # Length of beta hairpin
      B2 <- (l_prot-n_B_end)+1 # Start of C-term beta sheet
      Turn_start <- B1 + 1
      Turn_end <- B2 - 1
      turn_pdb <- pdb_file[Turn_start:Turn_end]
      dt[, PHI := list(unlist(turn_pdb$PHI))
        ][ , PSI := list(unlist(turn_pdb$PSI))
        ][ , ALPHA := list(unlist(turn_pdb$ALPHA))
        ][ , KAPPA := list(unlist(turn_pdb$KAPPA))
        ][ , TCO := list(unlist(turn_pdb$TCO))]
      # get mean accessibility for odd and even residues from start
      for(v in 1:2){
        
        
        if(v==1){ #set start and end index
          BBB <- 1
          startsies <- B1
          by_n = -2
          start_diff = -1
        }else{
          startsies <- B2
          BBB <- l_prot
          by_n = 2
          start_diff = 1
        }
        # get mean accessibility for different conditions
        # Mean accessibility for full beta sheet
        mean_acc <- 
          mean(c(pdb_file[seq(startsies,BBB,by_n)
                        ][!is.na(acc),acc],
                 pdb_file[seq(from = (startsies + start_diff),
                              to = BBB,
                              by = by_n)
                        ][!is.na(acc), acc]))
        # Mean accessibility for face 1 of beta sheet
        mean_acc1 <- 
          mean(pdb_file[seq(startsies,BBB,by_n)
                      ][!is.na(acc),acc])
        # Mean accessibility for face 2 of beta sheet
        mean_acc2 <- 
          mean(pdb_file[seq(from = (startsies + start_diff),
                            to = BBB,
                            by = by_n)
                      ][!is.na(acc),acc])
        
        # Logic to test if odd or even residues are less accessible, 
        # or both are equal. Assumption is that beta hairpins have
        # a polar (accessible) and nonpolar (less accessible/buried) face.
        # 2 indicates TURN - Polar - Hyd - Polar order, 1 is the opposite.
        if(mean_acc1 < mean_acc2){ 
          
          by1 <- 1
          low_acc <- mean_acc1
        } else if(mean_acc2 < mean_acc1) {
          by1 <- 2
          low_acc <- mean_acc2
        } else { 
          by1 <- 3
          low_acc <- mean_acc2
        }
        # Write accessibility info to dt
        if(v==1){ 
          dt[ , B1_inner:=by1][ , B1_mean_acc:=mean_acc][ , B1_low_acc:=low_acc]
        } else{
          dt[ , B2_inner:=by1][ , B2_mean_acc:=mean_acc][ , B2_low_acc:=low_acc]
        }
      }
      #empty vectors for later
      distances <- c() #will hold distances
      dice <- c() #will hold amino acid contact info
      cutoffs <- c() #holds distance info for contacting residues
      sites <- c() #holds info for contactict residues, temporarily
      B1_pos <- c() #holds info for Beta sheet 1 residues with contacting partners
      B2_pos <- c() #as above but for Beta sheet 2
      
      while (B1 > 0){ #keep looping to find contacts until reach last residue
        #get temp pdb_file of one B1 residue and all B2 residues for contact analysis
        temp_pdb <- pdb_file[c(B1,B2:l_prot)]
        #get distance between single B1 residue and all B2 residues
        temp_dist <- dist(temp_pdb[,c("X","Y","Z")])
        temp_dist <- temp_dist[1:nrow(temp_pdb)-1]
        minis <- which(temp_dist<=8) #store contacting residues
        distances <- c(distances,temp_dist) #append distances
        B1 <- B1-1 #increment B1
        if(length(minis) > 0){ # if residues are contacting, update lists
          cutoffs<-c(cutoffs, temp_dist[minis])
          sites<-c(sites, paste0(temp_pdb[1,AA],
                                temp_pdb[(minis+1),AA]))
          B1_pos <- c(B1_pos,
                    c(rep((B1+1), length(minis))))
          B2_pos <- c(B2_pos,c(minis))
        }
        
      }
      
      if(length(sites)>0){ #if residues contact, keep info. Otherwise add placeholder
        dice <- sites
      } else {
        dice <- "-"
        cutoffs <- "-"
        B1_pos <- "-"
        B2_pos <- "-"
      }
      #update dt
      dt[ , mean_dist:=mean(distances)
       ][ , sd_dist:=sd(distances)
       ][ , n_contacts:=length(sites)
       ]
      
      if(length(dice)>=4){ #if 4 or more residues contact, consider motif valid
        dt[ , touching_u:="YES"]
      }
        dt_out <- cbind(dt,
                      data.table(pairs=dice,
                                 distances=cutoffs,
                                 B1_pos=B1_pos,
                                 B2_pos=B2_pos))
    }
  return(dt_out)
}


#function to count occurrences of peptides
#newsies,beta,pos
aa_loc <- function(peptides.in,param2,param3){
  #lists to convert AA to full length names
  aa.list <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  three.list <- c("Alanine","Cysteine","Aspartic Acid","Glutamic Acid","Phenylalanine","Glycine","Histidine","Isoleucine","Lysine","Leucine","Methionine","Asparagine","Proline","Glutamine","Arginine","Serine","Threonine","Valine","Tryptophan","Tyrosine")
  if(nchar(param3)>1){
    numerals<-"NO"
    param3_n=as.numeric(substr(param3,2,2))
  } else{
    numerals<-"YES"
    param3_n<-as.numeric(param3)
  }
  #get the number of characters per residue
  peptides.in[,n_char:=nchar(eval(param2))]
  #output factor as positiontype_position (e.g. B1_O1,Turn_M1)
  peptides.in[,factor_var:=paste(c(param2,param3),collapse="_"),by=n_char]
  peptides.in[,dipept:="NA"]
  #separate by distance from each end
  if(numerals=="YES"){
    #peptides.in[(n_char)>=param3_n,dipept:=substr(eval(param2),(n_char-param3_n)+1,(n_char-param3_n)+1),by=ID]
    peptides.in[(n_char)>=param3_n,dipept:=substr(eval(param2),param3_n,param3_n),by=ID]
  } else{
    if(substr(param3,1,1)=="O"){
      peptides.in[(n_char/2)>=param3_n,dipept:=substr(eval(param2),param3_n,param3_n),by=ID]
    } else if(substr(param3,1,1)=="H"){
      peptides.in[(n_char/2)>=param3_n,dipept:=substr(eval(param2),(n_char+1)-param3_n,(n_char+1)-param3_n),by=ID]
    } else if(substr(param3,1,1)=="M"){
      if(param3_n==1){
        peptides.in[(n_char/2)==2.5,dipept:=substr(eval(param2),3,3),by=ID]
      }else if(param3_n==2){
        peptides.in[n_char/2==3.5,dipept:=substr(eval(param2),4,4),by=ID]
      }else if(param3_n==3){
        peptides.in[(n_char/2)==1.5,dipept:=substr(eval(param2),2,2),by=ID]
      }
      
    }
  }
  
  peptides.in<-peptides.in[dipept!="NA"]
  for(aa in 1:length(aa.list)){
    peptides.in[dipept==aa.list[aa],dipept:=three.list[aa],by=dipept]
  } 
  
  data.table(peptides.in %>% 
               select(ID,factor_var,n_char,touching_u, dipept,B1_inner,B2_inner)) -> 
    peptides.in
  return(peptides.in)
}

#function to count occurrences of amino acids
aa_count<-function(peptides.in,param2,param3){
  aa.list <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  three.list <- c("Alanine","Cysteine","Aspartic Acid","Glutamic Acid","Phenylalanine","Glycine","Histidine","Isoleucine","Lysine","Leucine","Methionine","Asparagine","Proline","Glutamine","Arginine","Serine","Threonine","Valine","Tryptophan","Tyrosine")
  peptides.in[,factor_var:=param3]
  for(aa in 1:length(aa.list)){
    peptides.in[,three.list[aa]:=stri_count(eval(param2),fixed=aa.list[aa])]

  }
  data.table(peptides.in %>% select(id,factor_var, Alanine:Tyrosine, beta_type)) -> peptides.in
  return(peptides.in)
}

#function to count occurence of dipeptides
dipept_count<-function(dt){
  #make sure there are no NA pairs
  dt <- dt[pairs!="-"]
  #make vector of all permutations for amino acid pairs
  aa.list <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  pbj <- permutations(20,2,aa.list,repeats.allowed = TRUE)
  pbj <- paste0(pbj[,1],pbj[,2])
  
  #make table of counts for each pair, 
  #normalized as observed frequency,
  #e.g. occurrence of pair at position/all pairs at position
  pairs <- dt[,.(pairs, beta_type)]
  pairs <- as.data.frame.matrix(t(table(pairs)))
  pairs$id <- rownames(pairs)
  runit <- colnames(pairs)
  
  # print(pairs)
  # print(rownames(pairs))
  data.table(pairs %>% 
               gather(AA,obs_freq,AA:YY))-> pairs
  pairs[, obs_freq := as.numeric(obs_freq)]
  
  pairs[,obs_freq:=as.numeric(obs_freq/sum(obs_freq)), by = id]
  #get total number of amino acids in pairs for normalization below
  all_aa <- sum(stri_count_regex(dt[,pairs],"[:alpha:]"))
  dt[, n_pairs := .N, by = beta_type]
  all_aa_table <- unique(dt[,.(beta_type,n_pairs)])
  all_aa_table[, id := beta_type
             ][, beta_type := NULL]
  #empty vector for ratios
  rats<-c()
  
  #loop through AA pairs for the Full condition to get 
  #expected frequencies
  for(runs in runit){
    temp_rat<-(sum(stri_count_fixed(
      dt[,pairs],substring(runs,1,1)))/all_aa)*
      (sum(stri_count_fixed(
        dt[,pairs],substring(runs,2,2)))/all_aa)
    rats<-c(rats,temp_rat)
    
  }
  
  #pairs is dt to hold expected counts for Full condition
  pairs <- merge(pairs, data.table(AA = runit, 
                                   exp_freq = rats))
  # print(pairs)
  # print(all_aa_table)
  pairs <- merge(pairs, all_aa_table, by = "id")
  pairs[obs_freq>0,ratio:=log2(obs_freq/exp_freq)
        ][obs_freq<=0,ratio:=0
          ][,factor_var:="Full"]
  #append DT's and return
  # pairs <- rbind(pairs)
  return(pairs)
}

ss_fasta2dt <- function(fasta_in, ss_names){
  
  # length(unlist(getSequence(
  #   fasta_in[str_replace(ss_names,":secstr",":sequence")],
  #   as.string = T))) %>%
  # length(unlist(getSequence(
  #   fasta_in[ss_names],
  #   as.string = T))) %>%
  data.table(
    "id" = str_remove(ss_names,":secstr"),
    "pdb_id" = str_remove(ss_names,":.*:.*"),
    "chain_id" = str_remove(str_remove(ss_names,":secstr"), ".*:"),
    "aa_seq" = unlist(
      getSequence(
        fasta_in[str_replace(ss_names,":secstr",":sequence")],
        as.string = T)),
    "sec_str" = unlist(
      getSequence(
        fasta_in[ss_names],
        as.string = T)))
  
}

filter_secondary_structure <- function(ss_dt,
                                       reg_pattern,
                                       use_regex = TRUE){
  as.data.table(
    ss_dt)[
      grepl(reg_pattern,
            sec_str,
            fixed = !(use_regex))]
}

dl_pdbs <- function(pdb, folder_out){
  wget_code <- paste0("wget -P ",
                      folder_out,
                      " https://files.rcsb.org/download/",
                      pdb,
                      ".cif.gz")
  system(wget_code)
}


sec_str_location <- 
  function(row_pos,
           dt,
           secstr,
           amino_acid,
           reg_pattern,
           col_prefix,
           use_regex = TRUE,
           trim_end = 0,
           trim_start = 0,
           index_col = TRUE,
           turn = FALSE){

    location_mat <-
      str_locate_all(
        as.character(dt[row_pos,get(secstr)]),
        reg_pattern)
    loc_start<- as.integer(
      unlist(
        location_mat[[1]][,1]
        )
      ) + trim_start
    loc_end<- as.integer(
      unlist(
        location_mat[[1]][,2]
        )
      ) - trim_end
    secondary <- c()
    primary <- c()
    
    for(t in 1:(length(loc_start))){
      secondary <- c(secondary,
                     substr(dt[row_pos, get(secstr)],
                            loc_start[t],
                            loc_end[t]))
      primary <- c(primary,
                   substr(dt[row_pos, get(amino_acid)],
                          loc_start[t],
                          loc_end[t]))
    }
    cns <- c(colnames(dt),
             paste0(col_prefix,
                  c("_aa","_ss")))
    dt2 <- 
      data.table(
        dt[row_pos],
        w = primary,
        x = secondary
        )
    if(isTRUE(turn)){
      n_1 <- substr(dt[row_pos, get(amino_acid)],1,(loc_start-1))
      n_1_rev <- paste(rev(unlist(strsplit(n_1,""))),collapse="")

      n_2 <- substr(dt[row_pos, get(amino_acid)],(loc_end+1),nchar(dt[row_pos, get(amino_acid)]))
      n_2_rev <- paste(rev(unlist(strsplit(n_2,""))),collapse="")
      dt2 <- data.table(
        dt2,
        n_1 = n_1,
        n_2 = n_2,
        n_1_rev = n_1_rev,
        n_2_rev = n_2_rev
      )
      cns <- c(cns,
               paste0(col_prefix,
                    c("_Beta1","_Beta2",
                      "_Beta1_rev","_Beta2_rev")))
    }
    if(isTRUE(index_col)){
      dt2 <- data.table(
        dt2,
        y = loc_start,
        z = loc_end
      )
      cns <- c(cns,
               paste0(col_prefix,
                      c("_start","_stop")))
      
    }
    colnames(dt2) <- cns
    return(dt2)
    
  }


aa_breakdown <- function(hairpins_dt){
  hydro_aa <- c("A","V","I","L","M")
  bulk_aa <- c("F","Y","W")
  special_aa <- c("C","G","P")
  polar_aa <- c("S","T","N","Q")
  pos_aa <- c("R","H","K")
  neg_aa <- c("D","E")
  hairpins_dt %>%
    filter(!grepl("3",beta_type)) %>%
    mutate(B1_correct = stri_reverse(subseq_Beta1),
           turn_length = nchar(subseq_aa),
           turn_type = paste0("turn_", turn_length, 
                              "_", beta_type)) %>%
    separate(B1_correct,
             paste0("B1_",c(1:14)),
             remove = FALSE,
             sep = c(1:13)
    ) %>%
    separate(subseq_Beta2,
             paste0("B2_",c(1:14)),
             remove = FALSE,
             sep = c(1:13)
    ) %>%
    separate(subseq_aa,
             paste0("turn_",c(1:5)),
             remove = FALSE,
             sep = c(1:5)
    ) %>%
    pivot_longer(c(paste0("B1_",c(1:14)),
                   paste0("B2_",c(1:14)),
                   paste0("turn_",c(1:5))),
                 names_to = "structure",
                 values_to = "amino_acid") %>%
    filter(amino_acid != "") %>%
    mutate(B1_inner = if_else(B1_inner == 1,
                              "Pol",
                              "Hydro"),
           B2_inner = if_else(B2_inner == 1,
                              "Pol",
                              "Hydro")) %>%
    mutate(B1_inner = paste0("N-term ",B1_inner),
           B2_inner = paste0("C-term ",B2_inner),
           turn_length = paste0("T ",turn_length)) %>%
    pivot_longer(c(turn_length, B1_inner, B2_inner),
                 names_to = "sub_thing",
                 values_to = "sub_pos") %>%
    filter(str_remove(structure, "_.*") ==
             str_remove(sub_thing, "_.*")) %>%
    tabyl(structure, amino_acid, sub_pos,
          show_missing_levels = FALSE)%>%
    adorn_totals("col", name = "TOTES") %>%
    adorn_percentages("row") %>%
    adorn_pct_formatting(digits = 3,
                         affix_sign = FALSE) %>%
    rbindlist(idcol = TRUE, fill = TRUE) %>%
    select(-X) %>%
    pivot_longer(cols = c("A":"Y"),
                 names_to = "AA",
                 values_to = "Perc") %>%
    mutate(color_code = str_remove(structure,"_[1-9]|_[1-9][1-9]"),
           turn_type = .id,
           .id = NULL,
           Perc = as.numeric(Perc)) %>%
    mutate(AA_group = case_when(
      AA %in% hydro_aa ~ "Hydro",
      AA %in% special_aa ~ "Special",
      AA %in% bulk_aa ~ "Bulky",
      AA %in% polar_aa ~ "Polar",
      AA %in% neg_aa ~ "Negative",
      TRUE ~ "Positive"
    ),
    AA_simpler = if_else(
      AA_group %in% c("Negative","Positive"),
      "Charged",
      AA_group
    ),
    AA_simpler = factor(AA_simpler,
                        levels = c(
                          "Bulky",
                          "Hydro",
                          "Polar",
                          "Charged",
                          "Special"
                        )),
    AA = fct_relevel(AA, names(color_pal)),
    structure = str_replace(structure,"B1_","N"),
    structure = str_replace(structure,"B2_","C"),
    structure = str_replace(structure,"turn_",""),
    structure = fct_relevel(structure, hairpin_levels),
    turn_type = fct_relevel(turn_type,
                            c("N-term Pol",
                              "N-term Hydro",
                              "T 1",
                              "T 2",
                              "T 3",
                              "T 4",
                              "T 5",
                              "C-term Pol",
                              "C-term Hydro")))
}

# Get counts for residues at specific turn positions for given subtypes
turn_breakdown <- function(hairpins_dt){
  hydro_aa <- c("A","V","I","L","M")
  bulk_aa <- c("F","Y","W")
  special_aa <- c("C","G","P")
  polar_aa <- c("S","T","N","Q")
  pos_aa <- c("R","H","K")
  neg_aa <- c("D","E")
  aa_list <- c("A","C","D","E","F",
               "G","H","I","K","L",
               "M","N","P","Q","R",
               "S","T","V","W","Y")
  hairpins_dt %>%
    filter(!grepl("3",beta_type)) %>%
    mutate(turn_length = nchar(subseq_aa),
           turn_size = paste0("turn_", turn_length, "_", turn_type)) %>%
    separate(subseq_aa,
             paste0("turn_",c(1:5)),
             remove = FALSE,
             sep = c(1:5)
    ) %>%
    pivot_longer(c(paste0("turn_",c(1:5))),
                 names_to = "structure",
                 values_to = "amino_acid") %>%
  filter(amino_acid != "") %>%
  mutate(turn_length = paste0("T ",turn_length)) %>%
  tabyl(structure, amino_acid, turn_size,
        show_missing_levels = FALSE) %>%
    adorn_totals("col", name = "TOTES") %>%
    rbindlist(idcol = TRUE, fill = TRUE) %>%
    pivot_longer(cols = all_of(c(aa_list)),
                 names_to = "AA",
                 values_to = "Perc") %>%
    mutate(color_code = str_remove(structure,"_[1-9]|_[1-9][1-9]"),
           turn_type = .id,
           .id = NULL,
           Perc = as.numeric(Perc),
           Perc = if_else(is.na(Perc),
                          0, Perc/TOTES)) %>%
    mutate(AA_group = case_when(
      AA %in% hydro_aa ~ "Hydro",
      AA %in% special_aa ~ "Special",
      AA %in% bulk_aa ~ "Bulky",
      AA %in% polar_aa ~ "Polar",
      AA %in% neg_aa ~ "Negative",
      TRUE ~ "Positive"
    ),
    AA_simpler = if_else(
      AA_group %in% c("Negative","Positive"),
      "Charged",
      AA_group
    ),
    AA_simpler = factor(AA_simpler,
                        levels = c(
                          "Bulky",
                          "Hydro",
                          "Polar",
                          "Charged",
                          "Special"
                        )),
    structure = str_remove(structure,"turn_"),
    turn_subtype = str_remove(turn_type, "turn_._"),
    full_turn_info = turn_type,
    turn_type = str_match(turn_type, "turn_."),
    AA = fct_relevel(AA, names(color_pal)))
 
}

# This function takes a formatted DT with PHI and PSI angles for residues in a turn
# region and classifies the type of turn based on existing nomenclature. Assumes 1 residue
# turns are actually gamma, 2 residue turns are Beta, and 3 residue turns are alpha turns 
# respectively.
classify_turns <- function(turn_dt){
  
  turn_number <- c(1:5)
  turn_names_1 <- c("Gamma_normal","Gamma_inverse")
  turn_names_2 <- c("Beta_I", "Beta_I'", "Beta_II", "Beta_II`", "Beta_VIa1", "Beta_VIa2", "Beta_VIb", "Beta_VIII")
  turn_names_4 <- c("Beta_I", "Beta_I'", "Beta_II", "Beta_II`", "Beta_VIa1", "Beta_VIa2", "Beta_VIb", "Beta_VIII")
  turn_names_5 <- c("Alpha_I_RS",
                    "Alpha_I_LS",
                    "Alpha_II_RS",
                    "Alpha_II_LS",
                    "Alpha_I_RU",
                    "Alpha_I_LU",
                    "Alpha_II_RU",
                    "Alpha_II_LU",
                    "Alpha_I_C")
  turn_names_3 <- turn_names_5
  turn_names <- list(turn_names_1,
                     turn_names_2,
                     turn_names_3,
                     turn_names_4,
                     turn_names_5)
  
  turn_dt <-  turn_dt %>%
    mutate(turn_type = 
             case_when(
               turn_length == 1 ~ "Gamma_other",
               turn_length == 2 ~ "Beta_IV",
               turn_length == 3 ~ "Alpha_Other",
               turn_length == 4 ~ "Beta_IV",
               turn_length == 5 ~ "Alpha_Other"))
  for( tn in  turn_number ){
    
    for( turn_classifier_position in 1:length(turn_names[[tn]])){
      name_out <- turn_names[[tn]][turn_classifier_position]
      
      # BETA
      if(tn == 4 | tn == 2 ){
        PHI_4_2_list <- c(-60, 60, -60, 60, -60, -120, -135, -60)[turn_classifier_position]
        PSI_4_2_list <- c(-30, 30, 120, -120, 120, 120, 135, -30)[turn_classifier_position]
        PHI_4_3_list <- c(-90, 90, 80, -80, -90, -60, -75, -120)[turn_classifier_position]
        PSI_4_3_list <- c(0, 0, 0, 0, 0, 0, 160, 120)[turn_classifier_position]
        
        PHI_4_2_low <- min(PHI_4_2_list-45, PHI_4_2_list+45)
        PHI_4_2_high <- max(PHI_4_2_list-45, PHI_4_2_list+45)
        PSI_4_2_low <- min(PSI_4_2_list-45, PSI_4_2_list+45)
        PSI_4_2_high <- max(PSI_4_2_list-45, PSI_4_2_list+45)
        PHI_4_3_low <- min(PHI_4_3_list-45, PHI_4_3_list+45)
        PHI_4_3_high <- max(PHI_4_3_list-45, PHI_4_3_list+45)
        PSI_4_3_low <- min(PSI_4_3_list-45, PSI_4_3_list+45)
        PSI_4_3_high <- max(PSI_4_3_list-45, PSI_4_3_list+45)
        if(tn == 4){
          turn_dt <- turn_dt %>%
            mutate(turn_type = if_else(
              turn_length == tn &
                PHI_2 >= PHI_4_2_low & PHI_2 <= PHI_4_2_high &
                PSI_2 >= PSI_4_2_low & PSI_2 <= PSI_4_2_high & 
                PHI_3 >= PHI_4_3_low & PHI_3 <= PHI_4_3_high &
                PSI_3 >= PSI_4_3_low & PSI_3 <= PSI_4_3_high,
              name_out,
              turn_type))
        } else {
          turn_dt <- turn_dt %>%
            mutate(turn_type = if_else(
              turn_length == tn &
                PHI_1 >= PHI_4_2_low & PHI_1 <= PHI_4_2_high &
                PSI_1 >= PSI_4_2_low & PSI_1 <= PSI_4_2_high & 
                PHI_2 >= PHI_4_3_low & PHI_2 <= PHI_4_3_high &
                PSI_2 >= PSI_4_3_low & PSI_2 <= PSI_4_3_high,
              name_out,
              turn_type))
        }
      # GAMMA  
      } else if (tn == 1){
        PHI_4_2_list <- c(75, -79)[turn_classifier_position]
        PSI_4_2_list <- c(-64, 69)[turn_classifier_position]
        
        PHI_4_2_low <- min(PHI_4_2_list-45, PHI_4_2_list+45)
        PHI_4_2_high <- max(PHI_4_2_list-45, PHI_4_2_list+45)
        PSI_4_2_low <- min(PSI_4_2_list-45, PSI_4_2_list+45)
        PSI_4_2_high <- max(PSI_4_2_list-45, PSI_4_2_list+45)
        turn_dt <- turn_dt %>%
          mutate(turn_type = if_else(
            turn_length == tn &
              PHI_1 >= PHI_4_2_low & PHI_1 <= PHI_4_2_high &
              PSI_1 >= PSI_4_2_low & PSI_1 <= PSI_4_2_high,
            name_out,
            turn_type))
      # ALPHA  
      } else if (tn == 5 | tn == 3) {
        PHI_2_list <- c(-60, 48, -69, 53, 59, -61, 54, -65, -103)[turn_classifier_position]
        PSI_2_list <- c(-29, 42, 129, -137, -157, 158, 39, -20, 143)[turn_classifier_position]
        PHI_3_list <- c(-72, 67, 88, -95, -67, 64, 67, -90, -85)[turn_classifier_position]
        PSI_3_list <- c(-29, 33, -16, 81, -29, 37, -5, 16, 2)[turn_classifier_position]
        PHI_4_list <- c(-96, 70, -91, 57, -68, 62, -125, 86, -54)[turn_classifier_position]
        PSI_4_list <- c(-20, 32, -32, 38, -39, 39, -34, 37, -39)[turn_classifier_position]
        
        PHI_2_low <- min(PHI_2_list-45, PHI_2_list+45)
        PHI_2_high <- max(PHI_2_list-45, PHI_2_list+45)
        PSI_2_low <- min(PSI_2_list-45, PSI_2_list+45)
        PSI_2_high <- max(PSI_2_list-45, PSI_2_list+45)
        PHI_3_low <- min(PHI_3_list-45, PHI_3_list+45)
        PHI_3_high <- max(PHI_3_list-45, PHI_3_list+45)
        PSI_3_low <- min(PSI_3_list-45, PSI_3_list+45)
        PSI_3_high <- max(PSI_3_list-45, PSI_3_list+45)
        PHI_4_low <- min(PHI_4_list-45, PHI_4_list+45)
        PHI_4_high <- max(PHI_4_list-45, PHI_4_list+45)
        PSI_4_low <- min(PSI_4_list-45, PSI_4_list+45)
        PSI_4_high <- max(PSI_4_list-45, PSI_4_list+45)
        
        if(tn == 5){
          turn_dt <- turn_dt %>%
            mutate(turn_type = if_else(
              turn_length == tn &
                PHI_2 >= PHI_2_low & PHI_2 <= PHI_2_high &
                PSI_2 >= PSI_2_low & PSI_2 <= PSI_2_high & 
                PHI_3 >= PHI_3_low & PHI_3 <= PHI_3_high &
                PSI_3 >= PSI_3_low & PSI_3 <= PSI_3_high &
                PHI_4 >= PHI_4_low & PHI_4 <= PHI_4_high &
                PSI_4 >= PSI_4_low & PSI_4 <= PSI_4_high,
              name_out,
              turn_type))
        } else {
          turn_dt <- turn_dt %>%
            mutate(turn_type = if_else(
              turn_length == tn &
                PHI_1 >= PHI_2_low & PHI_1 <= PHI_2_high &
                PSI_1 >= PSI_2_low & PSI_1 <= PSI_2_high & 
                PHI_2 >= PHI_3_low & PHI_2 <= PHI_3_high &
                PSI_2 >= PSI_3_low & PSI_2 <= PSI_3_high &
                PHI_3 >= PHI_4_low & PHI_3 <= PHI_4_high &
                PSI_3 >= PSI_4_low & PSI_3 <= PSI_4_high,
              name_out,
              turn_type))
        }
        
        
      }
      
    }
  }
  return(turn_dt)
}

