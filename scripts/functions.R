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
      #CHECK HERE!!!!!!!!!!!!
      pdb_main <- gsub("\\?[[:space:]]+\\?[[:space:]]+\\?[[:space:]]+\\?[[:space:]]+\\?[[:space:]]+\\?","?", pdb_main, fixed=FALSE)
      pdb_main <- unique(fread(pdb_main)[,!c("V1","V2","V7","V14","V15",
                                           "V16","V18","V20")])
      
      # if(length(colnames(pdb_main))>15){
      #   # print(pdb_ID)
      #   pdb_main <- pdb_main[,!c("V17","V19","V21","V23","V25")]
      #   # print(pdb_main)
      # }
      
      
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
          }, warning = function(w) {
            dt[ , touching_u:="NO DSSP"]
            return(dt)
          }, error = function(e) {
            dt[ , touching_u:="NO DSSP"]
            return(dt)
          }, finally = {
          })
      if(is.null(access$acc)){
        return(dt)
      }
 
      pdb_main <- merge(pdb_main,access[,c("Position","acc","Chain")],by=c("Position","Chain"))
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
    xyz <- as.data.table(t(colMeans(pdb_file[,c("X","Y","Z")])))
    
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

