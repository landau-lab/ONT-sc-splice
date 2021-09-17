## calculate log odds ratio for all clusters with >=2 junctions per cluster and total reads in genotyped cells >= 5
## permute by shuffling genotype  

library(Matrix)
library(parallel)
library(tidyverse)
library(matrixStats)

args = commandArgs(TRUE)
path.to.matrix = args[1]
path.to.genotype = args[2]
path.to.metadata = args[3]
nperm = args[4]
output.dir = args[5]

nperm = as.numeric(nperm)

# Test data
#path.to.matrix = "/gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/MDS_combined/1.Counts_matrix"
#path.to.genotype = "/gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/MDS_combined/2.Genotype_info"
#path.to.metadata = "/gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/MDS_combined/4.Annotated_metadata"


runJunctionPermutation = function(path.to.matrix, path.to.genotype, path.to.metadata, nperm){
  patient.names = c("MDS_P1", "MDS_P2", "MDS_P3")
  ## load metadata as list of files 
  setwd(path.to.metadata)
  metafiles = list.files(path.to.metadata)
  metadata.list = lapply(metafiles, function(x) read.csv(file = x))
  names(metadata.list)= patient.names
  for (patient in patient.names) {
    metadata.list[[patient]]$intron_junction = paste(metadata.list[[patient]]$intron_junction, metadata.list[[patient]]$clusterID, sep = ":")
    rownames(metadata.list[[patient]]) = metadata.list[[patient]]$intron_junction
  }
  
  #merge by junctions in both metadata and counts matrix
  #read in matrix as named list of matrix 
  setwd(path.to.matrix)
  files = list.files(path.to.matrix)
  mtx.list = lapply(files, function(x) read.table(file=x))
  names(mtx.list) = patient.names
  for (patient in patient.names){
    shared = intersect(rownames(mtx.list[[patient]]), rownames(metadata.list[[patient]]))
    mtx.list[[patient]] = mtx.list[[patient]][shared,-1]
    metadata.list[[patient]] = metadata.list[[patient]][shared,]
  }
  message("Matrix loaded")
  
  #load genotype information
  genotype.list = list()
  cell.type.list = list()
  
  setwd(path.to.genotype)
  pattern = c("_1", "_2", "_3")
  i = 1
  for (file in list.files(path.to.genotype)){
    geno = as.data.frame(read.table(file))
    geno = geno[which(geno$Cell.Assignment != "NA"),]
    geno = geno[,c("Genotype_1UMI", "Cell.Assignment")]
    rownames(geno) = gsub(rownames(geno), pattern = pattern[i], replacement = "")
    geno = geno[rownames(geno) %in% colnames(mtx.list[[i]]),]
    cell.type.list[[patient.names[i]]] = geno$Cell.Assignment
    geno$Genotype_1UMI = paste(geno$Genotype_1UMI, geno$Cell.Assignment, sep= "_")
    genotype = geno$Genotype_1UMI
    genotype = as.character(genotype)
    names(genotype) = rownames(geno)
    genotype.list[[patient.names[i]]] = genotype
    i = i+1
  }
  message("Genotype loaded")
  
  patient.OR = list()
  cell.totals = list()
  mtx.filt = list()
  metadata.filt = list()
  data.filt.list = list()
  clust.list = list()
  cell.type.update = list()
  
  for (patient in patient.names){
    
    patient.OR[[patient]] = list()
    cell.totals[[patient]] = list()
    mtx.filt[[patient]] = list()
    metadata.filt[[patient]] = list()
    data.filt.list[[patient]] = list()
    clust.list[[patient]] = list()
    cell.type.update[[patient]] = list()
    
    for (cell.type in unique(as.character(cell.type.list[[patient]]))){
      
      wt.geno = paste("WT", cell.type, sep = "_")
      mut.geno = paste("MUT", cell.type, sep = "_")
      genotype = genotype.list[[patient]][genotype.list[[patient]] %in% c(wt.geno, mut.geno)]
    
      wt = names(genotype.list[[patient]])[genotype.list[[patient]] == wt.geno]
      mut = names(genotype.list[[patient]])[genotype.list[[patient]] == mut.geno]
      
      if (length(wt) > 2 & length(mut) > 2) {
        
        #filter data for reads >= 5, clusters with junctions > 1 based on total observations 
        data = data.frame(intron_junction = metadata.list[[patient]]$intron_junction, clusterID = metadata.list[[patient]]$clusterID)
        data$obs.wt = rowSums(mtx.list[[patient]][,colnames(mtx.list[[patient]]) %in% wt])
        data$obs.mut = rowSums(mtx.list[[patient]][,colnames(mtx.list[[patient]]) %in% mut])
        data$total.reads.per.junction = data$obs.wt+data$obs.mut
        data = data[which(data$total.reads.per.junction >= 5),]
        
        #get list of unique clusters with n > 1
        clusters = data %>% group_by(clusterID) %>% tally()
        clust_list = clusters %>% filter(n != 1) %>% select(clusterID)
        clust.list[[patient]][[cell.type]] = clust_list
        
        if (length(clust_list$clusterID) > 2){
          
          cell.type.update[[patient]][[cell.type]] = cell.type
        
          #filter mtx and metadata for only junctions with read count >= 5
          data.filt = data[which(data$clusterID %in% clust_list$clusterID),]
          covered_junc = data.filt$intron_junction
          mtx.filt[[patient]][[cell.type]] = mtx.list[[patient]][covered_junc,]
          metadata.filt[[patient]][[cell.type]] = metadata.list[[patient]][covered_junc,]
          data.filt.list[[patient]][[cell.type]] = data.filt
          message("Data filtered")
          
          obs.ratio = list()
          
          #get log odds ratio for each cluster
          for(clust in clust_list$clusterID){
            subset = data.filt[which(data.filt$clusterID == clust),]
            subset$mut.all = sum(subset$obs.mut)
            subset$wt.all = sum(subset$obs.wt)
            subset$mut.remain = subset$mut.all - subset$obs.mut
            subset$wt.remain = subset$wt.all - subset$obs.wt
            rownames(subset) = subset$intron_junction
            for (junc in rownames(subset)){
              obs.ratio[junc] = log((subset[junc,"obs.mut"]/subset[junc,"mut.remain"])/(subset[junc,"obs.wt"]/subset[junc,"wt.remain"]))
            }
          }
          #get total cells and obs ratio
          obs.ratio.num = as.numeric(obs.ratio)
          names(obs.ratio.num) = names(obs.ratio)
          total.cells = sum(length(wt)+length(mut))
          cell.totals[[patient]][[cell.type]] = total.cells
          
          #calculate total cluster coverage 
          cluster.cov = data.filt %>% dplyr::group_by(clusterID) %>% dplyr::summarise(mut.cluster.cov = sum(obs.mut), wt.cluster.cov = sum(obs.wt))
          data.filt = left_join(data.filt, cluster.cov, by = "clusterID")
          
          #create data frame of final output
          patient.output = data.frame(obs.logOR.ratio = obs.ratio.num, total.cells = rep(total.cells, length(obs.ratio.num)), intron_junction = names(obs.ratio))
          patient.output = left_join(patient.output, data.filt, by = "intron_junction")
          rownames(patient.output) = patient.output$intron_junction
          colnames(patient.output) = lapply(colnames(patient.output), function(x) paste(patient,cell.type,x,sep="_"))
          patient.OR[[patient]][[cell.type]] = patient.output
        }
      }
    }
  }
  

  #reduce one output from all patients into one dataframe from three dataframes 
  #first collapse by patient
  all = list()
  for (patient in patient.names) {
    all[[patient]] = patient.OR[[patient]] %>% map(~ .x %>%
                                                     rownames_to_column('intron_junction'))
    all[[patient]] = all[[patient]] %>% reduce(full_join, by = "intron_junction")
  }
  
  #next collapse into giant dataframe 
  for (patient in patient.names){
    all[[patient]]$intron_junction = sub(":clu.*", "", all[[patient]]$intron_junction)
  }
  all = all %>% reduce(full_join, by = "intron_junction") 
  
  all.cell.types.list = list()
  for (patient in patient.names){
    all.cell.types.list[[patient]] = names(cell.type.update[[patient]])
  }
  all.cell.types = unique(unlist(all.cell.types.list))
  
  for (patient in patient.names){
    for (cell.type in all.cell.types) {
      if (!(cell.type %in% all.cell.types.list[[patient]])) {
        all[,paste(patient, cell.type, "obs.logOR.ratio", sep = "_")] = NA
        all[,paste(patient, cell.type, "total.cells", sep = "_")] = 0
        cell.totals[[patient]][[cell.type]] = as.numeric(0)
      }
    }
  }
  
  
  
  #calculate mean observed difference for each cell type 
  weight.obs.OR = list()
  for (cell.type in all.cell.types){
    mat.names = lapply(patient.names, function(x) paste(x, cell.type, "obs.logOR.ratio", sep = "_"))
    obs.mat = all[,unlist(mat.names)]
    
    cell.totals.update = as.data.frame(cell.totals)
    cell.total.names = lapply(patient.names, function(x) paste(x,cell.type, sep = "."))
    cell.total.sub = cell.totals.update[, unlist(cell.total.names)]
    weights = as.numeric(cell.total.sub[1,]/rowSums(cell.total.sub[1,]))
    
    weight.obs.OR[[cell.type]] = rowWeightedMeans(as.matrix(obs.mat), weights)
    names(weight.obs.OR[[cell.type]]) = all$intron_junction
  }
  
  message("Observed difference calculated")
  
  #initialize to calculate whether or not obs OR > shf OR
  temp = list()
  for (cell.type in all.cell.types){
    temp[[cell.type]] = rep(0,length(weight.obs.OR[[1]][[1]]))
  } # Create an empty vector to update in each iteration
  
  #initiate shuffled data frame, do this only once 
  shf.data.list = list()
  for (patient in patient.names){
    shf.data.list[[patient]] = list()
    for (cell.type in unique(as.character(cell.type.list[[patient]]))){
      shf.data.list[[patient]][[cell.type]] = data.frame(intron_junction = data.filt.list[[patient]][[cell.type]]$intron_junction, clusterID = data.filt.list[[patient]][[cell.type]]$clusterID) 
    }
  }
  
  #run through for every permutation
  for(x in 0:nperm){
    
    set.seed(x)
    
    #initialize list for each patient to create list of results 
    shf.patient.OR= list()
    shf.cell.totals = list()
    shf.cell.types = list()
    
    #cycle through each patient and shuffle genotypes before calculating OR 
    #use only filtered data frames based on previous filtering done 
    for (patient in patient.names) {
      
      shf.patient.OR[[patient]] = list()
      shf.cell.totals[[patient]] = list()
      shf.cell.types[[patient]] = list()
      
      for (cell.type in unique(as.character(cell.type.list[[patient]]))){
        wt.geno = paste("WT", cell.type, sep = "_")
        mut.geno = paste("MUT", cell.type, sep = "_")
        genotype = genotype.list[[patient]][genotype.list[[patient]] %in% c(wt.geno, mut.geno)]
        genotype = genotype[names(genotype) %in% colnames(mtx.filt[[patient]][[cell.type]])]
        
        if (length(genotype > 2)) {
        
          orig.names = names(genotype)
          shuffle.genotype = sample(genotype, size = length(genotype), replace = F)
          names(shuffle.genotype) = orig.names
          
          wt = names(shuffle.genotype)[shuffle.genotype == wt.geno]
          mut = names(shuffle.genotype)[shuffle.genotype == mut.geno]
          
          if (length(wt) > 2 & length(mut) > 2) {
            
            shf.data = shf.data.list[[patient]][[cell.type]]
            
            shf.data$shf.wt = rowSums(mtx.filt[[patient]][[cell.type]][,colnames(mtx.filt[[patient]][[cell.type]]) %in% wt])
            shf.data$shf.mut = rowSums(mtx.filt[[patient]][[cell.type]][,colnames(mtx.filt[[patient]][[cell.type]]) %in% mut])
            shf.data$total.reads.per.junction = shf.data$shf.wt+shf.data$shf.mut
            
            shf.ratio = list()
            
            if (length(clust.list[[patient]][[cell.type]]$clusterID) > 2){
              
              shf.cell.types[[patient]][[cell.type]] = cell.type
              
              #log OR for each junction after shuffling genotypes within each patient 
              for(clust in clust.list[[patient]][[cell.type]]$clusterID){
                subset = shf.data[which(shf.data$clusterID == clust),]
                subset$mut.all = sum(subset$shf.mut)
                subset$wt.all = sum(subset$shf.wt)
                subset$mut.remain = subset$mut.all - subset$shf.mut
                subset$wt.remain = subset$wt.all - subset$shf.wt
                rownames(subset) = subset$intron_junction
                for (junc in subset$intron_junction){
                  shf.ratio[junc] = log((subset[junc,"shf.mut"]/subset[junc,"mut.remain"])/(subset[junc,"shf.wt"]/subset[junc,"wt.remain"]))
                }
              }
              #message("odds ratio calculated")
              
              shf.ratio.num = as.numeric(shf.ratio)
              names(shf.ratio.num) = names(shf.ratio)
              shf.total.cells = sum(length(wt)+length(mut))
              shf.cell.totals[[patient]][[cell.type]] = shf.total.cells
              
              #create output data frame with results for each patient 
              shf.patient.output = data.frame(shf.logOR.ratio = shf.ratio.num, intron_junction = names(shf.ratio))
              rownames(shf.patient.output) = shf.patient.output$intron_junction
              colnames(shf.patient.output) = lapply(colnames(shf.patient.output), function(x) paste(patient,cell.type,x,sep="_"))
              shf.patient.OR[[patient]][[cell.type]] = shf.patient.output
            }
          }
        }
      }
    }
   
    
    #reduce one output from all patients into one dataframe from three dataframes 
    #first collapse by patient
    shf.all = list()
    for (patient in patient.names) {
      shf.all[[patient]] = shf.patient.OR[[patient]] %>% map(~ .x %>%
                                                       rownames_to_column('intron_junction'))
      shf.all[[patient]] = shf.all[[patient]] %>% reduce(full_join, by = "intron_junction")
    }
    
    #next collapse into giant dataframe 
    for (patient in patient.names){
      shf.all[[patient]]$intron_junction = sub(":clu.*", "", shf.all[[patient]]$intron_junction)
    }
    shf.all = shf.all %>% reduce(full_join, by = "intron_junction") 
    
    shf.cell.types.list = list()
    for (patient in patient.names){
      shf.cell.types.list[[patient]] = names(shf.cell.types[[patient]])
    }
    shf.cells = unique(unlist(shf.cell.types.list))
    
    for (patient in patient.names){
      for (cell.type in shf.cells) {
        if (!(cell.type %in% shf.cell.types.list[[patient]])) {
          shf.all[,paste(patient, cell.type, "shf.logOR.ratio", sep = "_")] = NA
          shf.all[,paste(patient, cell.type, "total.cells", sep = "_")] = 0
          shf.cell.totals[[patient]][[cell.type]] = as.numeric(0)
        }
      }
    }       

    #calculate mean observed difference 
    weight.shf.OR = list()
    for (cell.type in all.cell.types){
      mat.names = lapply(patient.names, function(x) paste(x, cell.type, "shf.logOR.ratio", sep = "_"))
      shf.mat = shf.all[,unlist(mat.names)]
      
      shf.cell.totals.update = as.data.frame(shf.cell.totals)
      shf.cell.total.names = lapply(patient.names, function(x) paste(x,cell.type, sep = "."))
      shf.cell.total.sub = shf.cell.totals.update[, unlist(shf.cell.total.names)]
      weights = as.numeric(shf.cell.total.sub[1,]/rowSums(shf.cell.total.sub[1,]))
      
      weight.shf.OR[[cell.type]] = rowWeightedMeans(as.matrix(shf.mat), weights)
      names(weight.shf.OR[[cell.type]]) = shf.all$intron_junction
    }
    
    
    for (cell.type in all.cell.types){
      shf.diff = abs(weight.obs.OR[[cell.type]]) > abs(weight.shf.OR[[cell.type]])
      temp[[cell.type]] = temp[[cell.type]] + shf.diff 
    }
    
    #progress(value = x, max.value = nperm, progress.bar = T)
    Sys.sleep(0.01)
    if(x == nperm) cat(" Permuted differences calculated")  
  }
  
  pvals = list()
  for (cell.type in all.cell.types){
    pvals[[cell.type]] = as.data.frame(1 - temp[[cell.type]]/(nperm + 1))
    colnames(pvals[[cell.type]]) = cell.type
  }
  
  pvals.out = pvals %>% map(~ .x %>% rownames_to_column('intron_junction'))
  pvals.out = pvals.out %>% reduce(full_join, by = "intron_junction")
  colnames(pvals.out) = c("intron_junction", lapply(all.cell.types, function(x) paste(x, "p_value", sep = "_")))
  final = left_join(pvals.out,all, by = "intron_junction")
  return(final)
}

out = runJunctionPermutation(path.to.matrix = path.to.matrix, path.to.genotype = path.to.genotype, path.to.metadata = path.to.metadata, nperm)
message("Done with permutations!")

message("merging with metadata")
##merge with metadata 
all.metadata = metadata.list %>% reduce(full_join, by = c("intron_junction", "startVerdict", "endVerdict", "gene", "start_min_exon_distance","end_min_exon_distance", 
                                                          "distanceToMostCommonAltAnnotatedStart", "distanceToMostCommonAltAnnotatedEnd", "mostCommonAltAnnotationEndQualityTest",
                                                          "mostCommonAltAnnotationStartQualityTest", "numSkippedExons", "verdict", "start", "end", "chr"))

all.metadata = all.metadata %>% select(intron_junction, gene, startVerdict, endVerdict, mostCommonAltAnnotationStartQualityTest, mostCommonAltAnnotationEndQualityTest, start_min_exon_distance, end_min_exon_distance, distanceToMostCommonAltAnnotatedStart, distanceToMostCommonAltAnnotatedEnd, numSkippedExons, verdict)
all.metadata = all.metadata %>% separate(intron_junction, c("chr", "start", "end", "clusterID"), sep= ":") %>% 
  mutate(intron_junction = paste(chr,start,end, sep = ":"))

out = left_join(out, all.metadata, by = "intron_junction")

message("writing output")
setwd(output.dir)
save(out, file = "MDS_combined_celltype_junction_permutation_logOR_R_output.Rdata")
message("Done!!")
