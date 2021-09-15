PseudotimeWindow = function(archrProject, pseudotime.column, treatment.column, treatment.levels, step, cluster.column = NULL, cluster.count = NULL){
  if(class(archrProject) != "ArchRProject"){
    stop("Specified object is not an ArchrProject")
  }
  if(!"MotifMatrix" %in% getAvailableMatrices(archrProject)){
    stop("No MotifMatrix found in ArchRProject")
  }
  if(!pseudotime.column %in% colnames(archrProject@cellColData)){
    stop("Specified pseudotime column not found in @cellColData slot")
  }
  if(!(treatment.column %in% colnames(archrProject@cellColData) & cluster.column %in% colnames(archrProject@cellColData))){
    stop("Specified treatment column or cluster column not found in @cellColData slot")
  }
  if(!sum(treatment.levels %in% unique(archrProject@cellColData[,treatment.column])) == 2){
    stop("Specified treatment levels not found in specified treatment column in @cellColData slot")
  }
  
  # Get motif deviation matrix
  mm <- getMatrixFromProject(archrProject, useMatrix = "MotifMatrix")
  m <- mm@assays@data$deviations
  
  # Order cells in pseudotime
  meta = as.data.frame(archrProject@cellColData)
  ncol(m[,!is.na(meta[,pseudotime.column])])
  m.sub = m[,!is.na(meta[,pseudotime.column])]
  m.meta = meta[!is.na(meta[,pseudotime.column]),]
  ord = rownames(m.meta)[order(m.meta[,pseudotime.column],decreasing = F)]
  m.meta = m.meta[ord,]
  m.sub = m[,ord]
  
  # Get sliding windows
  step.def = seq(from = 1, to = ncol(m.sub), by = step)
  define.windows = list()
  for(i in 1:(round(ncol(m.sub)/(step))-1)){
    define.windows[[i]] = c(step.def[i],step.def[i+10])
  }
  ind = lapply(define.windows, function(x) sum(is.na(x)) > 0)
  define.windows = define.windows[!unlist(ind)]
  
  # Get mean accessibility per window for each treatment
  access1 = lapply(define.windows, function(x){
    temp = m.sub[,c(x[1]:x[2])]
    meta.temp = m.meta[c(x[1]:x[2]),]
    wt = rownames(meta.temp)[which(as.character(meta.temp[,treatment.column]) == treatment.levels[2])]
    rowMeans(m.sub[,colnames(m.sub) %in% wt])
  })
  access2 = lapply(define.windows, function(x){
    temp = m.sub[,c(x[1]:x[2])]
    meta.temp = m.meta[c(x[1]:x[2]),]
    mut = rownames(meta.temp)[which(as.character(meta.temp[,treatment.column]) == treatment.levels[1])]
    rowMeans(m.sub[,colnames(m.sub) %in% mut])
  })
  access1 = t(do.call(rbind,access1))
  access2 = t(do.call(rbind,access2))
    
  # Get differential accessibility between treatments for each window
  win.access = lapply(define.windows, function(x){
    temp = m.sub[,c(x[1]:x[2])]
    meta.temp = m.meta[c(x[1]:x[2]),]
    wt = rownames(meta.temp)[which(as.character(meta.temp[,treatment.column]) == treatment.levels[2])]
    mut = rownames(meta.temp)[which(as.character(meta.temp[,treatment.column]) == treatment.levels[1])]
    message(paste0(length(wt)," ",treatment.levels[2]," cells"))
    message(paste0(length(mut)," ",treatment.levels[1]," cells"))
    message(paste0(x))
    out = rowMeans(m.sub[,colnames(m.sub) %in% mut]) - rowMeans(m.sub[,colnames(m.sub) %in% wt])
  })
  out = t(do.call(rbind, win.access))
  colnames(out) = as.character(1:ncol(out))
  
  # Get fraction of cells in treatment 2 vs treatment 1
  mut.fraction = unlist(lapply(define.windows, function(x){
    meta.temp = m.meta[c(x[1]:x[2]),]
    wt = rownames(meta.temp)[which(as.character(meta.temp[,treatment.column]) == treatment.levels[2])]
    mut = rownames(meta.temp)[which(as.character(meta.temp[,treatment.column]) == treatment.levels[1])]
    out = length(mut)/(length(wt)+ length(mut))
  }))
  
  # Get mean pseudotime value per window
  pseudotime.mean = unlist(lapply(define.windows, function(x){
    meta.temp = m.meta[c(x[1]:x[2]),]
    mean(meta.temp[,pseudotime.column])
  }))
  
  if(!is.null(cluster.column) & !is.null(cluster.count)){
    # Get fraction of HSCs in window
    hsc.fraction = unlist(lapply(define.windows, function(x){
      meta.temp = m.meta[c(x[1]:x[2]),]
      sum(meta.temp[,cluster.column] == cluster.count) / nrow(meta.temp)
    }))
  }
  
  annots = data.frame(Pseudotime = pseudotime.mean,
                      MutantFraction = mut.fraction,
                      HSCFraction = hsc.fraction)
  
  ord = sort(apply(out,1,function(x) which.max(x)), decreasing = F)
  out = out[names(ord),]
  
  return(list(Windows = out, MetaData = annots, Group1 = access1, Group2 = access2))
  
}
