perhomology<-function(x, show.progress=FALSE){
  
  #x must be an object of class "qsmToTable"
  if ("qsmToTable" %in% class(x)) {} else {stop("x must be a qsmToTable object")}
  
  if (mode(show.progress)!="logical"){stop("show.progress must be logical")}
  
  x$tree<-x$file
  RSlevels<-unique(x$tree)
  
  n<-length(RSlevels) #Number of trees in x
  
  #Compute persistent homology (geodesic distance function)
  
  if (show.progress==TRUE) {pb<-txtProgressBar(min=1, max=n, style=3)}
  
  results<-list()
  
  for (i in 1:n){#For each tree in the table
    
    if (show.progress==TRUE) {setTxtProgressBar(pb, i)}
    
    table<-x[x$tree==RSlevels[i],] #Create a table subset for computing persistent homology
    results[[i]]<-matrix(nrow=sum(table$apic=="true"), ncol=3) #Create matrix to store the results (birth and death of each axis)
    colnames(results[[i]])<-c("dimension", "birth", "death") #Rename columns
    table<-table[order(table$geodesic, decreasing=TRUE),]#Order table lines by geodesic values
    table$hzero<-rep(NA, nrow(table))#Create a new column in table to store the ID number of each zero order homology
    
    apicindex<-which(table$apic=="true")
    results[[i]][,1]<-rep(0, nrow(results[[i]]))
    
    for (l in 1:length(apicindex)){#For each apic point of a tree
      
      branch<-table$branch[apicindex[l]]
      parentbranch<-table$parentbranch[apicindex[l]]
      
      table$hzero[apicindex[l]]<-l
      
      results[[i]][l,2]<-table$geodesic[apicindex[l]]
      results[[i]][l,3]<-table$geodesic[apicindex[l]]-table$length[apicindex[l]]
      
      indexprec<-which(table$x2==table$x1[apicindex[l]] & table$y2==table$y1[apicindex[l]] & table$z2==table$z1[apicindex[l]])
      if (length(indexprec)>1) {indexprec<-indexprec[which(table$branch[indexprec]==branch | table$branch[indexprec]==parentbranch)]}

      if (length(indexprec)>0){
        
        branch<-table$branch[indexprec]
        parentbranch<-table$parentbranch[indexprec]
      
          while(is.na(table$hzero[indexprec])==TRUE){
            
            table$hzero[indexprec]<-l
            
            results[[i]][l,3]<-table$geodesic[indexprec]-table$length[indexprec]
            
            segment1<-which(table$x2==table$x1[indexprec] & table$y2==table$y1[indexprec] & table$z2==table$z1[indexprec])
            if (length(segment1)>1) {indexprec<-segment1[which(table$branch[segment1]==branch | table$branch[segment1]==parentbranch)]} else {indexprec<-segment1}
            if (length(indexprec)==0){break}
            branch<-table$branch[indexprec]
            parentbranch<-table$parentbranch[indexprec]}}}
    
    results[[i]]<-results[[i]][order(results[[i]][,3], decreasing=FALSE),]
    class(results[[i]])<-c("matrix", "barcode")}
  
  names(results)<-as.character(RSlevels)
  class(results)<-c("perhomology", "list")
  return(results)}