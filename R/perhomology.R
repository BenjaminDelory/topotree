perhomology<-function(x){
  
  #x must be an object of class "qsmToTable"
  if ("qsmToTable" %in% class(x)) {} else {stop("x must be a qsmToTable object")}
  
  #Compute persistent homology (geodesic distance function)
  
    results<-matrix(nrow=sum(x$apic=="true"), ncol=3) #Create matrix to store the results (birth and death of each axis)
    colnames(results)<-c("dimension", "birth", "death") #Rename columns
    x<-x[order(x$geodesic, decreasing=TRUE),]#Order x lines by geodesic values
    x$hzero<-rep(NA, nrow(x))#Create a new column in x to store the ID number of each zero order homology
    
    apicindex<-which(x$apic=="true")
    results[,1]<-rep(0, nrow(results))
    
    for (l in 1:length(apicindex)){#For each apic point of a tree
      
      branch<-x$branch[apicindex[l]]
      parentbranch<-x$parentbranch[apicindex[l]]
      
      x$hzero[apicindex[l]]<-l
      
      results[l,2]<-x$geodesic[apicindex[l]]
      results[l,3]<-x$geodesic[apicindex[l]]-x$length[apicindex[l]]
      
      indexprec<-which(x$x2==x$x1[apicindex[l]] & x$y2==x$y1[apicindex[l]] & x$z2==x$z1[apicindex[l]])
      if (length(indexprec)>1) {indexprec<-indexprec[which(x$branch[indexprec]==branch | x$branch[indexprec]==parentbranch)]}

      if (length(indexprec)>0){
        
        branch<-x$branch[indexprec]
        parentbranch<-x$parentbranch[indexprec]
      
          while(is.na(x$hzero[indexprec])==TRUE){
            
            x$hzero[indexprec]<-l
            
            results[l,3]<-x$geodesic[indexprec]-x$length[indexprec]
            
            segment1<-which(x$x2==x$x1[indexprec] & x$y2==x$y1[indexprec] & x$z2==x$z1[indexprec])
            if (length(segment1)>1) {indexprec<-segment1[which(x$branch[segment1]==branch | x$branch[segment1]==parentbranch)]} else {indexprec<-segment1}
            if (length(indexprec)==0){break}
            branch<-x$branch[indexprec]
            parentbranch<-x$parentbranch[indexprec]}}}
    
    results<-results[order(results[,3], decreasing=FALSE),]
    class(results)<-c("matrix", "barcode")
  
  return(results)}