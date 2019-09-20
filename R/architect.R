architect<-function(inputqsm){
  
    if (is.null(inputqsm)==FALSE){if ("qsmToTable" %in% class(inputqsm)){} else {stop("inputqsm must be a qsmToTable object")}}
  
    maxord<-max(inputqsm$order) #Tree trunk has order zero
  
    n<-length(unique(inputqsm$file))

    dataqsm<-data.frame(FileName=rep(NA, n), TBL=rep(NA, n), L0B=rep(NA, n), TN0B=rep(NA, n), TNLB=rep(NA, n), TLBL=rep(NA, n), D1LB=rep(NA, n), Height=rep(NA, n), Width=rep(NA, n))
    
    if (maxord>0){latbranch<-matrix(ncol=3*maxord, nrow=n)}
    
    files<-unique(inputqsm$file)
    
    k<-0
    
    for (i in 1:length(files)){
        
        xt<-inputqsm[inputqsm$file==files[i],] #x is a subset of the qsmToTable object
        
        k<-k+1
        
        #File name
        dataqsm$FileName[k]<-unique(xt$file)
        
        #TBL (total branch length)
        dataqsm$TBL[k]<-sum(xt$length)
        
        #L0B (total zero-order branch length)
        dataqsm$L0B[k]<-sum(xt$length[xt$order==0])

        #TN0R (number of zero-order branches)
        dataqsm$TN0B[k]<-length(unique(xt$branch[xt$order==0]))
        
        #TNLB (number of lateral branches)
        dataqsm$TNLB[k]<-length(unique(xt$branch[xt$order>0]))
        
        #TLBL (total lateral branch length)
        dataqsm$TLBL[k]<-sum(xt$length[xt$order>0])
        
        #D1LB (density of first-order branches on the zero-order branch)
        dataqsm$D1LB[k]<-length(unique(xt$branch[xt$order==1]))/dataqsm$L0B[k]
        
        #Latbranch (number, length, mean length)
        
        if (maxord>0){
          
          for (l in 1:maxord){
            
            latbranch[k,l]<-length(unique(xt$branch[xt$order==l])) #Number of lateral branches
            latbranch[k,l+maxord]<-sum(xt$length[xt$order==l]) #Length of lateral branches
            if (latbranch[k,l]==0) {latbranch[k,l+2*maxord]<-0} else {latbranch[k,l+2*maxord]<-latbranch[k,l+maxord]/latbranch[k,l]}}} #Mean length of lateral branches
        
        #Height
        dataqsm$Height[k]<-abs(max(xt$z2)-min(xt$z1))
        
        #Width
          
        widthy<-abs(max(xt$y2)-min(xt$y2))
        widthx<-abs(max(xt$x2)-min(xt$x2))
        dataqsm$Width[k]<-max(c(widthy, widthx))}
        
        # Results in a dataframe
        
        if (maxord>0){
          
          latbranch<-as.data.frame(latbranch)
          
          for (l in 1:maxord){
            colnames(latbranch)[l]<-paste("N", l, "LB", sep="")
            colnames(latbranch)[l+maxord]<-paste("L", l, "LB", sep="")
            colnames(latbranch)[l+2*maxord]<-paste("ML", l, "LB", sep="")}
          
        dataqsm<-data.frame(dataqsm[,1:6], as.data.frame(latbranch), dataqsm[,7:9])}
  
  return(dataqsm)}