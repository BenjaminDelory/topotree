qsmToTable<-function(inputqsm){

  if (mode(inputqsm)!="character"){stop("mode(inputqsm) must be character")}
  
  #Load QSM files
  
  filenames.qsm<-list.files(path=inputqsm, pattern="\\.txt$")
  path.qsm<-rep(inputqsm, length.out=length(filenames.qsm))
  filenamesqsm<-sub(x=filenames.qsm, pattern="\\.txt$", replacement="")
  message(paste("Number of qsm files in inputqsm:", length(filenames.qsm), sep=" "))
  
  if (length(filenames.qsm)==0){stop("There is no qsm file in inputqsm")}
  
  QSM<-lapply(paste(path.qsm, "/", filenames.qsm, sep=""), read.table, header=TRUE, sep=",")
  
  nodes<-0 #nodes is the number of rows (sum rows of each lie)

  for (i in 1:length(QSM)){nodes<-nodes+nrow(QSM[[i]])}
  
  #Construct table (1 line per segment)
  
  table<-matrix(nrow=nodes, ncol=17)
  #1 is file
  #2 is branch
  #3 is dbase
  #4 is dbasecum
  #5 is order
  #6 is bran
  #7 is apic
  #8 is x1
  #9 is y1
  #10 is z1
  #11 is x2
  #12 is y2
  #13 is z2
  #14 is length
  #15 is blength
  #16 is geodesic
  #17 is parentbranch
  
  n<-0 #n is the number of qsm files
  s<-0 #Count number of segments added to table
  
  for (i in 1:length(QSM)){ #For each tree
    
      n<-n+1 #Number of qsm files
      qsm<-QSM[[i]]
      branch<-0
      
      for (j in 1:nrow(qsm)){
        
        if (qsm[j,"parent"]==0){}
        
        else{
          
          s<-s+1
          
          #Order
          order<-qsm[j, "BranchOrder"]
          
          #Length
          length<-sqrt((qsm[qsm[j,"parent"], "start_1"]-qsm[j, "start_1"])^2+(qsm[qsm[j,"parent"], "start_2"]-qsm[j, "start_2"])^2+(qsm[qsm[j,"parent"], "start_3"]-qsm[j, "start_3"])^2)
          
          if (qsm[j, "branch"]>branch){
            
            bran<-1
            branch<-qsm[j, "branch"]
            cumsum<-length
            
            if (order==0){
              dbase<-0
              cumdbase<-0
              parentbranch<-0}
            
            else {
              
              parentbranch<-qsm[qsm[j,"parent"],"branch"]
              
              dbase<-table[which(table[,11]==qsm[qsm[j,"parent"],"start_1"] & table[,12]==qsm[qsm[j,"parent"],"start_2"] & table[,13]==qsm[qsm[j,"parent"],"start_3"]),15]
              cumdbase<-dbase+table[which(table[,11]==qsm[qsm[j,"parent"],"start_1"] & table[,12]==qsm[qsm[j,"parent"],"start_2"] & table[,13]==qsm[qsm[j,"parent"],"start_3"]),4]
              
              if (length(dbase)==0){
                dbase<-0
                cumdbase<-0}}}
          
          else{
            
            bran<-0
            cumsum<-cumsum+length}
      
          if (qsm[j,"extension"]==0) {apic<-1} else {apic<-0}
          
          x1<-qsm[qsm[j,"parent"], "start_1"]
          y1<-qsm[qsm[j,"parent"], "start_2"]
          z1<-qsm[qsm[j,"parent"], "start_3"]
          x2<-qsm[j, "start_1"]
          y2<-qsm[j, "start_2"]
          z2<-qsm[j, "start_3"]

          table[s,1:15]<-c(n, branch, dbase, cumdbase, order, bran, apic, x1, y1, z1, x2, y2, z2, length, cumsum)
          table[s,17]<-parentbranch
          }}}
  
  table[,16]<-table[,4]+table[,15] #geodesic distance
  index<-which(is.na(table[,1])==TRUE)
  table<-table[-index,-c(3,4)] #Remove lines and columns

  rownames(table)<-c(1:nrow(table))
  
  table<-as.data.frame(table)
  colnames(table)<-c("file","branch","order","bran","apic","x1","y1","z1","x2","y2","z2","length","blength","geodesic","parentbranch")
  table$file<-filenamesqsm[table$file]
  table$bran[table$bran==1]<-"true"
  table$bran[table$bran==0]<-"false"
  table$apic[table$apic==1]<-"true"
  table$apic[table$apic==0]<-"false"
  table<-table[,c(1,2,15,3:14)]
  
class(table)<-c("data.frame", "qsmToTable")
return(table)}