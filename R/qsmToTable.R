qsmToTable<-function(inputqsm, export=NULL){

  if (mode(inputqsm)!="character"){stop("inputqsm must be a character string")}
  if (is.null(export)==FALSE){if (mode(export)!="character"){stop("export must be a character string")}}

  #Load QSM file
  
  qsm<-read.table(inputqsm, header=TRUE, sep=",", dec=".")
  
  name<-sub(x=basename(inputqsm), pattern="\\.txt$", replacement="")
  
  #Construct table (1 line per segment)
  
  table<-matrix(nrow=nrow(qsm), ncol=17)
  
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
  
  s<-0 #Count number of segments added to table
  
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
              
              dbase<-table[which(table[,11]==qsm[qsm[j,"parent"],"start_1"] & table[,12]==qsm[qsm[j,"parent"],"start_2"] & table[,13]==qsm[qsm[j,"parent"],"start_3"] & table[,14]!=0), 15]
              cumdbase<-dbase+table[which(table[,11]==qsm[qsm[j,"parent"],"start_1"] & table[,12]==qsm[qsm[j,"parent"],"start_2"] & table[,13]==qsm[qsm[j,"parent"],"start_3"] & table[,14]!=0),4]
              
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

          table[s,1:15]<-c(1, branch, dbase, cumdbase, order, bran, apic, x1, y1, z1, x2, y2, z2, length, cumsum)
          table[s,17]<-parentbranch
          }}
  
  table[,16]<-table[,4]+table[,15] #geodesic distance
  index<-which(is.na(table[,1])==TRUE)
  table<-table[-index,-c(3,4)] #Remove lines and columns
  
  table<-as.data.frame(table)
  colnames(table)<-c("file","branch","order","bran","apic","x1","y1","z1","x2","y2","z2","length","blength","geodesic","parentbranch")
  table$file<-name
  table$bran[table$bran==1]<-"true"
  table$bran[table$bran==0]<-"false"
  table$apic[table$apic==1]<-"true"
  table$apic[table$apic==0]<-"false"
  table<-table[,c(1,2,15,3:14)]
  
  #Remove segments with length=0
  
  index<-rev(which(table$length==0))
  
  if (length(index)>0){
    
    for (i in 1:length(index)){
      
      if (table$apic[index[i]]=="false" & table$bran[index[i]]=="false" & table$length[index[i]]==0){
        
        message(paste("Segment (length=0) removed from branch ", table$branch[index[i]], " in ", table$file[index[i]], sep=""))
        table<-table[-index[i],]}
      
      if (table$apic[index[i]]=="true" & table$bran[index[i]]=="false" & table$length[index[i]]==0){
        
        message(paste("Apex segment (length=0) removed from branch ", table$branch[index[i]], " in ", table$file[index[i]], sep=""))
        table$apic[index[i]-1]<-"true"
        table<-table[-index[i],]}
      
      if (table$apic[index[i]]=="false" & table$bran[index[i]]=="true" & table$length[index[i]]==0){
        
        message(paste("Bran segment (length=0) removed from branch ", table$branch[index[i]], " in ", table$file[index[i]], sep=""))
        table$bran[index[i]+1]<-"true"
        table<-table[-index[i],]}
      
      if (table$apic[index[i]]=="true" & table$bran[index[i]]=="true" & table$length[index[i]]==0){
        
        message(paste("Bran/Apex segment (length=0) removed from branch ", table$branch[index[i]], " in ", table$file[index[i]], sep=""))
        table<-table[-index[i],]}}}
  
  rownames(table)<-c(1:nrow(table))
  
  #Extend branches when needed
  
  apicindex<-which(table$apic=="true")
  branindex<-which(table$bran=="true")
  
  for (i in 1:length(apicindex)){
    
    index<-which(table[branindex, "x1"]==table[apicindex[i], "x2"] & table[branindex, "y1"]==table[apicindex[i], "y2"] & table[branindex, "z1"]==table[apicindex[i], "z2"])
    
    if (length(index)>0){ #Extend apical branch
      
      table[nrow(table)+1,]<-table[apicindex[i],] #Add new segment
      table$apic[apicindex[i]]<-"false"
      table$bran[nrow(table)]<-"false"
      table$apic[nrow(table)]<-"true"
      table[nrow(table), c("x1", "y1", "z1")]<-table[apicindex[i], c("x2", "y2", "z2")]
      x1<-table[apicindex[i], "x1"]
      y1<-table[apicindex[i], "y1"]
      z1<-table[apicindex[i], "z1"]
      x2<-table[apicindex[i], "x2"]
      y2<-table[apicindex[i], "y2"]
      z2<-table[apicindex[i], "z2"]
      
      #Have to consider different scenarios
      
      diff<-c(x1-x2, y1-y2, z1-z2)
      
      if (sum(diff==0)==0){
      
          K<-((x2-x1)/(z2-z1))^2+((y2-y1)/(z2-z1))^2+1
          newz<-sqrt(0.001^2/K)+z2
          newx<-((x2-x1)/(z2-z1))*(newz-z2)+x2
          newy<-((y2-y1)/(z2-z1))*(newz-z2)+y2}
      
      if (sum(diff==0)==1){
        
        a<-which(diff==0)
        
        if (a==1){
          K<-((z2-z1)/(y2-y1))^2+1
          newx<-x2
          newy<-sqrt(0.001^2/K)+y2
          newz<-((z2-z1)/(y2-y1))*(newy-y2)+z2} 
        
        if (a==2){
          K<-((x2-x1)/(z2-z1))^2+1
          newy<-y2
          newz<-sqrt(0.001^2/K)+z2
          newx<-((x2-x1)/(z2-z1))*(newz-z2)+x2}
        
        if (a==3){
          K<-((x2-x1)/(y2-y1))^2+1
          newz<-z2
          newy<-sqrt(0.001^2/K)+y2
          newx<-((x2-x1)/(y2-y1))*(newy-y2)+x2}}
      
      if (sum(diff==0)==2){
        
        a<-which(diff!=0)
        
        if (a==1){
          newx<-x2+0.001
          newy<-y2
          newz<-z2}
        
        if (a==2){
          newy<-y2+0.001
          newx<-x2
          newz<-z2}
        
        if (a==3){
          newz<-z2+0.001
          newx<-x2
          newy<-y2}}
        
        length<-sqrt((newx-x2)^2+(newy-y2)^2+(newz-z2)^2)
        table[nrow(table), c("x2", "y2", "z2")]<-c(newx, newy, newz)
        table[nrow(table), "length"]<-length
        table[nrow(table), "blength"]<-table[nrow(table), "blength"]+length
        table[nrow(table), "geodesic"]<-table[nrow(table), "geodesic"]+length
      
      message(paste("Branch ", table$branch[apicindex[i]], " in ", table$file[apicindex[i]]," had to be extended", sep=""))}}
  
  table<-table[order(table$branch, table$blength), ]
  
  rownames(table)<-c(1:nrow(table))
  
  if (is.null(export)==FALSE){
    write.table(table, file=paste(export, "/", name, "_R.txt", sep=""), col.names = TRUE, row.names = FALSE, sep=",")}
  
class(table)<-c("data.frame", "qsmToTable")
return(table)}