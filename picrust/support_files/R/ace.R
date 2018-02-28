#!/usr/bin/env Rscript
#./find_ace.R  <input: reference_tree_with_internal_nodes_labelled> <input: functional_matrix_file> <ace_method> <output: count_table> <output: prob_table>

library(ape)
Args <- commandArgs(TRUE)

#load in the tree
tree<-read.tree(Args[1])

#load in the trait table
data<-read.delim(Args[2],check.names=FALSE,row.names=1)


asr_method=Args[3]
count_out_file=Args[4]
ci_out_file=Args[5]

#If tips in tree contain '_' then read.tree() places single quotes around these tip labels.
#This then causes sorting errors below since the rownames are different between the trait table and the tree.
#Fix this by putting quotes around any labels in the trait table that have a '_'.
for(i in grep("_",rownames(data))){
 rownames(data)[i]<-paste("'",rownames(data)[i],"'",sep="")
}


#order the trait table to match the tree tip labels
#The option "drop=FALSE" means that a dataframe will be returned even if there 
#is only 1 column (instead of a vector, which is the default).
data_ordered <- data[tree$tip.label , , drop=FALSE]

#do the actual ace reconsructions
reconstructions<-apply(data_ordered,2,ace,tree, type="continuous",method=asr_method)

#pull out only the ace node predictions
just_ace<-lapply(1:length(reconstructions),function(x) reconstructions[[x]]$ace)
names(just_ace)<-names(reconstructions)

#reformat the list into a matrix
just_ace_matrix<-do.call(cbind,just_ace)

#relabel the node names (ones created internally by ape) with the actual node labels in the tree
just_ace_matrix<-cbind(tree$node.label,just_ace_matrix)

#give a simple header label for the internal nodes
colnames(just_ace_matrix)[1]<-'nodes'

#Convert to a data frame
out_matrix<-data.frame(just_ace_matrix,check.names=FALSE)

#write to file
write.table(out_matrix,file=count_out_file,row.names=FALSE,quote=FALSE, sep="\t")


#extract 95% CI info
ci<-lapply(1:length(reconstructions),function(x) paste(round(reconstructions[[x]]$CI95[,1],digits=4),round(reconstructions[[x]]$CI95[,2],digits=4),sep="|"))
names(ci)<-names(reconstructions)                                                                                
ci_matrix<-do.call(cbind,ci)
ci_matrix<-cbind(tree$node.label,ci_matrix)

if(asr_method=="ML" || asr_method=="REML"){
  #extract information about brownian motion parameter
  sigma<-lapply(1:length(reconstructions),function(x) paste(round(reconstructions[[x]]$sigma2[1],digits=4),round(reconstructions[[x]]$sigma2[2],digits=4),sep="|"))
  names(sigma)<-names(reconstructions)                                                                                
  sigma_matrix<-do.call(cbind,sigma)
  sigma_matrix<-cbind(c('sigma'),sigma_matrix)
  
  #add it to the ci matrix
  ci_matrix<-rbind(ci_matrix,sigma_matrix)

  #get loglik value
  if(asr_method=="ML"){
    loglik<-lapply(1:length(reconstructions),function(x) round(reconstructions[[x]]$loglik,digits=4))                                                                                                                           
  }else{
    loglik<-lapply(1:length(reconstructions),function(x) round(reconstructions[[x]]$resloglik,digits=4))                                                                                                                                                                                        
  }

  names(loglik)<-names(reconstructions)
  loglik_matrix<-do.call(cbind,loglik)
  loglik_matrix<-cbind(c('loglik'),loglik_matrix)
  ci_matrix<-rbind(ci_matrix,loglik_matrix)

}
  
#just set the column name for the row names to someting arbtrary like 'nodes'
colnames(ci_matrix)[1]<-'nodes'

#output the data to file
out_matrix<-data.frame(ci_matrix,check.names=FALSE)                                                                                                                 
write.table(out_matrix,file=ci_out_file,row.names=FALSE,quote=FALSE, sep="\t")
