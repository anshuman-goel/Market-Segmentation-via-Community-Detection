#Loading Libraries
library(igraph)
library('proxy')

#Reading alpha value
args<-commandArgs(trailingOnly = TRUE)
alpha=as.numeric(args[1])

#Reading files and building graph
attrib<-read.csv("E:\\NCSU\\Semester 2\\Algorithms for Data Guided Business Intelligence\\Projects\\Market Segmentation via Community Detection\\06.Topic-7.Project-6.MarketSegmentation.AttributedGraphCommunityDetection\\data\\fb_caltech_small_attrlist.csv")
graph<-read.graph("E:\\NCSU\\Semester 2\\Algorithms for Data Guided Business Intelligence\\Projects\\Market Segmentation via Community Detection\\06.Topic-7.Project-6.MarketSegmentation.AttributedGraphCommunityDetection\\data\\fb_caltech_small_edgelist.txt","edgelist")


#Computing Similarity matrix
simA<<-as.matrix(simil(attrib,method = "cosine"))

#Phase 1
initial<<-1:vcount(graph)

phase1<-function(graph)
{
  iterations<-15
  community<-1:vcount(graph)
  new_mod<-modularity(graph,community)
  old_mod<--Inf
  while(old_mod!=new_mod && iterations>0)
  {
    for(i in 1:length(community)-1)
    {
      max_delta<--Inf
      max_j<-0
      mod<-modularity(graph,community)
      #print(mod)
      for(j in (i+1):length(community))
      {
          if(i!=j)
          {
            new_community<-rep(community)
            new_community[i]<-j
            #print(c(i,j,length(community)))
            delta<-(alpha*(modularity(graph,new_community)-mod)+(1-alpha)*simA[i,j])/length(unique(community))
            
            if(length(delta)!=0 && delta>max_delta && delta>0)
            {
              max_delta<-delta
              max_j<-j
            }
          }
      }
      if(max_j!=0)
      {
        community[i]<-max_j
      }

    }
    #print(community)
    old_mod<-new_mod
    new_mod<-modularity(graph,community)
    iterations<-iterations-1
    #break;
  }
  community
}

#Building Mapping to check for communities and building new graph with new vertices
mapping<-function(community)
{
  #print(community)
  for(i in 1:length(initial))
  {
    #print(initial[i])
    if(is.na(community[initial[i]]))
    {
      t<-initial[i]
      print(c(initial[i],community[initial[i]],t))
    }
    initial[i]<<-community[initial[i]]
  }

  count<-1
  mapping<-data.frame(matrix(0,nrow=length(unique(community)), ncol = 2))
  #print(mapping)
  for(i in 1:length(community))
  {
    flag<-FALSE
    if(count>1)
    {
      for(j in 1:(count-1))
      {
        if(community[i]==mapping[j,1])
        {
          flag<-TRUE
          break
        }
      }
      if(flag==FALSE)
      {
        mapping[count,1]<-community[i]
        mapping[count,2]<-count
        count<-count+1
      }
    }
    else
    {
      mapping[1,1]<-community[i]
      mapping[1,2]<-count
      count<-count+1
    }
  }
  for(i in 1:length(community))
  {
    for(j in 1:nrow(mapping))
    {
      if(community[i]==mapping[j,1])
      {
        community[i]<-mapping[j,2]
        break
      }
    }
  }
  #print(mapping[1,1])
  #print(length(initial))
  for(i in 1:length(initial))
  {
    #print(mapping)
    for(j in 1:(count-1))
    {
      #print(c(initial[i],mapping[j,1]))
      if(initial[i]==mapping[j,1])
      {
        initial[i]<<-mapping[j,2]
        break
      }
    }
  }
  #print(max(initial))
  community
}

#Updating Similarity Matrix
updatingsimA<-function(community)
{
  for(i in 1:length(community))
  {
    temp<-simA[i,community[i]]+simA[community[i],i]
    simA[i,community[i]]<<-temp
    simA[community[i],i]<<-temp
  }
}

#Phase 2
phase2<-function()
{
  #iterations<-5
  
  while(length(unique(initial))>9)
  {
    community<-phase1(graph)
    mapping<-mapping(community)
    print(length(unique(initial)))
    #print(unique(community))
    graph<-contract.vertices(graph,mapping)
    graph<-simplify(graph, remove.multiple = TRUE, remove.loops = TRUE)
    updatingsimA(community)
    #print(graph)
    #print(vcount(graph))
  }
  
  #Writing to the file
  
  fileName<-paste("E:\\NCSU\\Semester 2\\Algorithms for Data Guided Business Intelligence\\Projects\\Market Segmentation via Community Detection\\06.Topic-7.Project-6.MarketSegmentation.AttributedGraphCommunityDetection\\communities",alpha,sep="_")
  fileName<-paste(fileName,"txt",sep=".")
  fileptr<-file(fileName,"w")
  
  for(i in 1:max(initial))
  {
    ptr<-vector("numeric")
    for(j in 1:length(initial))
    {
      if(initial[j]==i){
        ptr<-append(ptr,j-1,after=length(ptr))
      }
    }
    cat(as.character(ptr),file=fileptr,sep = ",")
    cat("\n",file=fileptr)
  }
  close(fileptr)
}

phase2()