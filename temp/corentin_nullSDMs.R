## Both scripts may require the following libraries, and options to be loaded in R:
Sys.setenv(NOAWT=TRUE) #only on some macs before loading rJava
library(rJava) # required to run Maxent within dismo
#optional, to set the memory allocated to java, hence maxent (before loading dismo)
options(java.parameters = "-Xmx2g" )
library(dismo) # should also load automatically the required packages sp and raster
## 1) The main null model function:
mxt.nulltest(train, bg, n = 9, proj=bg, c.proj=TRUE, args = "addsamplestobackground", threshold = .1, abs.diff=TRUE, test=NULL, group=NULL, bg.group=NULL,save=NULL)
## The inputs for this function are:
## - train : a dataframe of the environmental values for the species calibration presences: one row for each presence, one column for each predictor. Note that categorical variable must be indicated with as.factor() .
## - n : the number of null replicates.
## - bg : the same as “train” but for the background points (the comparisons records for the model).
## - proj : a dataframe in the same format as “train” and “bg” which corresponds to the pixels over which the statistics will be calculated. By default this is the same as “bg” but it may be useful when the models are projected to a different region.
## - c.proj : TRUE, if “proj” includes both the calibration and evaluation records, FALSE, if it includes only the evaluation records (note that when FALSE is selected, AUCdiff will be meaningless).
## - args: a list of Maxent arguments.
## - abs.diff : TRUE or FALSE, whether the absolute value of AUCdiff should be returned.
## - threshold : the percentile of calibration records used to set the threshold to calculate the omission error rate.
## - test: the same as “train” but for the evaluation records. If it is not provided, a cross validation test procedure will be used.
## - group : (only if test is missing), a vector of length nrow(train) which indicates the partitions for the cross validation procedure. If omitted the presence records are assigned to partitions randomly (with a minimum of 10 records per partitions, and a maximum of 10 partitions).
## - bg.group : (only if test is missing) the same as “group” but for the background points. If omitted these are assigned to partitions randomly.
## - save: the path for the .csv file where the results should be saved (optional) after each Maxent run (for both the real data and the null replicates).
## The function returns an object with the following slots:
## @summary : a data.frame with the performance and significance of summary statistics.
## @random.reps : an array with the performance statistics for the null replicates.
## Note: when a save file path is not provided, the results of this function should be saved to a R object for further analyses.
## The function code:
mxt.nulltest = function(train, bg, n = 9, proj=bg, c.proj=TRUE, args = "addsamplestobackground", threshold = .1, abs.diff=TRUE, test=NULL, group=NULL, bg.group=NULL,save=NULL) {
  require(ROCR) # the library ’ROCR’ has to be previously installed
  # the following sub-function runs the Maxent models and calculates the evaluation statistics
  f.real= function(x, p, args, proj, c.proj, bg, train, test, abs.diff, threshold) {
    tmpdir=paste(tempdir(),runif(1,0,1),sep="/")
    dir.create(tmpdir, showWarnings = TRUE, recursive = FALSE)
    mod=maxent(x,p, args=args,path=tmpdir)
    test.bb = predict(mod,proj)
    if(c.proj==TRUE) train.bb=test.bb else train.bb=predict(mod,bg)
    testpp = predict(mod,test)
    trainpp = predict(mod,train)
    combinedtest = c(testpp, test.bb)
    labeltest = c(rep(1,length(testpp)), rep(0,length(test.bb)))
    combinedtrain = c(trainpp, train.bb)
    labeltrain = c(rep(1,length(trainpp)), rep(0,length(train.bb)))
    predtest =  prediction(combinedtest, labeltest)
    predtrain =  prediction(combinedtrain, labeltrain)
    AUC= function(x) {performance(x, "auc")@y.values[[1]]}
    AUCtest = AUC(predtest)
    AUCtrain = AUC(predtrain)
    AUCdiff= AUCtrain - AUCtest
    if(abs.diff==FALSE) AUCdiff=AUCdiff else AUCdiff=abs(AUCdiff)
    r=length(testpp[testpp<quantile(trainpp,threshold)])
    t=length(testpp)
    a=length(test.bb[test.bb>=quantile(trainpp,threshold)])/length(test.bb)
    bintest=binom.test(t-r,t,p=a,alternative = "g")
    OR_binom_p=bintest$p.value
    OR= r/t
    stats= c(AUCtrain, AUCtest, AUCdiff, OR, OR_binom_p)
    unlink(tmpdir,recursive=T,force=T)
    return(stats)}
  # the rest of the function iterates the sub-function above for each of the real and null models; it prints on screen and (optionally) saves the results after each iteration.
  if(missing(test)) {
    if(missing(group)) {
      if(nrow(train)<100) {k=trunc(nrow(train)/10)} else {k=10}
      group=kfold(train,k)} else {group=group
      k=length(unique(group))}
    if(missing(bg.group)) bg.group=kfold(bg, length(unique(group))) else bg.group=bg.group
    lbl=vector()
    for(i in 1:k) {lbl[[i]]=paste("Real",i,sep="_")}
    a=k-1
    for(i in 1:k) {for(j in 1:n) {lbl[[a+i+j]]=paste("Null",i,j,sep="_")}
      a=a+n-1}
    res = array(,dim=c(k*(n+1)+5,5), dimnames = list(c(lbl, "Mean_Real", "Mean_Null", "Std.Dev_Null", "Z_score", "p_value"), c("AUCtrain", "AUCtest", "AUCdiff", "OR", "OR_binom_p")))
    cat("Replicate",colnames(res),"  % Completion","\n", sep ="\t")
    null = array(,dim=c(n,5,k), dimnames = list(NULL, c("AUCtrain", "AUCtest", "AUCdiff", "OR", "OR_binom_p"), NULL))
    a=k-1
    for (i in 1:k) {
      k.train = train[group != i,]
      k.test = train[group == i,]
      k.bg = bg[bg.group != i,]
      x=rbind(k.train, k.bg)
      p=c(rep(1,nrow(k.train)),rep(0,nrow(k.bg)))
      reps.x= list()
      reps =list()
      for (j in 1:n) {
        s = sample(nrow(k.bg),nrow(k.train))
        reps[[j]]=k.bg[s,]
        if(any(args=="noaddsamplestobackground")) {
          reps.x[[j]] = rbind(k.bg[s,], k.bg)
          reps.p = p} else {
            reps.x[[j]] = rbind(k.bg[s,], k.bg[-s,])
            reps.p = c(rep(1,length(s)),rep(0,(nrow(k.bg)-length(s))))}
      }
      res[i,]=f.real(x, p, args, proj, c.proj, bg, train=k.train, test=k.test, abs.diff, threshold)
      cat(rownames(res)[i],res[i,],"  ",round((i*(n+1)-n)/((n+1)*k+1)*100,digits=2),"%","\n", sep ="\t")
      for (j in 1:n) {
        null[j,,i]=f.real(reps.x[[j]],p= reps.p, args, proj, c.proj, bg, train=reps[[j]], test=k.test, abs.diff, threshold)
        res[a+i+j,]=null[j,,i]
        cat(rownames(res)[a+i+j],res[a+i+j,],"  ",round((i*(n+1)-n+j)/((n+1)*k+1)*100,digits=2),"%","\n", sep ="\t")
        if(!missing(save)) write.csv(res,save)
      }
      a=a+n-1}
    res[k*(n+1)+1,]=apply(res[1:k,],2,mean)
    means.k.null = rowMeans(null, dims=2)
    res[k*(n+1)+2,]=colMeans(means.k.null)
    res[k*(n+1)+3,]= apply(means.k.null,2,sd)
    res[k*(n+1)+4,]=(res[k*(n+1)+1,]-res[k*(n+1)+2,])/res[k*(n+1)+3,]
    res[k*(n+1)+5,1:2]=(1-pnorm(res[k*(n+1)+4,1:2]))
    res[k*(n+1)+5,3:5]=(pnorm(res[k*(n+1)+4,3:5]))
    cat(rownames(res)[k*(n+1)+5],res[k*(n+1)+5,],"  100 %","\n", sep ="\t")
    real=res[-((k+1):(k*(n+1))),]
  } else {
    x=rbind(train, bg)
    p=c(rep(1,nrow(train)),rep(0,nrow(bg)))
    reps.x= list()
    reps =list()
    for (i in 1:n) {
      s = sample(nrow(bg),nrow(train))
      reps[[i]]=bg[s,]
      if(any(args=="noaddsamplestobackground")) {
        reps.x[[i]] = rbind(bg[s,], bg)
        reps.p = p} else {
          reps.x[[i]] = rbind(bg[s,], bg[-s,])
          reps.p = c(rep(1,length(s)),rep(0,(nrow(bg)-length(s))))
        }
    }
    lbl=vector()
    for (i in 1:n) {lbl[[i]]=paste("Null",i,sep="_")}
    res = array(,dim=c(n+5,5), dimnames = list(c("Real", lbl, "Mean_Null", "Std.Dev_Null", "Z_score", "p_value"), c("AUCtrain", "AUCtest", "AUCdiff", "OR", "OR_binom_p")))
    cat("Replicate",colnames(res),"  % Completion","\n", sep ="\t")
    res[1,]=f.real(x, p, args, proj, c.proj, bg, train=train, test=test, abs.diff, threshold)
    cat(rownames(res)[1],res[1,],"  1 %","\n", sep ="\t")
    for(i in 1:n) {res[i+1,]= f.real(reps.x[[i]], p= reps.p, args, proj, c.proj, bg, train=reps[[i]], test=test, abs.diff, threshold)
    cat(rownames(res)[i+1],res[i+1,],"  ", round((i+1)/(n+2)*100,digits=2), "%", "\n", sep ="\t")
    if(!missing(save)) write.csv(res,save)
    }
    null=res[2:(n+1),]
    res[n+2,]=apply(null,2,mean)
    res[n+3,]=apply(null,2,sd)
    res[n+4,]=(res[1,]-res[n+2,])/res[n+3,]
    res[n+5,1:2]=1-pnorm(res[n+4,1:2])
    res[n+5,3:5]=pnorm(res[n+4,3:5])
    cat(rownames(res)[n+5],res[n+5,],"  100 %","\n", sep ="\t")
    real=res[-(2:(n+1)),]
  }
  if(!missing(save)) write.csv(res,save)
  Results<- setClass("Results", representation(summary = "data.frame", random.reps = "array"))
  res=new("Results")
  res@summary = data.frame(real)
  res@random.reps = null
  return(res)}
## 2) The script to automate in parallel the optimization procedure:
# Create a list of maxent arguments corresponding to those settings 42 different settings:
a="addsamplestobackground"
b="noautofeature"
h="nohinge"
t="nothreshold"
p="noproduct"
LQ=lapply(seq(.5,3,.5), function(x) {c(a,b,p,t,h,paste("betamultiplier",x,sep="="))})
LQP=lapply(seq(.5,3,.5), function(x) {c(a,b,t,h,paste("betamultiplier",x,sep="="))})
LQT=lapply(seq(.5,3,.5), function(x) {c(a,b,p,h,paste("betamultiplier",x,sep="="))})
LQH=lapply(seq(.5,3,.5), function(x) {c(a,b,p,t,paste("betamultiplier",x,sep="="))})
LQPT=lapply(seq(.5,3,.5), function(x) {c(a,b,h,paste("betamultiplier",x,sep="="))})
LQPH=lapply(seq(.5,3,.5), function(x) {c(a,b,t,paste("betamultiplier",x,sep="="))})
LQTH=lapply(seq(.5,3,.5), function(x) {c(a,b,t,h,paste("betamultiplier",x,sep="="))})
LQPTH=lapply(seq(.5,3,.5), function(x) {c(a,b,paste("betamultiplier",x,sep="="))})
args=c(LQ,LQP,LQT,LQH,LQPT,LQPH,LQTH,LQPTH)
# Create a list of the files where the results for each combination of settings will be saved (replace “your_path/your_file” by the desired path and file name):
LQs=lapply(seq(.5,3,.5),function(x) {paste("your_path/your_file_LQ",x,".csv",sep="")})
LQPs=lapply(seq(.5,3,.5),function(x) {paste("your_path/your_file_LQP",x,".csv",sep="")})
LQTs=lapply(seq(.5,3,.5),function(x) {paste("your_path/your_file_LQT",x,".csv",sep="")})
LQHs=lapply(seq(.5,3,.5),function(x) {paste("your_path/your_file_LQH",x,".csv",sep="")})
LQPTs=lapply(seq(.5,3,.5),function(x) {paste("your_path/your_file_LQPT",x,".csv",sep="")})
LQPHs=lapply(seq(.5,3,.5),function(x) {paste("your_path/your_file_LQPH",x,".csv",sep="")})
LQTHs=lapply(seq(.5,3,.5),function(x) {paste("your_path/your_file_LQTH",x,".csv",sep="")})
LQPTHs=lapply(seq(.5,3,.5),function(x) {paste("your_path/your_file_LQPTH",x,".csv",sep="")})
save=c(LQs,LQPs,LQTs,LQHs,LQPTs,LQPHs,LQTHs,LQPTHs)
# Load the following libraries required for parallel execution:
library(foreach)
library(doSNOW)
# Make a cluster corresponding to how many parallel versions of R should be run (avoid using every single core in your computer):
cl <- makeCluster(6, "SOCK")
registerDoSNOW(cl)
# Run the foreach loop in parallel:
maxent.runs <- foreach(i = 1:length(args), .packages = c("dismo", "rJava","ROCR")) %dopar% {
  mxt.nulltest(train, bg, n = 1000, proj=bg, c.proj=TRUE, args = args[[i]], group=group, bg.group=bg.group,save=save[[i]])}
stopCluster(cl)
# The results should have been saved to the specified files but they can also be viewed in R, or as a summary with the following lines of code:
summary=lapply(maxent.runs, function(x) x@summary)
names(summary)=save
