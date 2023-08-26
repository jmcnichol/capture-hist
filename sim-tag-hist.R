# Tag history generation for the Jolly-Seber model with tag loss and retagging (including permanent tags)
# simulates capture and tag histories
# Jennifer McNichol
# 22 Nov 2021

# Here's an example. See below the function for some examples of test cases #
 n = 10
 nsample=3
 p = rep(1,nsample)
 phi = rep(1,nsample-1)
 lambda = matrix(rep(0.5,(nsample-1)*(nsample-1)),nrow=nsample-1)

 #test.taghist <- generate.taghist.retag(n,nsample,p,phi,lambda,0.5,T,123)
###################
 
 #'Notes:
 #'An individual that starts as single tagged must never have more than 1 tag.
 #'This function completely ignores loss on capture but it can be easily add in.


generate.taghist.retag <- function(n,nsample,p,phi,lambda,frac_double,brand=TRUE,seed){
  #'inputs: n = sample size
  #'        nsample = number of sample times 
  #'        p = capture prob
  #'        phi = survival prob
  #'        lambda = tag retention prob
  #'        frac_double = fraction of individuals double tagged
  #'        brand = TRUE for FALSE for branded or unbranded
  #'        seed = a seed
  #'output: taghist = matrix with n rows and 2*nsample columns of tag histories 
  #'        first = vector indicating the first time seen for each individual
  #'        last = vector of last time seen
  #'        retag = matrix indicating times when each individual was retagged
  #'        numtag = vector indicating the number of tags on each individual 
  
  #creation of vectors
  entrypoint=numeric(length=n)
  b=numeric(length=nsample) #probability of birth/immigration (sum to 1)
  bstar=numeric(length=nsample) # fraction entering the population
  first=numeric(length=n) #first entry time of individual i
  a=matrix(rep(0,n*nsample), nrow=n, ncol=nsample) #alive status of individual i at time j
  history=matrix(rep(0,n*nsample*2), nrow=n, ncol=(nsample*2))
  onetag=logical(n) #indicates if the individual only recieved one tag the first time it was tagged so we never give it two tags
  
  #initialization of parameters
  b[1:nsample]=1/nsample
  for(i in 1:nsample){
    bstar[i]=b[i]/sum(b[i:nsample])
  }
  # assuming no loss on capture for now
  loss_parm=0
  
  for(i in 1:n){
    set.seed(seed+i) # seed goes here to that everything stays the same for branded vs unbranded
    # Determine when the individual enters the population (before entrypoint[i])
    j=1
    alive=0
    while(alive==0){
      if(runif(1,0,1)<bstar[j]){
        alive=1
        entrypoint[i]=j
        break
      }#endif
      j=j+1
    }#endwhile
    # Determine when the individual was first captured
    l=entrypoint[i]
    a[i,entrypoint[i]]=1
    while((l<=nsample) & (alive==1)){
      if(runif(1,0,1)<p[l]){
        first[i]=l
        break
      }#endif
      if((l<nsample) & (first[i]==0)){
        if(runif(1,0,1)>phi[l]){
          alive=0
          #			a[i,l+1]=alive
        }#endif
      }#endif
      l=l+1
    }#endwhile
    
    # Determine tagging history for un branded
    # first tags
    if(first[i]>0 & (first[i]<nsample)){
      position=first[i]*2-1
      if(runif(1,0,1)<frac_double){ #is the individual double tagged?
        history[i,position]=1 
        history[i,position+1]=1
        tag1=1 # tag 1 status #tags retianed get 1
        tag2=1 #tag 2 status
      }else{
        history[i,position]=1 #single tag gets history'10'
        tag1=1
        tag2=0
        onetag[i] = TRUE
      }#endif
      
      # Determine if lost on capture
      loss=0
      if(runif(1,0,1)<loss_parm){ 
        loss=1
      }#endif
    
      alive=alive*(1-loss)
      
      for(k  in (first[i]+1):nsample){
        # Does the individual survive?
        if(runif(1,0,1)>phi[k-1]){
          alive=0
        }#endif
        # Does the individual retain its tags?	
        if(runif(1,0,1)>lambda[first[i],k-1] ){
          tag1=0 # tag is lost 
        }#endif
        if(runif(1,0,1)>lambda[first[i],k-1]){
          tag2=0
        }#endif
        
        if(tag1==2){ # if the retag from the last sample was retined it may still have the value 2 so reset it 
          tag1=1
        }#endif
        if(tag2==2){
          tag2=1
        }#endif
        
        # this part retagged if any tag was lost but we are only going to retag if all tags are lost. see next chunk.
        # retag the double tagged animlas making sure to only include animals who are captured
        #if(k<nsample){ # dont retag at the last sample time
        #  if (onetag[i]==FALSE) { # deal with single tags later along with the 00s since can only happen if theyre branded
        #    if(tag1==0 & tag2==1){
        #      tag1=2
        #    }#endif
        #    if(tag2==0 & tag1 ==1){
        #      tag2=2
        #    }#endif
        #  }#endif
        #}#endif  
        
        if (runif(1,0,1) <p[k]){ #is the individual captured?
          # retag branded individuals who lost all tags 
          if (brand == TRUE){
            if (onetag[i]==FALSE){
              if (tag1==0 & tag2==0 ){
                tag1=2
                tag2=2
              }#endif
            }else{
              if (tag1==0){
                tag1=2
                tag2=0
              }#endif
            }#endifelse
          }#endif
          
          if((alive==1) & ((tag1+tag2)>0)){ #animal is captured, alive and has at least one tag
            loss=0
            if(runif(1,0,1)<loss_parm){ #are they lost on capture?
              loss=1
            }#endif
            alive=alive*(1-loss)
            position=2*k-1
            history[i,position]=tag1
            history[i,position+1]=tag2
          }#endif
        }#endif
      }#endfor k 
    }#endif
    
    # First capture on last sampling time
    if(first[i]==nsample){
      position=first[i]*2-1
      if(runif(1,0,1)<frac_double){ #is the individual double tagged?
        history[i,position]=1 #double tagged gets history '11'
        history[i,position+1]=1
        tag1=1 # tag 1 status
        tag2=1 # tag 2 status
      }
      else{
        history[i,position]=1 #single tag gets history'10'
        tag1=1
      }#endif	
    }#endif
  }#endfor
  
  # Delete zero histories
  tag_history=matrix(0,nrow=n,ncol=2*nsample)
  index=0
  for(i in 1:length(first)){
    if(first[i]!=0){
      index=index+1
      tag_history[index,]=history[i,]
    }#endif
  }#endfor
  tag_history=tag_history[1:index,]

  ## processing stuff to put some of the outputs in the form we want ##
  #make vector to indicate the number of tags
  numtagfirst <- onetag
  numtagfirst[numtagfirst==TRUE] <- 1
  numtagfirst[numtagfirst==FALSE] <- 2
  
  
  #get last time seen
  #get every other column
  even.seq <- seq(from = 1, to = (nsample*2), by = 2)
  #add the two tag cols together for each sample time then find last non-zero time
  sum.tags <- sapply(even.seq,function(i) rowSums(tag_history[,i:(i+1)]))
  last <- max.col(sum.tags != 0, 'last')
  
  #get retag hist 
  retag_hist <- tag_history
  retag_hist[retag_hist==1] <- 0
  retag_hist[retag_hist==2] <- 1
  #only need one for each time captured
  retag_hist <- retag_hist[,even.seq]
  
  #dont want the deleted "never seen" histories to show up in first 
  first <- first[first != 0]
  
 
  
  # capture history 
  capt_hist <- sum.tags
  # 1 = captured, 0 = not captured 
  # using 5 since the sum of tags can never be 5
  capt_hist[capt_hist == 0] <- 5
  capt_hist[capt_hist != 5] <- 1
  capt_hist[capt_hist == 5] <- 0
  
  
  #make all the tag hist 1
  tag_history[tag_history==2] <- 1
  
  return(list(taghist = tag_history, caphist = capt_hist ,first = first, last = last, numtag = numtagfirst, retaghist = retag_hist))
}



###### Testing #####
 
#I did a bit of testing to make sure the data is generated as expected, here are some examples:
#Tags are always retained.
n = 10
nsample=3
p = rep(1,nsample)
phi = rep(1,nsample-1)
lambda = matrix(rep(1,(nsample-1)*(nsample-1)),nrow=nsample-1)
generate.taghist.retag(n,nsample,p,phi,lambda,0.5,T,123)




#If lambda = 0 and p = phi = 1 and animals are branded we get retagging of all tags after the first tagging.
n = 10
nsample=3
p = rep(1,nsample)
phi = rep(1,nsample-1)
lambda = matrix(rep(0,(nsample-1)*(nsample-1)),nrow=nsample-1)
generate.taghist.retag(n,nsample,p,phi,lambda,0.5,T,123)


#And if theyre not branded... we have to assume 00 means unseen.
n = 10
nsample=3
p = rep(1,nsample)
phi = rep(1,nsample-1)
lambda = matrix(rep(0,(nsample-1)*(nsample-1)),nrow=nsample-1)
generate.taghist.retag(n,nsample,p,phi,lambda,0.5,F,123)


#Check to make sure we didnt just replace the 00 with retag (ie. the history goes back though the 
#loop based on the retag it is determined if the new tags stay or are lost) *see row 7 in the branded vs unbranded*
n = 10
nsample=4
p = rep(1,nsample)
phi = rep(1,nsample-1)
lambda = matrix(rep(0.5,(nsample-1)*(nsample-1)),nrow=nsample-1)
generate.taghist.retag(n,nsample,p,phi,lambda,0.5,F,98765) # 00 in sample times 3 and 4 
generate.taghist.retag(n,nsample,p,phi,lambda,0.5,T,98765) # retag in smaple time 3 and the retained in 4




#Branded and p < 1, make sure we still have cases of unseen in branded animals.
n = 10
nsample=3
p = rep(0.5,nsample)
phi = rep(1,nsample-1)
lambda = matrix(rep(0,(nsample-1)*(nsample-1)),nrow=nsample-1)
generate.taghist.retag(n,nsample,p,phi,lambda,0.5,T,123)



#Here we can see that we have a case where a tagged individual lost all tags so was retagged then was not seen again.
n = 10
nsample=4
p = rep(0.5,nsample)
phi = rep(0.5,nsample-1)
lambda = matrix(rep(0,(nsample-1)*(nsample-1)),nrow=nsample-1)
generate.taghist.retag(n,nsample,p,phi,lambda,0.5,T,123)
