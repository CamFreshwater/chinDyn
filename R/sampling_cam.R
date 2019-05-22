






#abundance
n1 <-2000
n2 <- 500
n3 <- 1200

real_dist <- c(n1,n2, n3)/sum(c(n1,n2, n3))
#fishing areas limits

fl <- c(10,20,30)





#mean distribution at time in a range from 1 to 30 of all possible positions

mupos2 <- seq(1,30,3)
mupos1 <- rep(c(5,15,25), each=5)
mupos3 <- rep(c(3,13,23), each=5)

# how patchy is the distribution of one stock? decrease sd to make it more patchy
sdpos <- 1.3

posX1 <- list()
posX2 <- list()
posX3  <- list()

Narea1  <- list()
Narea2  <- list()

sampl1  <- list()
sampl2  <- list()


#sampling days

samday <- 1:5

#loop over time
for(i in 1:10){

	#distribution of fish over space gradiaent

	posX1[[i]] <- pnorm(1:30 +.5, mupos1[i],sd=sdpos) - pnorm(1:30 -.5, 
	                                                          mupos1[i], sd=sdpos)
	posX2[[i]] <- pnorm(1:30 +.5, mupos2[i],sd=sdpos)-pnorm(1:30 -.5, 
	                                                        mupos2[i], sd=sdpos)
	posX3[[i]] <- pnorm(1:30 +.5, mupos3[i],sd=sdpos)-pnorm(1:30 -.5,
	                                                        mupos3[i], sd=sdpos)


	Narea1[[i]] <- c(sum(posX1[[i]][1:10] *n1),sum(posX2[[i]][1:10] *n2),sum(posX3[[i]][1:10] *n3))
	
	Narea2[[i]] <- c(sum(posX1[[i]][11:20] *n1),sum(posX2[[i]][11:20] *n2),sum(posX3[[i]][11:20] *n3))


	if(i %in%  samday ){

		sampl1[[i]] <- c(rmultinom(1,10, prob=Narea1[[i]]/sum(Narea1[[i]])))
		#sampl2[[i]] <- c(rmultinom(1,10, prob=Narea2[[i]]/sum(Narea2[[i]])))

		n1 <- n1 - sampl1[[i]][1]
		n2 <- n2 - sampl1[[i]][2]
		n3 <- n3 - sampl1[[i]][3]

	}

}




sampled_dist1 <- apply(matrix(unlist(sampl1),ncol=3, byrow=T),2,sum)/sum(apply(matrix(unlist(sampl1),ncol=3, byrow=T),2,sum))
#sampled_dist2 <- apply(matrix(unlist(sampl2),ncol=3, byrow=T),2,sum)/sum(apply(matrix(unlist(sampl2),ncol=3, byrow=T),2,sum))

real_dist







