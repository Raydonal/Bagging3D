##########################################################################
### PACOTES


require(shapes)
require(Anthropometry)
require(cluster)
library(clue)
require(CircStats)
require(aricode)
require(parallel)
require(dendextend)


##############################################################################
# Algorithms

### K-means 

LloydShapes=function (array3D, numClust, algSteps = 10, niter = 10, stopCr = 1e-04, 
                      simul, verbose) 
{
  
  
  
  time_iter <- list()
  comp_time <- c()
  list_asig_step <- list()
  list_asig <- list()
  vect_all_rate <- c()
  initials <- list()
  ll <- 1:numClust
  dist <- matrix(0, dim(array3D)[3], numClust)
  if (verbose) {
    print(Sys.time())
  }
  time_ini <- Sys.time()
  vopt <- 1e+08
  for (iter in 1:niter) {
    obj <- list()
    meanshapes <- 0
    meanshapes_aux <- 0
    asig <- 0
    mean_sh <- list()
    n <- dim(array3D)[3]
    if (verbose) {
      cat("New iteration:")
      print(iter)
      cat("Optimal value with which this iteration starts:")
      print(vopt)
    }
    initials[[iter]] <- sample(1:n, numClust, replace = FALSE)
    if (verbose) {
      cat("Initial values of this iteration:")
      print(initials[[iter]])
    }
    meanshapes <- array3D[, , initials[[iter]]]
    meanshapes_aux <- array3D[, , initials[[iter]]]
    for (step in 1:algSteps) {
      for (h in 1:numClust) {
        dist[, h] = apply(array3D[, , 1:n], 3, riemdist, 
                          y = meanshapes[, , h])
      }
      asig = max.col(-dist)
      for (h in 1:numClust) {
        if (table(asig == h)[2] == 1) {
          meanshapes[, , h] = array3D[, , asig == h]
          mean_sh[[step]] <- meanshapes
        }
        else {
          meanshapes[, , h] = procGPA(array3D[, , asig == 
                                                h], distances = TRUE, pcaoutput = TRUE)$mshape
          mean_sh[[step]] <- meanshapes
        }
      }
      obj[[step]] <- c(0)
      for (l in 1:n) {
        obj[[step]] <- obj[[step]] + dist[l, asig[l]]^2
      }
      obj[[step]] <- obj[[step]]/n
      list_asig_step[[step]] <- asig
      if (verbose) {
        paste(cat("Clustering of the Nstep", step, 
                  ":\n"))
        print(table(list_asig_step[[step]]))
      }
      if (verbose) {
        if (iter <= 10) {
          paste(cat("Objective function of the Nstep", 
                    step))
          print(obj[[step]])
        }
      }
      if (step > 1) {
        aux <- obj[[step]]
        aux1 <- obj[[step - 1]]
        if (((aux1 - aux)/aux1) < stopCr) {
          break
        }
      }
    }
    obj1 <- 0
    for (l in 1:n) {
      obj1 <- obj1 + dist[l, asig[l]]^2
    }
    obj1 <- obj1/n
    if (obj1 > min(unlist(obj))) {
      if (min(unlist(obj)) < vopt) {
        vopt <- min(unlist(obj))
        if (verbose) {
          cat("optimal")
          print(vopt)
        }
        optim_obj <- which.min(unlist(obj))
        copt <- mean_sh[[optim_obj]]
        asig_opt <- list_asig_step[[optim_obj]]
      }
    }
    else if (obj1 < vopt) {
      vopt <- obj1
      if (verbose) {
        cat("optimal")
        print(vopt)
      }
      optim_obj <- which.min(unlist(obj))
      copt <- mean_sh[[optim_obj]]
      asig_opt <- list_asig_step[[optim_obj]]
    }
    time_iter[[iter]] <- Sys.time()
    if (iter == 1) {
      comp_time[1] <- difftime(time_iter[[iter]], time_ini, 
                               units = "mins")
      if (verbose) {
        cat("Computational time of this iteration: \n")
        print(time_iter[[iter]] - time_ini)
      }
    }
    else {
      comp_time[iter] <- difftime(time_iter[[iter]], time_iter[[iter - 
                                                                  1]], units = "mins")
      if (verbose) {
        cat("Computational time of this iteration: \n")
        print(time_iter[[iter]] - time_iter[[iter - 1]])
      }
    }
    if (verbose) {
      cat("Optimal clustering of this iteration: \n")
    }
    optim_obj <- which.min(unlist(obj))
    list_asig[[iter]] <- list_asig_step[[optim_obj]]
    if (verbose) {
      print(table(list_asig[[iter]]))
    }
    if (simul) {
      as1 <- table(list_asig[[iter]][1:(n/2)])
      as2 <- table(list_asig[[iter]][seq(n/2 + 1, n)])
      if (max(as1) != n/2 & max(as2) != n/2) {
        suma <- min(as1) + min(as2)
        all_rate <- 1 - suma/n
      }
      else if ((max(as1) == n/2 & max(as2) != n/2) || (max(as1) != 
                                                       n/2 & max(as2) == n/2)) {
        minim <- min(min(as1), min(as2))
        all_rate <- 1 - minim/n
      }
      else if (max(as1) == n/2 & max(as2) == n/2) {
        all_rate <- 1
      }
      vect_all_rate[iter] <- all_rate
      if (verbose) {
        cat("Optimal allocation rate in this iteration:")
        print(all_rate)
      }
    }
  }
  if (simul) {
    dimnames(copt) <- NULL
    return(list(asig = asig_opt, cases = copt, vopt = vopt, 
                compTime = comp_time, AllRate = vect_all_rate, initials = initials))
  }
  else {
    return(list(asig = asig_opt, cases = copt, vopt = vopt, 
                initials = initials))
  }
}


cl_bag.LloydShapes=function(array3D, B, numClust, algorithm = "LloydShapes",  parameters = 
                              list(algSteps=10,niter=10,stopCr = 1e-04,   simul=FALSE, verbose=FALSE)) 
{
  
  
  
  
  cl_boot.shape=function (array3D, B, numClust, algorithm =  "LloydShapes", 
                          parameters = list(), resample = TRUE) 
  {
    
    
    cl_predict.shape=function (out, newdata = NULL, type = c("class_ids", 
                                                             "memberships")) 
    {
      
      rxdist.shape=function (A, B) {
        
        resultado = data.frame(matrix(0, ncol = dim(B)[3], nrow = dim(A)[3]))
        for(index_point in seq(dim(A)[3])){
          
          #distances = c()
          for(index in seq(dim(B)[3])){
            resultado[index_point,index] = riemdist(A[,,index_point], B[,,index])
            
          }
        }
        resultado
        
        
      }
      
      as_cl_class_ids_or_membership.shape=function (x, type = c("class_ids", "memberships")) 
      {
        type <- match.arg(type)
        if (type == "class_ids") {
          if (is.matrix(x)) {
            as.cl_class_ids(.structure(max.col(x), names = rownames(x)))
          }
          else as.cl_class_ids(x)
        }
        else as.cl_membership(x)
      }
      
      
      d = rxdist.shape(newdata, out$cases)
      as_cl_class_ids_or_membership.shape(max.col(-d), type)
    }
    
    
    clusterings = if (!resample) {
      array3D = rep.int(list(array3D), B)
      eval(as.call(c(list(as.name("lapply"), array3D, algorithm), 
                     if (!is.null(numClust)) list(numClust), parameters)))
    }
    else {
      replicate(B, expr = {
        algorithm = match.fun(algorithm)
        ind = sample(dim(array3D)[3], replace = TRUE)
        train =  array3D[,,ind]
        
        out = eval(as.call(c(list(algorithm, train), if (!is.null(numClust)) list(numClust), 
                             parameters)))
        clue:::as.cl_partition(cl_predict.shape(out, array3D, "memberships"))
      }, simplify = FALSE)
    }
    cl_ensemble(list = clusterings)
  }
  
  
  n_of_classes.shape=function (x) {
    
    cl_class_ids.shape=function (x) as.cl_class_ids(x$asig)
    
    length(unique(cl_class_ids.shape(x) ))
  }
  
  
  cl_membership.shape=function (x, k = n_of_classes.shape(x)) {
    cl_membership_from_class_ids.shape=function (x, k , meta = NULL) 
    {
      x = factor(x)
      n_of_objects = length(x)
      n_of_classes = nlevels(x)
      M = matrix(0, n_of_objects, k)
      M[cbind(seq_len(n_of_objects), as.numeric(x))] = 1
      M[is.na(x), ] = NA
      
      clue:::.make_cl_membership(M, n_of_classes, TRUE, meta)
    }
    
    cl_membership_from_class_ids.shape(x$asig, k)
  }
  
  
  algorithm = match.fun(algorithm)
  ini1=Sys.time()
  reference = eval(as.call(c(list(algorithm, array3D), if (!is.null(numClust)) list(numClust), 
                             parameters)))
  end1=Sys.time()
  time1=end1-ini1
  ini2=Sys.time()
  clusterings = cl_boot.shape(array3D, B, numClust, algorithm, parameters, resample = TRUE)
  k1= n_of_classes.shape(reference)
  k2 =sapply(clusterings, n_of_classes)
  k=max(k1,k2)
  
  M_ref = cl_membership.shape(reference, k)
  M = matrix(0, NROW(M_ref), k)
  for (b in seq_len(B)) {
    mem = cl_membership(clusterings[[b]], k)
    ind = solve_LSAP(crossprod(M_ref, mem), maximum = TRUE)
    M = M + mem[, ind]
  }
  end2=Sys.time()
  time2=end2-ini2
  return(list(bagg=clue:::as.cl_partition(cl_membership(as.cl_membership(M/B), k)), original=reference, tempo1=time1,tempo2=time2))
}



##########################################################################
#Anisotropy

cubos.generator.aniso=function(n_sample, dispersion, semente){
  
  ## seed to ensure the reproducibility of the experiment
  set.seed(semente)  
  
  
  param1=seq(-pi,pi)
  param2=sample(param1,1)
  param3=sample(1:2049^3,1)
  
  
  ## Function to rotate
  rotate=function(n , mean, k, dados)
  {
    angle=rvm(n=1, mean, k)
    
    x=c()
    y=c()
    z=c()
    x1=c()
    y1=c()
    z1=c()
    dados.rot=c()
    
    for(i in seq(dim(dados)[3])){
      x[i]=list(dados[,1,i])
      y[i]=list(dados[,2,i])
      z[i]=list(dados[,3,i])
      
      x1[i]=list( x[[i]]*cos(angle) - y[[i]]*sin(angle) )
      y1[i]=list( x[[i]]*sin(angle) + y[[i]]*cos(angle) )
      z1[i]=list( z[[i]] )
      
      
      dados.rot[i]=list(cbind(x1[[i]],y1[[i]],z1[[i]]))
      
      
    }
    return(dados.rot)
  }
  
  library(Anthropometry)
  
  ##### ANISOTROPY
  
  ## Cubes 8 Landmarks
  Ms_cube = cube8landm
  
  # Number of landmarks
  k_cube = dim(Ms_cube)[1]
  
  # Configuration matrix dimension
  vars_cube = k_cube * dim(Ms_cube)[2]
  
  # Dispersion: sigma^2 = (0.01,0.25,16,36,64,100)
  sigma_cube = dispersion^2
  
  # Covariance matrix
  unsk=as.vector(rep(1,k_cube))
  unsk.tran=t(unsk)
  gamma=4
  unsm=as.vector(rep(1,3))
  unsm.tran=t(unsm)
  
  Sigma_cube = kronecker(sigma_cube * ( (unsk%*%unsk.tran) + ((gamma-1)*diag(k_cube)))
                         , ((unsm%*%unsm.tran)+ (gamma-1)*diag(3) )/(gamma^2))
  
  
  
  
  
  # Sample size
  n_cube = n_sample
  
  ## Group of cubes:
  cubos.simul = mvtnorm::rmvt(n_cube,Sigma_cube,df=99)[,c(1 : k_cube * dim(Ms_cube)[2] 
                                                          - 2, 1 : k_cube * dim(Ms_cube)[2] - 1, 1 : k_cube * dim(Ms_cube)[2])]
  cubos.Simul = as.vector(Ms_cube) + t(cubos.simul) 
  
  
  # Configuration matrices for each cube in the group
  cluster = array(cubos.Simul, dim = c(k_cube, dim(Ms_cube)[2], n_cube))
  
  # Rotating the cubes
  rotate.cubes=rotate(1 , mean=param3, k=param2, cluster)
  
  # Groups of rotated cubes
  cl=array(unlist(rotate.cubes), dim = c(k_cube,3,n_cube))
  
  return(list(cubos=cluster, cubos.rotate=cl))
}

parallep.generator.aniso=function(n_sample, dispersion, semente){
  
  ## seed to ensure the reproducibility of the experiment
  set.seed(semente)  
  
  
  param1=seq(-pi,pi)
  param2=sample(param1,1)
  param3=sample(1:2049^3,1)
  
  
  ## Function to rotate
  rotate=function(n , mean, k, dados)
  {
    angle=rvm(n=1, mean, k)
    
    x=c()
    y=c()
    z=c()
    x1=c()
    y1=c()
    z1=c()
    dados.rot=c()
    
    for(i in seq(dim(dados)[3])){
      x[i]=list(dados[,1,i])
      y[i]=list(dados[,2,i])
      z[i]=list(dados[,3,i])
      
      x1[i]=list( x[[i]]*cos(angle) - y[[i]]*sin(angle) )
      y1[i]=list( x[[i]]*sin(angle) + y[[i]]*cos(angle) )
      z1[i]=list( z[[i]] )
      
      
      dados.rot[i]=list(cbind(x1[[i]],y1[[i]],z1[[i]]))
      
      
    }
    return(dados.rot)
  }
  
  library(Anthropometry)
  
  ## Parallelepiped 8 Landmarks
  Ms_paral = parallelep8landm
  
  # Number of landmarks
  k_paral = dim(Ms_paral)[1] 
  
  
  # Configuration matrix dimension
  vars_paral = k_paral * dim(Ms_paral)[2] 
  
  # Dispersion: sigma^2 = (0.01,0.25,16,36,64,100)
  sigma_paral = dispersion^2
  
  # Covariance matrix
  unsk=as.vector(rep(1,k_paral))
  unsk.tran=t(unsk)
  gamma=4
  unsm=as.vector(rep(1,3))
  unsm.tran=t(unsm)
  
  
  Sigma_paral = kronecker(  sigma_paral * ( (unsk%*%unsk.tran) + ((gamma-1)*diag(k_paral)) )
                            , 
                            ( (unsm%*%unsm.tran)+ (gamma-1)*diag(3) )/(gamma^2) )
  
  
  
  
  
  # Sample size
  n_paral = n_sample
  
  ## Group of parallelepipeds:
  simu2_paral = mvtnorm::rmvt(n_paral, Sigma_paral, df = 99)[,c(1 : k_paral * 
                                                                  dim(Ms_paral)[2] - 2, 1 : k_paral * dim(Ms_paral)[2] - 
                                                                  1, 1 : k_paral * dim(Ms_paral)[2])]
  Simu2_paral = as.vector(Ms_paral) + t(simu2_paral) 
  
  
  
  # Configuration matrices of each parallelepiped in the group
  cluster = array(Simu2_paral, dim = c(k_paral, dim(Ms_paral)[2], n_paral))
  
  # Rotating the parallelepipeds
  rotate.para=rotate(1 , mean=param3, k=param2, cluster)
  
  # Group of rotated parallelepipeds
  cl=array(unlist(rotate.para), dim = c(k_paral,3,n_paral))
  
  return(list(parallep=cluster, parallep.rotate=cl))
}



######################################################################

## Parameters for simulation

simulation=50
n_cube=50
n_paral=50
B=100

###################################################
cl <- makeCluster(3)

clusterExport(cl, varlist = c("cl_bag.LloydShapes", "LloydShapes", "B"))

clusterEvalQ(cl, {
  library(shapes)
  library(Anthropometry)
  library(shapes)
  library(cluster)
  library(clue)
})



###################################################################################

cubes3=list()
parallep3=list()
Cluster3=list()
rand3.kmeans.aniso=list()
rand.bagg3.kmeans.aniso=list()
FM3.kmeans.aniso=list()
FM.bagg3.kmeans.aniso=list()
time.without.bagg3=list()
time.with.bagg3=list()

## Known clustering
verd3=c(rep(1,n_cube),rep(2,n_paral))


## Dispersion
sigma3=3

set.seed(2021)


for (c in 1:simulation) {
  
  cubes3[[c]]=cubos.generator.aniso(n_cube,sigma3, c)$cubos.rotate
  
  
  parallep3[[c]]=parallep.generator.aniso(n_paral,sigma3,c)$parallep.rotate
  
  
  ## Set formed by the two groups
  
  
  Cluster3[[c]] = abind::abind(cubes3[[c]],parallep3[[c]])  
  
}

## Application of methods

inicio.parallelization3=Sys.time()
bagg3.kmeans.aniso=parLapply(cl, Cluster3, function(x) cl_bag.LloydShapes(x, B=B, numClust=2, algorithm = "LloydShapes",  parameters = 
                                                                          list(algSteps=10,niter=10,stopCr = 1e-04,simul=FALSE, verbose=FALSE)))
fim.parallelization3=Sys.time()
time.parallelization3=fim.parallelization3-inicio.parallelization3




for (c in 1:simulation) {
  
  # Rand Index without bagging
  rand3.kmeans.aniso[[c]]=RI(bagg3.kmeans.aniso[[c]]$original$asig, verd3)
  
  
  
  # Rand Index with bagging
  rand.bagg3.kmeans.aniso[[c]]=RI(max.col(as.matrix(bagg3.kmeans.aniso[[c]]$bagg$.Data)), verd3)
  
  
  # FMI without bagging
  FM3.kmeans.aniso[[c]]=FM_index_R(bagg3.kmeans.aniso[[c]]$original$asig, verd3)
  
  
  # FMI with bagging
  FM.bagg3.kmeans.aniso[[c]]=FM_index_R(max.col(as.matrix(bagg3.kmeans.aniso[[c]]$bagg$.Data)), verd3)
  
  # Times
  time.without.bagg3[[c]]=bagg3.kmeans.aniso[[c]]$tempo1
  
  time.with.bagg3[[c]]=bagg3.kmeans.aniso[[c]]$tempo2
  
  
}

export.kmeans.aniso3=cbind(rand3.kmeans.aniso,rand.bagg3.kmeans.aniso, FM3.kmeans.aniso,  FM.bagg3.kmeans.aniso, time.without.bagg3, time.with.bagg3, time.parallelization3
)

write.table(export.kmeans.aniso3, file = "kmeans_N100_aniso3_B100_simul50.csv", sep = ",", 
            na = "", quote = FALSE, row.names = F)






###################################################################################

cubes6=list()
parallep6=list()
Cluster6=list()
rand6.kmeans.aniso=list()
rand.bagg6.kmeans.aniso=list()
FM6.kmeans.aniso=list()
FM.bagg6.kmeans.aniso=list()
time.without.bagg6=list()
time.with.bagg6=list()

## Known clustering
verd6=c(rep(1,n_cube),rep(2,n_paral))


## Dispersion
sigma6=6

set.seed(2021)

for (c in 1:simulation) {
  
  cubes6[[c]]=cubos.generator.aniso(n_cube,sigma6, c)$cubos.rotate
  
  
  parallep6[[c]]=parallep.generator.aniso(n_paral,sigma6,c)$parallep.rotate
  
  
  
  ## Set formed by the two groups
  
  
  Cluster6[[c]] = abind::abind(cubes6[[c]],parallep6[[c]])  
  
}

inicio.parallelization6=Sys.time()

## Application of methods
bagg6.kmeans.aniso=parLapply(cl, Cluster6, function(x) cl_bag.LloydShapes(x, B=B, numClust=2, algorithm = "LloydShapes",  parameters = 
                                                                          list(algSteps=10,niter=10,stopCr = 1e-04,simul=FALSE, verbose=FALSE)))
fim.parallelization6=Sys.time()
time.parallelization6=fim.parallelization6-inicio.parallelization6



for (c in 1:simulation) {
  
  # Rand Index without bagging
  rand6.kmeans.aniso[[c]]=RI(bagg6.kmeans.aniso[[c]]$original$asig, verd6)
  
  
  # Rand Index with bagging
  rand.bagg6.kmeans.aniso[[c]]=RI(max.col(as.matrix(bagg6.kmeans.aniso[[c]]$bagg$.Data)), verd6)
  
  
  # FMI without bagging
  FM6.kmeans.aniso[[c]]=FM_index_R(bagg6.kmeans.aniso[[c]]$original$asig, verd6)
  
  
  
  # FMI with bagging
  FM.bagg6.kmeans.aniso[[c]]=FM_index_R(max.col(as.matrix(bagg6.kmeans.aniso[[c]]$bagg$.Data)), verd6)
  
  
  
  time.without.bagg6[[c]]=bagg6.kmeans.aniso[[c]]$tempo1
  
  time.with.bagg6[[c]]=bagg6.kmeans.aniso[[c]]$tempo2
  
  
}

export.kmeans.aniso6=cbind(rand6.kmeans.aniso,rand.bagg6.kmeans.aniso, FM6.kmeans.aniso,  FM.bagg6.kmeans.aniso, time.without.bagg6, time.with.bagg6, time.parallelization6)

write.table(export.kmeans.aniso6, file = "kmeans_N100_aniso6_B100_simul50.csv", sep = ",", 
            na = "", quote = FALSE, row.names = F)




###################################################################################

cubes9=list()
parallep9=list()
Cluster9=list()
rand9.kmeans.aniso=list()
rand.bagg9.kmeans.aniso=list()
FM9.kmeans.aniso=list()
FM.bagg9.kmeans.aniso=list()
time.without.bagg9=list()
time.with.bagg9=list()

## Known clustering
verd9=c(rep(1,n_cube),rep(2,n_paral))


## Dispersion
sigma9=9

set.seed(2021)


for (c in 1:simulation) {
  
  cubes9[[c]]=cubos.generator.aniso(n_cube,sigma9, c)$cubos.rotate
  
  
  parallep9[[c]]=parallep.generator.aniso(n_paral,sigma9,c)$parallep.rotate
  
  
  
  ## Set formed by the two groups
  
  
  Cluster9[[c]] = abind::abind(cubes9[[c]],parallep9[[c]])  
  
}  

## Application of methods

inicio.parallelization9=Sys.time()
bagg9.kmeans.aniso=parLapply(cl, Cluster9, function(x) cl_bag.LloydShapes(x, B=B, numClust=2, algorithm = "LloydShapes",  parameters = 
                                                                          list(algSteps=10,niter=10,stopCr = 1e-04,simul=FALSE, verbose=FALSE)))
fim.parallelization9=Sys.time()
time.parallelization9=fim.parallelization9-inicio.parallelization9


for (c in 1:simulation) {
  
  # Rand Index without bagging
  rand9.kmeans.aniso[[c]]=RI(bagg9.kmeans.aniso[[c]]$original$asig, verd9)
  
  
  # Rand Index with bagging
  rand.bagg9.kmeans.aniso[[c]]=RI(max.col(as.matrix(bagg9.kmeans.aniso[[c]]$bagg$.Data)), verd9)
  
  
  
  # FMI without bagging
  FM9.kmeans.aniso[[c]]=FM_index_R(bagg9.kmeans.aniso[[c]]$original$asig, verd9)
  
  
  # FMI with bagging
  FM.bagg9.kmeans.aniso[[c]]=FM_index_R(max.col(as.matrix(bagg9.kmeans.aniso[[c]]$bagg$.Data)), verd9)
  
  
  time.without.bagg9[[c]]=bagg9.kmeans.aniso[[c]]$tempo1
  
  time.with.bagg9[[c]]=bagg9.kmeans.aniso[[c]]$tempo2
  
  
}

export.kmeans.aniso9=cbind(rand9.kmeans.aniso,rand.bagg9.kmeans.aniso, FM9.kmeans.aniso,  FM.bagg9.kmeans.aniso, time.without.bagg9, time.with.bagg9,time.parallelization9)

write.table(export.kmeans.aniso9, file = "kmeans_N100_aniso9_B100_simul50.csv", sep = ",", 
            na = "", quote = FALSE, row.names = F)


stopCluster(cl)


