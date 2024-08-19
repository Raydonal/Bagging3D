
##############################################################################
require(shapes)
require(Anthropometry)
require(cluster)
require(clue)
require(CircStats)
require(fastkmedoids)
require(dendextend)
require(aricode)

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

### CLARANS 
clarans_shapes=function(dados, num.clusters, type=c("riem", "proc.full", "proc.partial"),
                        numlocal, maxneighbor, seed)
{
  proc.full = function (x, y, reflect = FALSE) 
  {
    require(shapes)
    out = sin(riemdist(x, y, reflect = reflect))
    out
  }
  proc.partial= function (x, y, reflect = FALSE) 
  {
    require(shapes)
    out = sqrt(2) * sqrt(abs(1 - cos(riemdist(x, y, reflect = reflect))))
    out
  }
  riem=function (x, y, reflect = FALSE) 
  {
    require(shapes)
    out = riemdist(x, y, reflect = reflect)
    out
  }  
  
  metric = function(type) {
    if(type == "proc.full") {
      return (proc.full)
    } else if(type == "proc.partial") {
      return (proc.partial)
    } else if(type == "riem") {
      return (riem)
    } else { 
      return (0)
    }
  }  
  
  distance=metric(type)
  matrixdist=function (x, distance_fcn) 
  {
    distance_from_idxs <- function(idxs) {
      i1 <- idxs[1]
      i2 <- idxs[2]
      distance_fcn(x[,,i1 ], x[, ,i2])
    }
    size <- dim(x)[3]
    d <- apply(utils::combn(size, 2), 2, distance_from_idxs)
    attr(d, "Size") <- size
    xnames <- rownames(x)
    if (!is.null(xnames)) {
      attr(d, "Labels") <- xnames
    }
    attr(d, "Diag") <- FALSE
    attr(d, "Upper") <- FALSE
    class(d) <- "dist"
    d
  }
  
  
  
  dissi=as.vector(matrixdist(dados, distance))
  N=dim(dados)[3]
  
  
  adapt=fastclarans(rdist=dissi, n=N, k=num.clusters, numlocal=numlocal, maxneighbor=maxneighbor, seed=seed)
  
  best.medoids=adapt@medoids
  clustering=adapt@assignment
  cases=dados[,,best.medoids]
  
  return(list(best.medoids=best.medoids, clustering=clustering, cases=cases))
  
}

cl_bag.claransshapes=function(dados, B, num.clusters, algorithm = "clarans_shapes",  parameters = 
                                list(numlocal, maxneighbor, seed, type = "riem")) 
{
  
  
  
  
  cl_boot.shape=function (dados, B, num.clusters, algorithm =  "clarans_shapes", 
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
      dados = rep.int(list(dados), B)
      eval(as.call(c(list(as.name("lapply"), dados, algorithm), 
                     if (!is.null(num.clusters)) list(num.clusters), parameters)))
    }
    else {
      replicate(B, expr = {
        algorithm = match.fun(algorithm)
        ind = sample(dim(dados)[3], replace = TRUE)
        train =  dados[,,ind]
        
        out = eval(as.call(c(list(algorithm, train), if (!is.null(num.clusters)) list(num.clusters), 
                             parameters)))
        clue:::as.cl_partition(cl_predict.shape(out, dados, "memberships"))
      }, simplify = FALSE)
    }
    cl_ensemble(list = clusterings)
  }
  
  
  n_of_classes.shape=function (x) {
    
    cl_class_ids.shape=function (x) as.cl_class_ids(x$clustering)
    
    length(unique(cl_class_ids.shape(x) ))
  }
  
  
  cl_membership.shape=function (x, num.clusters = n_of_classes.shape(x)) {
    cl_membership_from_class_ids.shape=function (x, num.clusters , meta = NULL) 
    {
      x = factor(x)
      n_of_objects = length(x)
      n_of_classes = nlevels(x)
      M = matrix(0, n_of_objects, num.clusters)
      M[cbind(seq_len(n_of_objects), as.numeric(x))] = 1
      M[is.na(x), ] = NA
      
      clue:::.make_cl_membership(M, n_of_classes, TRUE, meta)
    }
    
    cl_membership_from_class_ids.shape(x$clustering, num.clusters)
  }
  
  
  
  algorithm = match.fun(algorithm)
  ini1=Sys.time()
  reference = eval(as.call(c(list(algorithm, dados), if (!is.null(num.clusters)) list(num.clusters), 
                             parameters)))
  end1=Sys.time()
  time1=end1-ini1
  ini2=Sys.time()
  clusterings = cl_boot.shape(dados, B, num.clusters, algorithm, parameters, resample = TRUE)
  k1= n_of_classes.shape(reference)
  k2 =sapply(clusterings, n_of_classes)
  num.clusters=max(k1,k2)
  
  M_ref = cl_membership.shape(reference, num.clusters)
  M = matrix(0, NROW(M_ref), num.clusters)
  for (b in seq_len(B)) {
    mem = cl_membership(clusterings[[b]], num.clusters)
    ind = solve_LSAP(crossprod(M_ref, mem), maximum = TRUE)
    M = M + mem[, ind]
  }
  end2=Sys.time()
  time2=end2-ini2
  return(list(bagg=clue:::as.cl_partition(cl_membership(as.cl_membership(M/B), num.clusters)), original=reference, tempo1=time1,tempo2=time2))
}


### HILL CLIMBING 
hill_climbing_shapes <- function(dados, num.clusters, max_iterations) {
  
  
  value.criterion = function(dados, num.clusters, candidates) {
    matrix.dist <- matrix(0, dim(dados)[3], num.clusters)
    
    for (index in seq_along(candidates)) {
      matrix.dist[, index] <- apply(dados[, , 1:dim(dados)[3]], 3, riemdist, 
                                    dados[,,candidates[index]])
    }
    
    belong <- max.col(-matrix.dist)
    estimacao <- sum(apply(matrix.dist, 1, min))
    
    return(list(belong = belong, estimacao = estimacao, candidates = candidates))
  }
  
  
  n = dim(dados)[3]
  

  initial = sample(1:dim(dados)[3], num.clusters)  
  
  best_criterion = value.criterion(dados, num.clusters, candidates=initial)


  
  for (iteration in 1:max_iterations) {
    
    new_solution = initial
    random_index = sample(1:num.clusters, 1)
    new_neighbor = sample(1:n, 1)
    new_solution[random_index]=new_neighbor
    
    while(any(duplicated(new_solution))){
      random_index = sample(1:num.clusters, 1)
      new_neighbor = sample(1:n, 1)
      new_solution[random_index]=new_neighbor
    }
    
    new_criterion = value.criterion(dados, num.clusters, candidates=new_solution)
    
    if (new_criterion$estimacao < best_criterion$estimacao) {
      best_criterion=new_criterion
      
    }
  }
  
  
  
  return(list(best.center=best_criterion$candidates,
              best.criterion=best_criterion$estimacao , best.clustering=best_criterion$belong, cases=dados[,,best_criterion$candidates]))
}


cl_bag.hill_climbing_shapes=function(dados, B, num.clusters, algorithm = "hill_climbing_shapes",  parameters = 
                                       list(max_iterations = 300)) 
{
  
  
  cl_boot.shape=function (dados, B, num.clusters, algorithm =  "hill_climbing_shapes", 
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
      dados = rep.int(list(dados), B)
      eval(as.call(c(list(as.name("lapply"), dados, algorithm), 
                     if (!is.null(num.clusters)) list(num.clusters), parameters)))
    }
    else {
      replicate(B, expr = {
        algorithm = match.fun(algorithm)
        ind = sample(dim(dados)[3], replace = TRUE)
        train =  dados[,,ind]
        
        out = eval(as.call(c(list(algorithm, train), if (!is.null(num.clusters)) list(num.clusters), 
                             parameters)))
        clue:::as.cl_partition(cl_predict.shape(out, dados, "memberships"))
      }, simplify = FALSE)
    }
    cl_ensemble(list = clusterings)
  }
  
  
  n_of_classes.shape=function (x) {
    
    cl_class_ids.shape=function (x) as.cl_class_ids(x$best.clustering)
    
    length(unique(cl_class_ids.shape(x) ))
  }
  
  
  cl_membership.shape=function (x, num.clusters = n_of_classes.shape(x)) {
    cl_membership_from_class_ids.shape=function (x, num.clusters , meta = NULL) 
    {
      x = factor(x)
      n_of_objects = length(x)
      n_of_classes = nlevels(x)
      M = matrix(0, n_of_objects, num.clusters)
      M[cbind(seq_len(n_of_objects), as.numeric(x))] = 1
      M[is.na(x), ] = NA
      
      clue:::.make_cl_membership(M, n_of_classes, TRUE, meta)
    }
    
    cl_membership_from_class_ids.shape(x$best.clustering, num.clusters)
  }
  
  
  
  algorithm = match.fun(algorithm)
  ini1=Sys.time()
  reference = eval(as.call(c(list(algorithm, dados), if (!is.null(num.clusters)) list(num.clusters), 
                             parameters)))
  end1=Sys.time()
  time1=end1-ini1
  ini2=Sys.time()
  clusterings = cl_boot.shape(dados, B, num.clusters, algorithm, parameters, resample = TRUE)
  k1= n_of_classes.shape(reference)
  k2 =sapply(clusterings, n_of_classes)
  num.clusters=max(k1,k2)
  
  M_ref = cl_membership.shape(reference, num.clusters)
  M = matrix(0, NROW(M_ref), num.clusters)
  for (b in seq_len(B)) {
    mem = cl_membership(clusterings[[b]], num.clusters)
    ind = solve_LSAP(crossprod(M_ref, mem), maximum = TRUE)
    M = M + mem[, ind]
  }
  end2=Sys.time()
  time2=end2-ini2
  return(list(bagg=clue:::as.cl_partition(cl_membership(as.cl_membership(M/B), num.clusters)), original=reference, tempo1=time1,tempo2=time2))
}


###############################################################################
## Brains data set
data(brains)
Cluster=brains$x


## Known clustering
verd=brains$handed


set.seed(2023)

## Application of methods

## Kmeans
bagg1.kmeans=cl_bag.LloydShapes(Cluster, B=100, numClust=2, algorithm = "LloydShapes",  parameters = 
                                  list(algSteps=10,niter=10,stopCr = 1e-04,   simul=FALSE, verbose=FALSE)) 




## Clarans
bagg1.clarans=cl_bag.claransshapes(Cluster, B=100, num.clusters=2, algorithm = "clarans_shapes",  parameters = 
                                     list(numlocal=2, maxneighbor=0.0125, seed=2023, type = "riem")) 



## Hill
bagg1.hill=cl_bag.hill_climbing_shapes(Cluster, B=100, num.clusters=2, algorithm = "hill_climbing_shapes",  parameters = 
                                         list(max_iterations = 300)) 

### Rand

rand.kmeans=RI(bagg1.kmeans$original$asig, verd)
rand.bagg.kmeans=RI(max.col(as.matrix(bagg1.kmeans$bagg$.Data)), verd)


rand.clarans=RI(bagg1.clarans$original$clustering, verd)
rand.bagg.clarans=RI(max.col(as.matrix(bagg1.clarans$bagg$.Data)), verd)


rand.hill=RI(bagg1.hill$original$best.clustering, verd)
rand.bagg.hill=RI(max.col(as.matrix(bagg1.hill$bagg$.Data)), verd)



### FM
FM.kmeans=FM_index_R(bagg1.kmeans$original$asig, verd)[1]
FM.bagg.kmeans=FM_index_R(max.col(as.matrix(bagg1.kmeans$bagg$.Data)), verd)[1]


FM.clarans=FM_index_R(bagg1.clarans$original$clustering, verd)[1]
FM.bagg.clarans=FM_index_R(max.col(as.matrix(bagg1.clarans$bagg$.Data)), verd)[1]


FM.hill=FM_index_R(bagg1.hill$original$best.clustering, verd)[1]
FM.bagg.hill=FM_index_R(max.col(as.matrix(bagg1.hill$bagg$.Data)), verd)[1]



############

#Rand
rand.kmeans
rand.bagg.kmeans

rand.clarans
rand.bagg.clarans

rand.hill
rand.bagg.hill


# Relative gain Rand
(rand.bagg.kmeans-rand.kmeans)/(rand.kmeans)*100

(rand.bagg.clarans-rand.clarans)/(rand.clarans)*100

(rand.bagg.hill-rand.hill)/(rand.hill)*100


#FM

FM.kmeans
FM.bagg.kmeans

FM.clarans
FM.bagg.clarans

FM.hill
FM.bagg.hill


# Relative gain FM
(FM.bagg.kmeans-FM.kmeans)/(FM.kmeans)*100

(FM.bagg.clarans-FM.clarans)/(FM.clarans)*100

(FM.bagg.hill-FM.hill)/(FM.hill)*100



################## Plots

library(ggplot2)
library(gridExtra)
library(patchwork)

### PCA

out=procGPA(Cluster,scale=FALSE)
round(out$percent,2) #percentage of variability explained by the PC


# K-means
dad = data.frame(x = out$rawscores[,1], y = out$rawscores[,2])
dad$grupo_verdadeiro = verd
dad$grupo_metodo1 = bagg1.kmeans$original$asig
dad$grupo_metodo2 = max.col(as.matrix(bagg1.kmeans$bagg$.Data))

saveRDS(dad, file="brain.rds")


g1=ggplot(dad, aes(x = x, y = y, color = factor(grupo_verdadeiro), shape = factor(grupo_metodo1))) +
  geom_point(size = 4,alpha = 4) +
  scale_color_manual(values = c("#FF0000", "#469bd2")) +  # Cores para os grupos verdadeiros
  scale_shape_manual(values = c(1, 9)) +  # Formas para os grupos do método
  theme_bw() +
  labs(color = "True labels", shape = "Algorithm labels")+
  theme(panel.grid = element_line(color = "#dcdcdc", size = 0.1))


g2=ggplot(dad, aes(x = x, y = y, color = factor(grupo_verdadeiro), shape = factor(grupo_metodo2))) +
  geom_point(size = 4,alpha = 4) +
  scale_color_manual(values = c("#FF0000", "#469bd2")) +  # Cores para os grupos verdadeiros
  scale_shape_manual(values = c(1, 9)) +  # Formas para os grupos do método
  theme_bw() +
  labs(color = "True labels", shape = "Algorithm labels")+
  theme(panel.grid = element_line(color = "#dcdcdc", size = 0.1))


g1 = g1 + ggtitle("K-means")
g2 = g2 + ggtitle("K-means Bagging")



combined <- g1 + g2 & theme(legend.position = "right")

p1 <- combined + plot_layout(guides = "collect")


ggsave('Fig8a.pdf', p1, width = 10, height = 4)





# clarans
dad = data.frame(x = out$rawscores[,1], y = out$rawscores[,2])
dad$grupo_verdadeiro = verd
dad$grupo_metodo1 = bagg1.clarans$original$clustering
dad$grupo_metodo2 = max.col(as.matrix(bagg1.clarans$bagg$.Data))


g3=ggplot(dad, aes(x = x, y = y, color = factor(grupo_verdadeiro), shape = factor(grupo_metodo1))) +
  geom_point(size = 4,alpha = 4) +
  scale_color_manual(values = c("#FF0000", "#469bd2")) +  # Cores para os grupos verdadeiros
  scale_shape_manual(values = c(5, 13)) +  # Formas para os grupos do método
  theme_bw() +
  labs(color = "True labels", shape = "Algorithm labels")+
  theme(panel.grid = element_line(color = "#dcdcdc", size = 0.1))


g4=ggplot(dad, aes(x = x, y = y, color = factor(grupo_verdadeiro), shape = factor(grupo_metodo2))) +
  geom_point(size = 4,alpha = 4) +
  scale_color_manual(values = c("#FF0000", "#469bd2")) +  # Cores para os grupos verdadeiros
  scale_shape_manual(values = c(5, 13)) +  # Formas para os grupos do método
  theme_bw() +
  labs(color = "True labels", shape = "Algorithm labels")+
  theme(panel.grid = element_line(color = "#dcdcdc", size = 0.1))



g3 = g3 + ggtitle("CLARANS")
g4 = g4 + ggtitle("CLARANS Bagging")

combined <- g3 + g4 & theme(legend.position = "right")

p1 <- combined + plot_layout(guides = "collect")


ggsave('Fig8b.pdf', p1, width = 10, height = 4)

# Hill Climbing
dad = data.frame(x = out$rawscores[,1], y = out$rawscores[,2])
dad$grupo_verdadeiro = as.numeric(verd)
dad$grupo_metodo1 = bagg1.hill$original$best.clustering
dad$grupo_metodo2 = sort(max.col(as.matrix(bagg1.hill$bagg$.Data)))


g5=ggplot(dad, aes(x = x, y = y, color = factor(grupo_verdadeiro), shape = factor(grupo_metodo1))) +
  geom_point(size = 4,alpha = 4) +
  scale_color_manual(values = c("#FF0000", "#469bd2")) +  # Cores para os grupos verdadeiros
  scale_shape_manual(values = c(2, 8)) +  # Formas para os grupos do método
  theme_bw() +
  labs(color = "True labels", shape = "Algorithm labels")+
  theme(panel.grid = element_line(color = "#dcdcdc", size = 0.1))


g6=ggplot(dad, aes(x = x, y = y, color = factor(grupo_verdadeiro), shape = factor(grupo_metodo2))) +
  geom_point(size = 4,alpha = 4) +
  scale_color_manual(values = c("#FF0000", "#469bd2")) +  # Cores para os grupos verdadeiros
  scale_shape_manual(values = c(2, 8)) +  # Formas para os grupos do método
  theme_bw() +
  labs(color = "True labels", shape = "Algorithm labels")+
  theme(panel.grid = element_line(color = "#dcdcdc", size = 0.1))

g5 = g5 + ggtitle("Hill Climbing")
g6 = g6 + ggtitle("Hill Climbing Bagging")


combined <- g5 + g6 & theme(legend.position = "right")

p1 <- combined + plot_layout(guides = "collect")

ggsave('Fig8c.pdf', p1, width = 10, height = 4)
