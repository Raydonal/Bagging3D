##########################################################################
### PACOTES


require(shapes)
require(Anthropometry)
require(cluster)
require(clue)
require(CircStats)
require(dendextend)
require(aricode)
library(parallel)



##############################################################################
# Algorithms

### HILL CLIMBING
hill_climbing_shapes = function(dados, num.clusters, max_iterations) {
  
  
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




##########################################################################
#Isotropy

cubos.generator.iso=function(n_sample, dispersion, semente){
  
  
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
  
  
  #library(Anthropometry)
  ##### ISOTROPY
  
  ## Cubes 8 Landmarks
  Ms_cube = array(c(0,0,0,0,10,10,10,10,0,10,0,10,10,10,0,0,10,10,0,0,10,0,10,0), dim = c(8, 3, 1))
  
  # Number of landmarks
  k_cube = dim(Ms_cube)[1]
  
  # Configuration matrix dimension
  vars_cube = k_cube * dim(Ms_cube)[2]
  
  # Dispersion: sigma^2 = (0.01,0.25,16,36,64,100)
  sigma_cube = dispersion^2
  
  # Covariance matrix
  Sigma_cube = diag(sigma_cube,vars_cube)
  
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

parallep.generator.iso=function(n_sample, dispersion, semente){
  
  
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
  
  
  #library(Anthropometry)
  
  ##### ISOTROPY
  
  ## Parallelepiped 8 Landmarks
  Ms_paral = array(  c(0,10,0,10,10,10,0,0,0,0,0,0,20,20,20,20,10,10,0,0,10,0,10,0), dim = c(8, 3, 1))
  
  # Number of landmarks
  k_paral = dim(Ms_paral)[1] 
  
  
  # Configuration matrix dimension
  vars_paral = k_paral * dim(Ms_paral)[2] 
  
  # Dispersion: sigma^2 = (0.01,0.25,16,36,64,100)
  sigma_paral = dispersion^2
  
  # Covariance matrix
  Sigma_paral = diag(sigma_paral,vars_paral)
  
  # Sample size
  n_paral=n_sample
  
  # Group of parallelepipeds:
  
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

clusterExport(cl, varlist = c("hill_climbing_shapes", "cl_bag.hill_climbing_shapes", "B"))

clusterEvalQ(cl, {
  library(shapes)
  library(Anthropometry)
  library(cluster)
  library(clue)
})







#####################################################################
set.seed(2021)

# Dispersion
sigma3=3


## Known clustering
verd3=c(rep(1,n_cube),rep(2,n_paral))


cubes3=list()
parallep3=list()
Cluster3=list()
rand3.hill.iso=list()
rand.bagg3.hill.iso=list()
FM3.hill.iso=list()
FM.bagg3.hill.iso=list()
time.without.bagg3=list()
time.with.bagg3=list()

for (c in 1:simulation){
  
  cubes3[[c]]=cubos.generator.iso(n_cube,sigma3,c)$cubos.rotate
  
  parallep3[[c]]=parallep.generator.iso(n_paral,sigma3,c)$parallep.rotate
  
  ## Set formed by the two groups
  
  Cluster3[[c]] = abind::abind(cubes3[[c]],parallep3[[c]])  
  
}

## Application of methods

inicio.parallelization3=Sys.time()
bagg3.hill.iso=parLapply(cl, Cluster3, function(x) cl_bag.hill_climbing_shapes(x, B=B, num.clusters=2, algorithm = "hill_climbing_shapes",  parameters = 
                                                                                 list(max_iterations = 300)))
fim.parallelization3=Sys.time()
time.parallelization3=fim.parallelization3-inicio.parallelization3

for (c in 1:simulation){
  
  # Rand Index without bagging
  
  rand3.hill.iso[[c]]=RI(bagg3.hill.iso[[c]]$original$best.clustering, verd3)
  
  # Rand Index with bagging
  
  rand.bagg3.hill.iso[[c]]=RI(max.col(as.matrix(bagg3.hill.iso[[c]]$bagg$.Data)), verd3)
  
  
  # FMI without bagging
  
  FM3.hill.iso[[c]]=FM_index_R(bagg3.hill.iso[[c]]$original$best.clustering, verd3)[1]
  
  # FMI with bagging
  
  FM.bagg3.hill.iso[[c]]=FM_index_R(max.col(as.matrix(bagg3.hill.iso[[c]]$bagg$.Data)), verd3)[1]
  
  # Times
  
  time.without.bagg3[[c]]=bagg3.hill.iso[[c]]$tempo1
  
  time.with.bagg3[[c]]=bagg3.hill.iso[[c]]$tempo2
  
}


export.hill.iso3=cbind(rand3.hill.iso,rand.bagg3.hill.iso, FM3.hill.iso, FM.bagg3.hill.iso,time.without.bagg3,time.with.bagg3,  time.parallelization3)

write.table(export.hill.iso3, file = "hill_N100_iso3_B100_simul50.csv", sep = ",", 
            na = "", quote = FALSE, row.names = F)





#####################################################################
set.seed(2021)

# Dispersion
sigma6=6

## Known clustering
verd6=c(rep(1,n_cube),rep(2,n_paral))

cubes6=list()
parallep6=list()
Cluster6=list()
rand6.hill.iso=list()
rand.bagg6.hill.iso=list()
FM6.hill.iso=list()
FM.bagg6.hill.iso=list()
time.without.bagg6=list()
time.with.bagg6=list()

for (b in 1:simulation){
  
  cubes6[[b]]=cubos.generator.iso(n_cube,sigma6,b)$cubos.rotate
  
  parallep6[[b]]=parallep.generator.iso(n_paral,sigma6,b)$parallep.rotate
  
  ## Set formed by the two groups
  
  Cluster6[[b]] = abind::abind(cubes6[[b]],parallep6[[b]])  
}

## Application of methods

inicio.parallelization6=Sys.time()
bagg6.hill.iso=parLapply(cl, Cluster6, function(x) cl_bag.hill_climbing_shapes(x, B=B, num.clusters=2, algorithm = "hill_climbing_shapes",  parameters = 
                                                                                 list(max_iterations = 300)))
fim.parallelization6=Sys.time()
time.parallelization6=fim.parallelization6-inicio.parallelization6

for (b in 1:simulation){
  
  # Rand Index without bagging
  
  rand6.hill.iso[[b]]=RI(bagg6.hill.iso[[b]]$original$best.clustering, verd6)
  
  # Rand Index with bagging
  
  rand.bagg6.hill.iso[[b]]=RI(max.col(as.matrix(bagg6.hill.iso[[b]]$bagg$.Data)), verd6)
  
  # FMI without bagging
  
  FM6.hill.iso[[b]]=FM_index_R(bagg6.hill.iso[[b]]$original$best.clustering, verd6)[1]
  
  # FMI with bagging
  
  FM.bagg6.hill.iso[[b]]=FM_index_R(max.col(as.matrix(bagg6.hill.iso[[b]]$bagg$.Data)), verd6)[1]
  
  # Times
  
  time.without.bagg6[[b]]=bagg6.hill.iso[[b]]$tempo1
  
  time.with.bagg6[[b]]=bagg6.hill.iso[[b]]$tempo2
  
}


export.hill.iso6=cbind(rand6.hill.iso,rand.bagg6.hill.iso, FM6.hill.iso, FM.bagg6.hill.iso, time.without.bagg6,time.with.bagg6,  time.parallelization6)

write.table(export.hill.iso6, file = "hill_N100_iso6_B100_simul50.csv", sep = ",", 
            na = "", quote = FALSE, row.names = F)




#####################################################################
set.seed(2021)

# Dispersion

sigma9=9


## Known clustering
verd9=c(rep(1,n_cube),rep(2,n_paral))

cubes9=list()
parallep9=list()
Cluster9=list()
rand9.hill.iso=list()
rand.bagg9.hill.iso=list()
FM9.hill.iso=list()
FM.bagg9.hill.iso=list()
time.without.bagg9=list()
time.with.bagg9=list()

for (a in 1:simulation){
  
  cubes9[[a]]=cubos.generator.iso(n_cube,sigma9,a)$cubos.rotate
  
  
  parallep9[[a]]=parallep.generator.iso(n_paral,sigma9,a)$parallep.rotate
  
  ## Set formed by the two groups
  
  Cluster9[[a]] = abind::abind(cubes9[[a]],parallep9[[a]])  
}  
  ## Application of methods
  
  inicio.parallelization9=Sys.time()
  bagg9.hill.iso=parLapply(cl, Cluster9, function(x) cl_bag.hill_climbing_shapes(x, B=B, num.clusters=2, algorithm = "hill_climbing_shapes",  parameters = 
                                                                                   list(max_iterations = 300)))
  fim.parallelization9=Sys.time()
  time.parallelization9=fim.parallelization9-inicio.parallelization9 
  
for (a in 1:simulation){
    
  # Rand Index without bagging
  
  rand9.hill.iso[[a]]=RI(bagg9.hill.iso[[a]]$original$best.clustering, verd9)
  
  # Rand Index with bagging
  
  rand.bagg9.hill.iso[[a]]=RI(max.col(as.matrix(bagg9.hill.iso[[a]]$bagg$.Data)), verd9)
  
  
  # FMI without bagging
  
  FM9.hill.iso[[a]]=FM_index_R(bagg9.hill.iso[[a]]$original$best.clustering, verd9)[1]
  
  # FMI with bagging
  
  FM.bagg9.hill.iso[[a]]=FM_index_R(max.col(as.matrix(bagg9.hill.iso[[a]]$bagg$.Data)), verd9)[1]
  
  #Times
  
  time.without.bagg9[[a]]=bagg9.hill.iso[[a]]$tempo1
  
  time.with.bagg9[[a]]=bagg9.hill.iso[[a]]$tempo2
  
}


export.hill.iso9=cbind(rand9.hill.iso, rand.bagg9.hill.iso, FM9.hill.iso, FM.bagg9.hill.iso, time.without.bagg9, time.with.bagg9, time.parallelization9)

write.table(export.hill.iso9, file = "hill_N100_iso9_B100_simul50.csv", sep = ",", 
            na = "", quote = FALSE, row.names = F)


stopCluster(cl)

