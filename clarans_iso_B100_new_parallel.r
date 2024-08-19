##########################################################################
### PACOTES


require(shapes)
require(Anthropometry)
require(cluster)
require(clue)
require(CircStats)
require(fastkmedoids)
require(aricode)
library(parallel)
require(dendextend)



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

clusterExport(cl, varlist = c("clarans_shapes", "cl_bag.claransshapes", "B"))

clusterEvalQ(cl, {
  library(shapes)
  library(Anthropometry)
  library(shapes)
  library(cluster)
  library(clue)
  library(CircStats)
  library(aricode)
  library(fastkmedoids)
})




#####################################################################
set.seed(2021)

sigma3=3


## Known clustering
verd3=c(rep(1,n_cube),rep(2,n_paral))


cubes3=list()
parallep3=list()
Cluster3=list()
rand3.clarans.iso=list()
rand.bagg3.clarans.iso=list()
FM3.clarans.iso=list()
FM.bagg3.clarans.iso=list()
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
  bagg3.clarans.iso=parLapply(cl, Cluster3, function(x) cl_bag.claransshapes(x, B=B, num.clusters=2,algorithm = "clarans_shapes",  parameters = 
                                                list(numlocal=2, maxneighbor=0.0125, seed=2021, type = "riem")))
  fim.parallelization3=Sys.time()
  time.parallelization3=fim.parallelization3-inicio.parallelization3

for (c in 1:simulation){
    
  
# Rand Index without bagging

  rand3.clarans.iso[[c]]=RI(bagg3.clarans.iso[[c]]$original$clustering, verd3)


# Rand Index with bagging

  rand.bagg3.clarans.iso[[c]]=RI(max.col(as.matrix(bagg3.clarans.iso[[c]]$bagg$.Data)), verd3)


# FMI without Bagging

  FM3.clarans.iso[[c]]=FM_index_R(bagg3.clarans.iso[[c]]$original$clustering, verd3)

  
# FMI with Bagging

  FM.bagg3.clarans.iso[[c]]=FM_index_R(max.col(as.matrix(bagg3.clarans.iso[[c]]$bagg$.Data)), verd3)

  
## Times
  
  time.without.bagg3[[c]]=bagg3.clarans.iso[[c]]$tempo1
  
  time.with.bagg3[[c]]=bagg3.clarans.iso[[c]]$tempo2
  
  
}


export.clarans.iso3=cbind(rand3.clarans.iso,rand.bagg3.clarans.iso, FM3.clarans.iso,  FM.bagg3.clarans.iso,
                          time.without.bagg3, time.with.bagg3, time.parallelization3)

write.table(export.clarans.iso3, file = "clarans_N100_iso3_B100_simul50.csv", sep = ",", 
            na = "", quote = FALSE, row.names = F)





#####################################################################
set.seed(2021)

sigma6=6


## Known clustering
verd6=c(rep(1,n_cube),rep(2,n_paral))


cubes6=list()
parallep6=list()
Cluster6=list()
rand6.clarans.iso=list()
rand.bagg6.clarans.iso=list()
FM6.clarans.iso=list()
FM.bagg6.clarans.iso=list()
time.without.bagg6=list()
time.with.bagg6=list()

for (c in 1:simulation){
  
  cubes6[[c]]=cubos.generator.iso(n_cube,sigma6,c)$cubos.rotate
  
  parallep6[[c]]=parallep.generator.iso(n_paral,sigma6,c)$parallep.rotate
  
  
  ## Set formed by the two groups
  
  Cluster6[[c]] = abind::abind(cubes6[[c]],parallep6[[c]])  
  
}

  ## Application of methods
  
  inicio.parallelization6=Sys.time()
  bagg6.clarans.iso=parLapply(cl, Cluster6, function(x) cl_bag.claransshapes(x, B=B, num.clusters=2,algorithm = "clarans_shapes",  parameters = 
                                                                               list(numlocal=2, maxneighbor=0.0125, seed=2021, type = "riem")))
  fim.parallelization6=Sys.time()
  time.parallelization6=fim.parallelization6-inicio.parallelization6

    
for (c in 1:simulation){
    
  # Rand Index without bagging
  
  rand6.clarans.iso[[c]]=RI(bagg6.clarans.iso[[c]]$original$clustering, verd6)
  
  
  # Rand Index with bagging
  
  rand.bagg6.clarans.iso[[c]]=RI(max.col(as.matrix(bagg6.clarans.iso[[c]]$bagg$.Data)), verd6)
  
  
  
  # FMI without Bagging
  
  FM6.clarans.iso[[c]]=FM_index_R(bagg6.clarans.iso[[c]]$original$clustering, verd6)
  
  
  # FMI with Bagging
  
  FM.bagg6.clarans.iso[[c]]=FM_index_R(max.col(as.matrix(bagg6.clarans.iso[[c]]$bagg$.Data)), verd6)
  
  
  ## Times
  
  time.without.bagg6[[c]]=bagg6.clarans.iso[[c]]$tempo1
  
  time.with.bagg6[[c]]=bagg6.clarans.iso[[c]]$tempo2
  
  
}


export.clarans.iso6=cbind(rand6.clarans.iso,rand.bagg6.clarans.iso, FM6.clarans.iso,  FM.bagg6.clarans.iso,
                          time.without.bagg6, time.with.bagg6, time.parallelization6)

write.table(export.clarans.iso6, file = "clarans_N100_iso6_B100_simul50.csv", sep = ",", 
            na = "", quote = FALSE, row.names = F)








#####################################################################
set.seed(2021)

sigma9=9


## Known clustering
verd9=c(rep(1,n_cube),rep(2,n_paral))


cubes9=list()
parallep9=list()
Cluster9=list()
rand9.clarans.iso=list()
rand.bagg9.clarans.iso=list()
FM9.clarans.iso=list()
FM.bagg9.clarans.iso=list()
time.without.bagg9=list()
time.with.bagg9=list()


for (c in 1:simulation){
  
  cubes9[[c]]=cubos.generator.iso(n_cube,sigma9,c)$cubos.rotate
  
  parallep9[[c]]=parallep.generator.iso(n_paral,sigma9,c)$parallep.rotate
  
  
  ## Set formed by the two groups
  
  Cluster9[[c]] = abind::abind(cubes9[[c]],parallep9[[c]])  
  
}

  ## Application of methods
  
  inicio.parallelization9=Sys.time()
  bagg9.clarans.iso=parLapply(cl, Cluster9, function(x) cl_bag.claransshapes(x, B=B, num.clusters=2,algorithm = "clarans_shapes",  parameters = 
                                                                               list(numlocal=2, maxneighbor=0.0125, seed=2021, type = "riem")))
  fim.parallelization9=Sys.time()
  time.parallelization9=fim.parallelization9-inicio.parallelization9


for (c in 1:simulation){
    
  # Rand Index without bagging
  
  rand9.clarans.iso[[c]]=RI(bagg9.clarans.iso[[c]]$original$clustering, verd9)
  
  
  # Rand Index with bagging
  
  rand.bagg9.clarans.iso[[c]]=RI(max.col(as.matrix(bagg9.clarans.iso[[c]]$bagg$.Data)), verd9)
  
  
  # FMI without Bagging
  
  FM9.clarans.iso[[c]]=FM_index_R(bagg9.clarans.iso[[c]]$original$clustering, verd9)
  
  
  # FMI with Bagging
  
  FM.bagg9.clarans.iso[[c]]=FM_index_R(max.col(as.matrix(bagg9.clarans.iso[[c]]$bagg$.Data)), verd9)
  
  
  ## Times
  
  time.without.bagg9[[c]]=bagg9.clarans.iso[[c]]$tempo1
  
  time.with.bagg9[[c]]=bagg9.clarans.iso[[c]]$tempo2
  
  
}


export.clarans.iso9=cbind(rand9.clarans.iso,rand.bagg9.clarans.iso, FM9.clarans.iso,  FM.bagg9.clarans.iso,
                          time.without.bagg9, time.with.bagg9, time.parallelization9)

write.table(export.clarans.iso9, file = "clarans_N100_iso9_B100_simul50.csv", sep = ",", 
            na = "", quote = FALSE, row.names = F)


stopCluster(cl)



