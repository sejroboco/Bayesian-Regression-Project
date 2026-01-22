# TP0
# ENSAI
# Library utiles  
library(mvtnorm) # pour simuler des vecteurs gaussiens
library(mnormt)
library(invgamma) # pour simuler une inverse gamma

####################################################### 
# Simulations des données que l'on étudiera 
simuvarexpl <- matrix(rep(0,300*200),ncol=300,nrow=200)
for (i in 1 :300) simuvarexpl[,i] <- runif(200,-5,5)
#simuvarexpl = as.data.frame(simuvarexpl)
simuvarexpl=as.matrix(simuvarexpl) 
trueind <- c(10,20,30,40,50)
beta <- c(1,-1,2,-2,3)
ysimu <- simuvarexpl[,trueind]%*% beta + rnorm(200,0,2)

###############################################################
# Fonction de prédiction 
predictions <- function(matableTest,muchap,betachap)
{
ychap <- muchap * rep(1,dim(matableTest)[1])+ as.matrix(matableTest[,]) %*% betachap
return(ychap)
}

########################################################################
########################################################################
# Question 2
##############Simulation d'un coefficient aléatoire
# on a calculé la loi de beta | y (c'est une gaussienne)
# on voit dans la variance a posteriori que le bayésien fonctionne comme du ridge 
# la matrice devient inversible 
# on va la simuler (ici le Gibbs n'a qu'une seule variable !!!)
# Initialisaton des paramétres 
s2eps=4  # supposé connu ici 
s2beta= 1 # a priori ici éa influence le résultat 
# un s2beta trop petit va trop contraindre les beta 
# et ils auront du mal é s'éloigner de zéro
# pour aller vite on prend des tailles de burn in et d'itérations petites
niter=1000
nburn=200

###############################################################
# ci dessous la fonction qui simule le beta a posteriori 

BayesBasic  <- function(y,X,sigma2beta,sigma2eps,nbiter,nburn)
{
p <- dim(X)[2]
n <- dim(X)[1]
identite=diag(1,p,p)
### resultats a garder
resbeta <- matrix(rep(0,p*(nbiter-nburn)),nrow=p,
ncol=(nbiter-nburn))
### initialisation
beta <- rep(0,p)
### iterations
for (iter in 1 :nbiter)
{
print(iter)
Sigmabeta <- solve(t(X)%*%X/sigma2eps + identite/sigma2beta)
beta <- as.numeric(rmnorm(1,Sigmabeta%*%t(X)%*%y/sigma2eps,
Sigmabeta))
if (iter > nburn)
{
resbeta[,iter-nburn] <- beta
}
}
return(resbeta)
}
#### Illustrations 
# si on prend s2beta trop grand on a des beta trop variables
# et inversement, ils se rapprochent de zéro si sigmabeta est petit
res=BayesBasic(ysimu,simuvarexpl,1,4,500,200)
res=BayesBasic(ysimu,simuvarexpl,200,4,300,10)
# => trop dispersé 
res=BayesBasic(ysimu,simuvarexpl,0.1,4,300,10)
plot(sort(abs(rowMeans(res))))
plot(density(res[50,]))
plot(density(res[51,]))
# si on prend s2beta = 0.1 alors tous les coef vont étre trés resserrés 
# et les vrais seront sous estimés en valeur absolue 
# Inversement si on prend s2beta= 1000 alors les coef vont étre plus dispersés 
# et ils seront surestimés en valeur absolue  
niter=500
nburn=300
# attention on devrait prendre seulement les 100 premiéres observations (training)
resbayes=BayesBasic(ysimu[1:100],simuvarexpl[1:100,],s2beta,s2eps,niter,nburn)
plot(density(resbayes[1,])) # bruit 
plot(density(resbayes[5,])) # bruit 
mean(resbayes[5,])
plot(density(resbayes[40,])) # vrai coef
mean(resbayes[40,]) # vrai coef 
plot(density(resbayes[50,])) # vrai coef
mean(resbayes[50,]) # vrai coef
#rowMeans(resBayes)
#colMeans(t(resbayes))
plot(sort(abs(colMeans(t(resbayes))))) # on voit les coef qui se dégagent 
# on voit les 5 coef sortir du lot 

########################################################################
########################################################################
# Question 3
#### Algorithme avec variance de Zellner
# Pour extraire une  matrice inversible on peut prendre les premiéres variables
# on peut aussi se référer é l'indice de conditionnement
# Par exemple calculer les valeurs propres
X=as.matrix(simuvarexpl[1:100,]) 
valpX=eigen(t(X)%*%X, symmetric =TRUE, only.values = TRUE)
valpX
valmax=max(valpX$values)
valmin=min(valpX$values)
indcond=valmax/valmin
indcond
# on trouve ici quelque chose de négatif car ce sont les approximations 
# (normalement cela vaut zéro) 
# Prenons une sous matrice de taille 50*50 par exemple 
indalea=sample(c(1:300),50,replace=FALSE)
X2=X[,indalea] 


# c = coefficiet d'échelle 
# L = coefficient ridge
BayesZellner  <- function(y,X,c ,L, sigma2beta,sigma2eps,nbiter,nburn)
{
p <- dim(X)[2]
n <- dim(X)[1]
identite=diag(1,p,p)
### resultats à garder
resbeta <- matrix(rep(0,p*(nbiter-nburn)),nrow=p,
ncol=(nbiter-nburn))
### initialisation
beta <- rep(0,p)
### iterations
for (iter in 1 :nbiter)
{
print(iter)
Sigmabeta <- c*solve(t(X)%*%X/sigma2eps + L*identite/sigma2beta)
beta <- as.numeric(rmnorm(1,Sigmabeta%*%t(X)%*%y/sigma2eps,
Sigmabeta))
if (iter > nburn)
{
resbeta[,iter-nburn] <- beta
}
}
return(resbeta)
}
# Algorithme avec Zelnner (+ Ridge é => question supplémentaire) 
s2eps=4
s2beta=200
# si on prend s2beta = 0.1 alors tous les coef vont étre trés resserrés 
# et les vrais seront sous estimés en valeur absolue 
# Inversement si on prend s2beta= 1000 alors les coef vont étre plus dispersés 
# et ils seront surestimés en valeur absolue  
niter=300
nburn=100
c=1
L=100
# attention ici on prend seulement les 100 premiéres observations (training)
resbayes=BayesZellner(ysimu[1:100],as.matrix(simuvarexpl[1:100,]),c,L, s2beta,s2eps,niter,nburn)
######################################################################
plot(density(resbayes[1,])) # bruit 
plot(density(resbayes[5,])) # bruit 
mean(resbayes[5,])
plot(density(resbayes[40,])) # vrai coef
mean(resbayes[40,]) # vrai coef 
plot(density(resbayes[50,])) # vrai coef
mean(resbayes[50,]) # vrai coef
#colMeans(t(resbayes))
plot(sort(abs(colMeans(t(resbayes))))) # on voit les coef qui se dégagent 
# Rermarque : on pourrait ajouter un coefficient ridge pour stabiliser l'inverse.


#######################################################
#######################################################
# Question 4
# Algorithme avec  s2 aléatoire é 
# prenons sigma2eps aléatoire de loi InvGamma(c,d) 
# On montre que sa loi a posteriori est une InvGamma de paramétre 
# (c+n/2,d+(Y-mu-X beta)'(Y-mu-X beta)/2)
# on programme le Gibbs ainsi : 
library(invgamma)
BayesBasic2  <- function(y,X,sigma2beta,c,d,nbiter,nburn)
{
p <- dim(X)[2]
n <- dim(X)[1]
identite=diag(1,p,p)
### resultats a garder
resbeta <- matrix(rep(0,p*(nbiter-nburn)),nrow=p,ncol=(nbiter-nburn)) # paramétre vectoriel
reseps <- rep(0,(nbiter-nburn)) # paramétre univarié 
### initialisation
beta <- rep(0,p)
sigma2eps <- var(y)   # on intialise avec l'empirique
### iterations
for (iter in 1 :nbiter)
{
print(iter)
Sigmabeta <- solve(t(X)%*%X/as.numeric(sigma2eps) + identite/sigma2beta)
# premiére marge conditionnelle sachant s2eps
beta <- as.numeric(rmnorm(1,Sigmabeta%*%t(X)%*%y/as.numeric(sigma2eps),Sigmabeta))
s=t(y-X%*%beta)%*%(y-X%*%beta)
# deuxiéme loi conditionnelle sachant beta
sigma2eps <- rinvgamma(1,c+n/2,d+s/2)
if (iter > nburn)
{
resbeta[,iter-nburn] <- beta
reseps[iter-nburn] <- sigma2eps
}
}
resu=list(resbeta,reseps)
return(resu)
}
############ Illustrations 
# on fixe c=1=d pour l'inverse Gamma 
c=d=3  # => la variance du prior vaut l'infini !
# moyenne inverse gamma(c,d) = c/(d-1)
# variance inverse gamma(c,d) = c^2/((d-1)^2(d-2)) 
sigma2beta=10 # a priori sur la variance de beta
res2=BayesBasic2(ysimu,simuvarexpl,sigma2beta,c,d,500,200)
resbeta=res2[[1]]
reseps=res2[[2]] 
resbeta[,1]
plot(density(resbeta[1,]))
plot(density(resbeta[50,]))
plot(density(resbeta[40,]))
plot(density(reseps))
mean(resbeta[50,]) # moyenne a posteriori qui est l'estimateur bayésien
mean(reseps) # moyenne a posteriori 
plot(sort(abs(colMeans(t(resbeta)))))


########################################################################
########################################################################
# Question 5
# Ici on est trés proche du Bayes A qui sera étudié au TP1
# On passe... 

########################################################################


########################################################################
# Question 6
########################### # Régression bayésienne et ABC 
# Illustration avec un modéle de scale mixture
# Simulation du modéle
# on pose 
# Y_i = mu_i + X_i beta_i + e_i 
# avec e ~ N(0,s)
# mu ~uniforme (-a,a)
# beta_i ~ N(m,v) 
#  
#####################################################
# Les paramétres 
n=200 # nbre observations (simulations)
d=1 # taille beta 
m=6
v=1
s=1
m0=0
a=5
b=5
x <- matrix(rep(0,d*n),ncol=d,nrow=n)
for (i in 1:d){
x[,i] <- runif(n,-5,5)}
x <- as.data.frame(x)
X=as.matrix(x) 
beta <- rnorm(n,m,v) # chaque beta est aléatoire 
m0 <- runif(n,-a,a)
s <- rchisq(n,b)
y <- m0 + X%*%beta + rnorm(n,0,s)  # nos y simulés 



############################################################################################
# ABC basique 
############################################################################################
ABCbasic  = function(y, X, a, b, V, seuil, K)
{ 
d = dim(X)[2]
identite=diag(rep(1,d)) # matrice identité =1 (dimension =1)
identn=diag(rep(1,n))    # matrice identité n x n 
initbeta=rep(0,(d*K))
resbeta=matrix(initbeta,ncol=K)
resmu=rep(0,K)
ress2=rep(0,K)
n=length(y)
compteur=0
for (j in 1:(K-1))
{
mu=runif(1,-a,a) # a priori uniforme 
s2=rchisq(1,b) # a priori sur la variance du bruit 
s=sqrt(s2)
beta=rmnorm(1,0,V*identite ) # a priori sur beta 
m=mu+X%*%beta  # moyenne des y sachant les paramétres
res = rmnorm(1,m,s*identn) # simulation du y sachant les paramétres 
if (sqrt(sum((y-res)^2)/n) < seuil)  # acceptation rejet 
{
compteur = compteur+1  
resbeta[,compteur] =beta 
ress2[compteur] =s2 
resmu[compteur] = mu 
}
}
accept=sum(resmu!=0) # nbre d'acceptation 
betafinal=resbeta[,c(1:accept)]
resmufinal=resmu[c(1:accept)]
ress2final=ress2[c(1:accept)]
return(list(mu=resmufinal,beta=betafinal,s2=ress2final, taux=100*accept/(K-1)))
}
##################################################
a=5
b=5
V=40
seuil=10
K=5000
res = ABCbasic(y, X, a, b, V, seuil, K)
res$taux
mean(res$mu)
mean(res$s2)
mean(res$beta)
plot(density(res$beta))
# le déséquilibre peut venir du fait que l'on présente plus souvent des beta 
# é gauche qu'é droite de 5 
# il s'estompe si on diminue le seuil ou si on augmente V
# il y a un pb d'identifiabilité avec mu !!! 
# il faut enlever mu ou alors reporter la moyenne de beta sur mu !!! 

