##########################################################################################
### Retromed Spain-portugal Diversity - Landscape at the Metacommunity pondscape level ##############
##########################################################################################

setwd("~/Desktop/H2020_Retromed")

save.image("~/Desktop/H2020_Retromed/Retromed_Space.RData")
library(dplyr)
library(tibble)
library(ggplot2)
library(gridExtra)
###############################################################################################
# Vilavona
Vilavona<-cbind(VILAVONA,100)
Vilavona<-Vilavona[,-4]
colnames(Vilavona)[8]<-"J"
Vilavona<-Vilavona[,c(1:3, 5,7,6,4,8)]
M.distance.Vilavona<-as.matrix(dist(Vilavona[,2:3]))
Meta.comm.Vilavona<-Meta.comm.Albera
J.Vilavona<-Vilavona[,8]
#D50s<-seq(50,3000,,50)
D50s<-10^seq(log10(10),log10(3000),,50)
Gr<-Vilavona$grado
Cc<-Vilavona$cc
Bt<-Vilavona$bc
N<-nrow(M.distance.Vilavona)

############# 0 #############

out.d50.obs.landscape.Vilavona<-NULL
for (d in D50s){
  out.t<-iteration.model(replicas=112, Meta.pool=Meta.comm.Vilavona, m.pool=0.01, Js=J.Vilavona, id.module=Vilavona[,7],
                         M.dist=M.distance.Vilavona, D50=d, m.max=1,
                         id.fixed=NULL, D50.fixed=1000, m.max.fixed=1, comm.fixed=Meta.comm.Vilavona,
                         Lottery=F, it=2000, prop.dead.by.it=0.05, id.obs=1:ncol(M.distance.Vilavona))
  out.t<-resume.out(out.t)
  S<-out.t[[1]][11:(10+N)]
  Be.all<-out.t[[1]][(10+N+1):(10+N+N)]
  gamma<-out.t[[1]][10]
  print(c(length(d),length(S), length(Be.all), length(gamma),length(Bt), length(Gr), length(Cc)))
  out.t2<-cbind(d,S, Be.all, gamma,Bt, Gr, Cc)
  
  print(tail(out.t2))
  out.d50.obs.landscape.Vilavona<-rbind(out.d50.obs.landscape.Vilavona, out.t2)
}

a.vilavona<- out.d50.obs.landscape.Vilavona
save.image("~/Dropbox/H2020/RETROMED/Frontiers_Diciembre_2022/Frontiers_Respuestas/Retromed_H2020_febrero_2023.RData")
############# Figures Vilavona
#########################################3
### Plot results: alpha and beta diversity in species` dispersal gradient
### also considering the gradient in local ponds isolation
head(a)
#### ALPHA
ggplot(data.frame(a))+aes(d,S, col=Gr)+geom_point()+scale_x_log10()
ggplot(data.frame(a[which(a[,6]<quantile(a[,6],.05)),]))+aes(d,S, col=Gr)+geom_point()+scale_x_log10()
ii.min.g<-which(a[,7]<quantile(a[,7],.05)); 
ii.max.g<-which(a[,7]>quantile(a[,7],.95))
d<-a[ii.min.g,1]; S<-a[ii.min.g,2];  gg<-a[ii.min.g,6]; cc<-a[ii.min.g,7]; B<-a[ii.min.g,3]

ggplot(data.frame(d,S))+aes(d,S, col=gg*10)+geom_point()+theme_bw()
M.hill<-nls(S~S0+(Smax-S0)*(d^q)/(d50^q+d^q), start = list(S0=4, Smax=30, q=3, d50=200)   )
summary(M.hill)
ph<-coefficients(M.hill)
f.hill<-function(x) ph[1]+(ph[2]-ph[1])*(x^ph[3])/(ph[4]^ph[3]+x^ ph[3])

d2<-a[ii.max.g,1]; S2<-a[ii.max.g,2];  gg2<-a[ii.max.g,6]; cc2<-a[ii.max.g,7]; B2<-a[ii.max.g,3]
ggplot(data.frame(d2,S2))+aes(d,S, col=gg*10)+geom_point()+theme_bw()
M.hill.2<-nls(S2~S0+(Smax-S0)*(d2^q)/(d50^q+d2^q), start = list(S0=4, Smax=32, q=3, d50=100)   )
summary(M.hill.2)
ph2<-coefficients(M.hill.2)
f.hill.2<-function(x) ph2[1]+(ph2[2]-ph2[1])*(x^ph2[3])/(ph2[4]^ph2[3]+x^ ph2[3])
f.hill.delta<-function(x)(ph2[1]+(ph2[2]-ph2[1])*(x^ph2[3])/(ph2[4]^ph2[3]+x^ ph2[3]))/(ph[1]+(ph[2]-ph[1])*(x^ph[3])/(ph[4]^ph[3]+x^ ph[3]))
S.ratio<-ggplot(data.frame(a))+aes(d)+geom_blank()+geom_function(fun=f.hill.delta, col="red")+
  theme_classic()+labs(x="dispersal",y="Richness ratio")+
  scale_x_log10()+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=6,face="bold"))


S.o.Vilavona<-ggplot(data.frame(a))+aes(d,S, col=Gr)+geom_point()+
  labs(title="Vilavona (146 ponds)", y="Local Richness", x=NULL)+
  geom_function(fun=f.hill, col="white",  size=1.5)+ 
  geom_function(fun=f.hill, col="red", linetype="solid", size=.7)+
  geom_function(fun=f.hill.2, col="white",  size=1.5)+ 
  geom_function(fun=f.hill.2, col="salmon", linetype="solid", size=.7)+
  scale_x_log10()+theme_bw()
S.o.Vilavona<-S.o.Vilavona+annotation_custom(ggplotGrob(S.ratio), xmin = log10(400), xmax = log10(3000),  ymin = 1, ymax = 27)
#  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"))
S.o.Vilavona


#### BETA ###############################################


ii.min.g<-which(a[,6]<quantile(a[,6],.05)); ggplot(data.frame(a[ii.min.g,]))+aes(d,Be.all, col=Gr)+geom_point()+scale_x_log10()
ii.max.g<-which(a[,6]>quantile(a[,6],.95)); ggplot(data.frame(a[ii.max.g,]))+aes(d,Be.all, col=Gr)+geom_point()+scale_x_log10()

M.hill.b<-nls(B~S0+(Smax-S0)*(d^q)/(d50^q+d^q), start = list(S0=4.396e-01, Smax=9.179e-01, q=-2.56, d50=250) )
summary(M.hill.b)
ph.b<-coefficients(M.hill.b)
f.hill.b<-function(x) ph.b[1]+(ph.b[2]-ph.b[1])*(x^ph.b[3])/(ph.b[4]^ph.b[3]+x^ ph.b[3])

M.hill.2b<-nls(B2~S0+(Smax-S0)*(d2^q)/(d50^q+d2^q), start = list(S0=0, Smax=10, q=-3, d50=150)   )
summary(M.hill.2b)
ph2b<-coefficients(M.hill.2b)
f.hill.2b<-function(x) ph2b[1]+(ph2b[2]-ph2b[1])*(x^ph2b[3])/(ph2b[4]^ph2b[3]+x^ ph2b[3])
f.hill.deltab<-function(x)(ph2b[1]+(ph2b[2]-ph2b[1])*(x^ph2b[3])/(ph2b[4]^ph2b[3]+x^ ph2b[3]))/(ph.b[1]+(ph.b[2]-ph.b[1])*(x^ph.b[3])/(ph.b[4]^ph.b[3]+x^ ph.b[3]))
B.ratio<-ggplot(data.frame(a))+aes(d)+geom_blank()+geom_function(fun=f.hill.deltab, col="red")+
  theme_classic()+labs(x="dispersal",y="Beta ratio")+
  scale_x_log10()+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=6,face="bold"))
#B.ratio

B.o.Vilavona<-ggplot(data.frame(a))+aes(d,Be.all, col=Gr)+geom_point()+
  labs( y="Beta", x=NULL, title =" ")+
  scale_x_log10()+theme_bw()+
  geom_function(fun=f.hill.b, col="white",  size=1.5)+ 
  geom_function(fun=f.hill.b, col="red", linetype="solid", size=.7)+
  geom_function(fun=f.hill.2b, col="white",  size=1.5)+ 
  geom_function(fun=f.hill.2b, col="salmon", linetype="solid", size=.7)
B.o.Vilavona<-B.o.Vilavona + annotation_custom(ggplotGrob(B.ratio), xmin = log10(500), xmax = log10(3000),  ymin = .53, ymax = .99)


grid.arrange(S.o.Vilavona, B.o.Vilavona, ncol=2)# FIGURE 1

#################################################################################################################
# Vilavona NULL MODEL
# sequence of D50
random_vs_real_D50.Vilavona<-H2020_Random_Landscape_local_ponds.D50(n.random = 140,Meta.pool=Meta.comm.Vilavona, m.pool=0.01,  Js=Vilavona[,8],
                                                                 id.module=Vilavona[,7],
                                                                 D50.min=1, D50.max=2000, m.max=1,M.X.Y=Vilavona[,2:3],
                                                                 id.fixed=0, D50.fixed=100, m.max.fixed=1, comm.fixed=Meta.comm.Vilavona,
                                                                 Lottery=F, it=0, prop.dead.by.it=0, id.obs=1:nrow(Vilavona))

Vilavona.random<-plotea.null.landscape.D50(random_vs_real_D50.Vilavona)
random_vs_real_D50.Vilavona_Febrero_2023<-random_vs_real_D50.Vilavona

#################################################################################################################
proportion.removed<-seq(0.01,0.95,,30)
#############################################################################################################################
##############################
### Remove from lower to higher centrality for different dispersal levels
############################################################
ppp.random_d50.500.Vilavona<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Albera, m.pool=0.01,
                                            M.xy.J=Vilavona, free.cores=2, JJ_in=8,
                                            D50=500, m.max=1, loss.by="random",by_column_n=8,
                                            id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                            Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Gr.d50.500.Vilavona<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Albera, m.pool=0.01,
                                          M.xy.J=Vilavona, free.cores=2, JJ_in=8,
                                          D50=500, m.max=1, loss.by="by.column",by_column_n=4,
                                          id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                          Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)


#ppp.R.Bt.d50.500.Vilavona<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=28, Meta.pool=Meta.comm.Albera, m.pool=0.01,
#                                    M.xy.J=Vilavona, free.cores=2, JJ_in=8,
#                                    D50=500, m.max=1, loss.by="by.column",by_column_n=6,
#                                    id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
#                                    Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Cc.d50.500.Vilavona<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Albera, m.pool=0.01,
                                          M.xy.J=Vilavona, free.cores=2, JJ_in=8,
                                          D50=500, m.max=1, loss.by="by.column",by_column_n=5,
                                          id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                          Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

#ppp.R.Module.d50.500.Vilavona<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=28, Meta.pool=Meta.comm.Albera, m.pool=0.01,
#                                        M.xy.J=Vilavona, free.cores=2, JJ_in=8,
#                                        D50=500, m.max=1, loss.by="by.column",by_column_n=7,
#                                        id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
#                                        Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)


ppp.R.500.Vilavona<-list("random" = ppp.random_d50.500.Vilavona, "Degree"=ppp.R.Gr.d50.500.Vilavona, "Closennes" = ppp.R.Cc.d50.500.Vilavona)

#par(mfcol=c(2,4))
#plotea(a = ppp.R.500.Vilavona[[1]], b = ppp.R.500.Vilavona)
#ppp->ppp.random_d50.500
###################################################################################3



ppp.random_d50.1000.Vilavona<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Albera, m.pool=0.01,
                                             M.xy.J=Vilavona, free.cores=2, JJ_in=8,
                                             D50=1000, m.max=1, loss.by="random",by_column_n=8,
                                             id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                             Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Gr.d50.1000.Vilavona<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Albera, m.pool=0.01,
                                           M.xy.J=Vilavona, free.cores=2, JJ_in=8,
                                           D50=1000, m.max=1, loss.by="by.column",by_column_n=4,
                                           id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                           Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)


#ppp.R.Bt.d50.1000<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=28, Meta.pool=Meta.comm.Albera, m.pool=0.01,
#                                     M.xy.J=Vilavona, free.cores=2, JJ_in=8,
#                                     D50=1000, m.max=1, loss.by="by.column",by_column_n=6,
#                                     id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
#                                     Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Cc.d50.1000.Vilavona<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Albera, m.pool=0.01,
                                           M.xy.J=Vilavona, free.cores=2, JJ_in=8,
                                           D50=1000, m.max=1, loss.by="by.column",by_column_n=5,
                                           id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                           Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

#ppp.R.Module.d50.1000<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=28, Meta.pool=Meta.comm.Albera, m.pool=0.01,
#                                         M.xy.J=Vilavona, free.cores=2, JJ_in=8,
#                                         D50=1000, m.max=1, loss.by="by.column",by_column_n=7,
#                                         id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
#                                         Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)


ppp.R.1000.Vilavona<-list("random" = ppp.random_d50.1000.Vilavona, "Degree"=ppp.R.Gr.d50.1000.Vilavona, "Closennes" = ppp.R.Cc.d50.1000.Vilavona)
#plotea(a = ppp.R.1000.Vilavona[[1]], b = ppp.R.1000)
#ppp->ppp.random_d50.500



############

ppp.random_d50.100.Vilavona<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Albera, m.pool=0.01,
                                            M.xy.J=Vilavona, free.cores=2, JJ_in=8,
                                            D50=100, m.max=1, loss.by="random",by_column_n=8,
                                            id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                            Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Gr.d50.100.Vilavona<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Albera, m.pool=0.01,
                                          M.xy.J=Vilavona, free.cores=2, JJ_in=8,
                                          D50=100, m.max=1, loss.by="by.column",by_column_n=4,
                                          id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                          Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)


#ppp.R.Bt.d50.100.Vilavona<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=28, Meta.pool=Meta.comm.Albera, m.pool=0.01,
#                                    M.xy.J=Vilavona, free.cores=2, JJ_in=8,
#                                    D50=100, m.max=1, loss.by="by.column",by_column_n=6,
#                                    id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
#                                    Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Cc.d50.100.Vilavona<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Albera, m.pool=0.01,
                                          M.xy.J=Vilavona, free.cores=2, JJ_in=8,
                                          D50=100, m.max=1, loss.by="by.column",by_column_n=5,
                                          id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                          Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

#ppp.R.Module.d50.100<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=28, Meta.pool=Meta.comm.Albera, m.pool=0.01,
#                                        M.xy.J=Vilavona, free.cores=2, JJ_in=8,
#                                        D50=100, m.max=1, loss.by="by.column",by_column_n=7,
#                                        id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
#                                        Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)


ppp.R.100.Vilavona<-list("random" = ppp.random_d50.100.Vilavona, "Degree"=ppp.R.Gr.d50.100.Vilavona, "Closennes" = ppp.R.Cc.d50.100.Vilavona)
#plotea(a = ppp.R.100[[1]], b = ppp.R.100)
#ppp->ppp.random_d50.500
###################################################################################3



ppp.random_d50.2000.Vilavona<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Albera, m.pool=0.01,   # Random loss as refference
                                             M.xy.J=Vilavona, free.cores=2, JJ_in=8,
                                             D50=2000, m.max=1, loss.by="random",by_column_n=8,
                                             id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                             Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Gr.d50.2000.Vilavona<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Albera, m.pool=0.01,      # Degree
                                           M.xy.J=Vilavona, free.cores=2, JJ_in=8,
                                           D50=2000, m.max=1, loss.by="by.column",by_column_n=4,
                                           id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                           Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)


#ppp.R.Bt.d50.2000<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=28, Meta.pool=Meta.comm.Albera, m.pool=0.01,      # Betweenness
#                                     M.xy.J=Vilavona, free.cores=2, JJ_in=8,
#                                     D50=2000, m.max=1, loss.by="by.column",by_column_n=6,
#                                     id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
#                                     Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Cc.d50.2000.Vilavona<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Albera, m.pool=0.01,      # Closeness
                                           M.xy.J=Vilavona, free.cores=2, JJ_in=8,
                                           D50=2000, m.max=1, loss.by="by.column",by_column_n=5,
                                           id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                           Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

#ppp.R.Module.d50.2000<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=28, Meta.pool=Meta.comm.Albera, m.pool=0.01,   # Modules
#                                         M.xy.J=Vilavona, free.cores=2, JJ_in=8,
#                                         D50=2000, m.max=1, loss.by="by.column",by_column_n=7,
#                                         id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
#                                         Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)


ppp.R.2000.Vilavona<-list("random" = ppp.random_d50.2000.Vilavona, "Degree"=ppp.R.Gr.d50.2000.Vilavona, "Closennes" = ppp.R.Cc.d50.2000.Vilavona)


############################################################################################################
# D50 with the maximum effect on biodiversity
###################################################################################3
ii.max.diff<-which(f.hill.delta.Vilavona(seq(1,3000,,10000))==max(f.hill.delta.Vilavona(seq(1,3000,,10000))))
d50.max.S.ratio<-seq(1,3000,,10000)[ii.max.diff]
ii.max.diff.b<-which(f.hill.deltab.Vilavona(seq(1,3000,,10000))==min(f.hill.deltab.Vilavona(seq(1,3000,,10000))))
d50.max.B.ratio<-seq(1,3000,,10000)[ii.max.diff.b]

d50.max.diver<-(d50.max.S.ratio+d50.max.B.ratio)/2


ppp.random_d50.d50.max.diver.Vilavona<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Vilavona, m.pool=0.01,   # Random loss as refference
                                                      M.xy.J=Vilavona, free.cores=2, JJ_in=8,
                                                      D50=d50.max.diver, m.max=1, loss.by="random",by_column_n=8,
                                                      id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                                      Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Gr.d50.d50.max.diver.Vilavona<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Vilavona, m.pool=0.01,      # Degree
                                                    M.xy.J=Vilavona, free.cores=2, JJ_in=8,
                                                    D50=d50.max.diver, m.max=1, loss.by="by.column",by_column_n=4,
                                                    id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                                    Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Cc.d50.d50.max.diver.Vilavona<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Vilavona, m.pool=0.01,      # Closeness
                                                    M.xy.J=Vilavona, free.cores=2, JJ_in=8,
                                                    D50=d50.max.diver, m.max=1, loss.by="by.column",by_column_n=5,
                                                    id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                                    Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Cc.d50.d50.max.diver.Vilavona<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Vilavona, m.pool=0.01,      # Closeness
                                                    M.xy.J=Vilavona, free.cores=2, JJ_in=8,
                                                    D50=d50.max.diver, m.max=1, loss.by="by.column",by_column_n=5,
                                                    id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                                    Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.d50.max.diver.Vilavona<-list("random" = ppp.random_d50.d50.max.diver.Vilavona, "Degree"=ppp.R.Gr.d50.d50.max.diver.Vilavona, "Closennes" = ppp.R.Cc.d50.d50.max.diver.Vilavona)

###########################################################################################################3
#### Figure 3 plot previous results
par(mfcol=c(2,5), lwd=1.3, cex=0.7, mar=c(2.5,2.5,.5,.5))
plotea(a = ppp.R.d50.max.diver.Vilavona[[1]], b = ppp.R.d50.max.diver.Vilavona)
plotea(a = ppp.R.100.Vilavona[[1]], b = ppp.R.100.Vilavona)
plotea(a = ppp.R.500.Vilavona[[1]], b = ppp.R.500.Vilavona)
plotea(a = ppp.R.1000.Vilavona[[1]], b = ppp.R.1000.Vilavona)
plotea(a = ppp.R.2000.Vilavona[[1]], b = ppp.R.2000.Vilavona)



land.lost<-plotea.gg(a = ppp.R.d50.max.diver.Vilavona[[1]], b = ppp.R.d50.max.diver.Vilavona)
grid.arrange(S.o.Vilavona, B.o.Vilavona, Vilavona.random[[1]], Vilavona.random[[2]], land.lost[[1]],land.lost[[2]],ncol=6, widths=c(1.5,1.5,.6,.6, .6,.6))

plotea.gg(a = ppp.R.Cc.d50.d50.max.diver.Vilavona[[1]], b = ppp.R.Cc.d50.d50.max.diver.Vilavona, titulo = "Vilavona")
######################################################################################################

