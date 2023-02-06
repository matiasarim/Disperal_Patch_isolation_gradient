##########################################################################################
### Retromed Spain-portugal Diversity - Landscape at the Metacommunity pondscape level ##############
##########################################################################################

#setwd("~/Desktop/H2020_Retromed")

#save.image("~/Desktop/H2020_Retromed/Retromed_Space.RData")
library(dplyr)
library(tibble)
library(ggplot2)
library(gridExtra)
# Pnae
#xy.Pnae_centralidad <- read.delim2("~/Desktop/H2020_Pnae_Spain/xy.Pnae_centralidad.txt")

#xy.Pnae_centralidad<-xy.Pnae_centralidad[-which(xy.Pnae_centralidad$grado.Pnae>0.2),]
#xy.Pnae_centralidad<-cbind(xy.Pnae_centralidad,100)
#colnames(xy.Pnae_centralidad)[8]<-"J"
#M.distance.Pnae<-as.matrix(dist(xy.Pnae_centralidad[,2:3]))
#min(lower.tri( M.distance.Pnae))
#Meta.comm.Pnae<-sort(as.matrix(Vector_spp_abund[,4]),decreasing = T)/sum(sort(as.matrix(Vector_spp_abund[,4]),decreasing = T))
#plot(log(Meta.comm.Pnae))
#J.Pnae<-xy.Pnae_centralidad[,8]
##########################################################################################
# Pnae
Pnae<-cbind(PNAE,100)
Pnae<-Pnae[,-4]
colnames(Pnae)[8]<-"J"
Pnae<-Pnae[,c(1:3, 5,7,6,4,8)]
M.distance.Pnae<-as.matrix(dist(Pnae[,2:3]))
Meta.comm.Pnae<-Meta.comm
J.Pnae<-Pnae[,8]
#D50s<-seq(50,3000,,50)
D50s<-10^seq(log10(10),log10(3000),,50)
Gr<-Pnae$grado
Cc<-Pnae$cc
Bt<-Pnae$bc
N<-nrow(M.distance.Pnae)

############# 0 #############

out.d50.obs.landscape.Pnae<-NULL
for (d in D50s){
  out.t<-iteration.model(replicas=112, Meta.pool=Meta.comm.Pnae, m.pool=0.01, Js=J.Pnae, id.module=Pnae[,7],
                         M.dist=M.distance.Pnae, D50=d, m.max=1,
                         id.fixed=NULL, D50.fixed=1000, m.max.fixed=1, comm.fixed=Meta.comm.Pnae,
                         Lottery=F, it=2000, prop.dead.by.it=0.05, id.obs=1:ncol(M.distance.Pnae))
  out.t<-resume.out(out.t)
  S<-out.t[[1]][11:(10+N)]
  Be.all<-out.t[[1]][(10+N+1):(10+N+N)]
  gamma<-out.t[[1]][10]
  print(c(length(d),length(S), length(Be.all), length(gamma),length(Bt), length(Gr), length(Cc)))
  out.t2<-cbind(d,S, Be.all, gamma,Bt, Gr, Cc)
  
  print(tail(out.t2))
  out.d50.obs.landscape.Pnae<-rbind(out.d50.obs.landscape.Pnae, out.t2)
}

a.pnae<- out.d50.obs.landscape.Pnae
save.image("~/Dropbox/H2020/RETROMED/Frontiers_Diciembre_2022/Frontiers_Respuestas/Retromed_H2020_febrero_2023.RData")
############# Figures Pnae
#########################################3
### Plot results: alpha and beta diversity in species` dispersal gradient
### also considering the gradient in local ponds isolation
head(a)
#### ALPHA
ggplot(data.frame(a))+aes(d,S, col=Gr)+geom_point()+scale_x_log10()
ggplot(data.frame(a[which(a[,6]<quantile(a[,6],.05)),]))+aes(d,S, col=Gr)+geom_point()+scale_x_log10()
ii.min.g<-which(a[,7]<quantile(a[,7],.005)); 
ii.max.g<-which(a[,7]>quantile(a[,7],.995))
d<-a[ii.min.g,1]; S<-a[ii.min.g,2];  gg<-a[ii.min.g,6]; cc<-a[ii.min.g,7]; B<-a[ii.min.g,3]

ggplot(data.frame(d,S))+aes(d,S, col=gg*10)+geom_point()+theme_bw()
M.hill.pnae<-nls(S~S0+(Smax-S0)*(d^q)/(d50^q+d^q), start = list(S0=11.5, Smax=25, q=3, d50=500)   )
summary(M.hill.pnae)
ph.pnae<-coefficients(M.hill.pnae)
f.hill.pnae<-function(x) ph.pnae[1]+(ph.pnae[2]-ph.pnae[1])*(x^ph.pnae[3])/(ph.pnae[4]^ph.pnae[3]+x^ ph.pnae[3])

d2<-a[ii.max.g,1]; S2<-a[ii.max.g,2];  gg2<-a[ii.max.g,6]; cc2<-a[ii.max.g,7]; B2<-a[ii.max.g,3]
ggplot(data.frame(d2,S2))+aes(d2,S2, col=gg*10)+geom_point()+theme_bw()
M.hill.2.pnae<-nls(S2~S0+(Smax-S0)*(d2^q)/(d50^q+d2^q), start = list(S0=0, Smax=30, q=5, d50=100)   )
summary(M.hill.2.pnae)
ph2.pnae<-coefficients(M.hill.2.pnae)
f.hill.2.pnae<-function(x) ph2.pnae[1]+(ph2.pnae[2]-ph2.pnae[1])*(x^ph2.pnae[3])/(ph2.pnae[4]^ph2.pnae[3]+x^ ph2.pnae[3])
f.hill.delta.pnae<-function(x)(ph2.pnae[1]+(ph2.pnae[2]-ph2.pnae[1])*(x^ph2.pnae[3])/(ph2.pnae[4]^ph2.pnae[3]+x^ ph2.pnae[3]))/(ph.pnae[1]+(ph.pnae[2]-ph.pnae[1])*(x^ph.pnae[3])/(ph.pnae[4]^ph.pnae[3]+x^ ph.pnae[3]))
S.ratio<-ggplot(data.frame(a))+aes(d)+geom_blank()+geom_function(fun=f.hill.delta.pnae, col="red")+
  theme_classic()+labs(x="dispersal",y="Richness ratio")+
  scale_x_log10()+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=6,face="bold"))


S.o.Pnae<-ggplot(data.frame(a))+aes(d,S, col=Gr)+geom_point()+
  labs( y="Local Richness", x=NULL, title = "Pnae (255 ponds)")+
  geom_function(fun=f.hill.pnae, col="white",  size=1.5)+ 
  geom_function(fun=f.hill.pnae, col="red", linetype="solid", size=.7)+
  geom_function(fun=f.hill.2.pnae, col="white",  size=1.5)+ 
  geom_function(fun=f.hill.2.pnae, col="salmon", linetype="solid", size=.7)+
  scale_x_log10()+theme_bw()
S.o.Pnae<-S.o.Pnae+annotation_custom(ggplotGrob(S.ratio), xmin = log10(520), xmax = log10(2900),  ymin = 1, ymax = 22)
#  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"))
S.o.Pnae


#### BETA ###############################################


ii.min.g<-which(a[,6]<quantile(a[,6],.005)); ggplot(data.frame(a[ii.min.g,]))+aes(d,Be.all, col=Gr)+geom_point()+scale_x_log10()
ii.max.g<-which(a[,6]>quantile(a[,6],.995)); ggplot(data.frame(a[ii.max.g,]))+aes(d,Be.all, col=Gr)+geom_point()+scale_x_log10()

M.hill.b.pnae<-nls(B~S0+(Smax-S0)*(d^q)/(d50^q+d^q), start = list(S0=4.396e-01, Smax=9.179e-01, q=-2.56, d50=250) )
summary(M.hill.b.pnae)
ph.b.pnae<-coefficients(M.hill.b.pnae)
f.hill.b.pnae<-function(x) ph.b.pnae[1]+(ph.b.pnae[2]-ph.b.pnae[1])*(x^ph.b.pnae[3])/(ph.b.pnae[4]^ph.b.pnae[3]+x^ ph.b.pnae[3])

M.hill.2b.pnae<-nls(B2~S0+(Smax-S0)*(d2^q)/(d50^q+d2^q), start = list(S0=0, Smax=10, q=-1, d50=200)   )
summary(M.hill.2b.pnae)
ph2b.pnae<-coefficients(M.hill.2b.pnae)
f.hill.2b.pnae<-function(x) ph2b.pnae[1]+(ph2b.pnae[2]-ph2b.pnae[1])*(x^ph2b.pnae[3])/(ph2b.pnae[4]^ph2b.pnae[3]+x^ ph2b.pnae[3])
f.hill.deltab.pnae<-function(x)(ph2b.pnae[1]+(ph2b.pnae[2]-ph2b.pnae[1])*(x^ph2b.pnae[3])/(ph2b.pnae[4]^ph2b.pnae[3]+x^ ph2b.pnae[3]))/(ph.b.pnae[1]+(ph.b.pnae[2]-ph.b.pnae[1])*(x^ph.b.pnae[3])/(ph.b.pnae[4]^ph.b.pnae[3]+x^ ph.b.pnae[3]))
B.ratio.pnae<-ggplot(data.frame(a))+aes(d)+geom_blank()+geom_function(fun=f.hill.deltab.pnae, col="red")+
  theme_classic()+labs( x="dispersal",y="Beta ratio")+
  scale_x_log10()+
  theme(axis.text=element_text(size=5),axis.title=element_text(size=6,face="bold"))
#B.ratio

B.o.Pnae<-ggplot(data.frame(a))+aes(d,Be.all, col=Gr)+geom_point()+
  labs(y="Beta", x=NULL, title = "")+
  scale_x_log10()+theme_bw()+
  geom_function(fun=f.hill.b.pnae, col="white",  size=1.5)+ 
  geom_function(fun=f.hill.b.pnae, col="red", linetype="solid", size=.7)+
  geom_function(fun=f.hill.2b.pnae, col="white",  size=1.5)+ 
  geom_function(fun=f.hill.2b.pnae, col="salmon", linetype="solid", size=.7)
B.o.Pnae<-B.o.Pnae + annotation_custom(ggplotGrob(B.ratio.pnae), xmin = log10(510), xmax = log10(3100),  ymin = .7, ymax = .99)


grid.arrange(S.o.Pnae, B.o.Pnae, ncol=2)# FIGURE 1

#################################################################################################################
# Pnae NULL MODEL
# sequence of D50
random_vs_real_D50.Pnae<-H2020_Random_Landscape_local_ponds.D50(n.random = 140,Meta.pool=Meta.comm.Pnae, m.pool=0.01,  Js=Pnae[,8],
                                                           id.module=Pnae[,7],
                                                           D50.min=1, D50.max=2000, m.max=1,M.X.Y=Pnae[,2:3],
                                                           id.fixed=0, D50.fixed=100, m.max.fixed=1, comm.fixed=Meta.comm.Pnae,
                                                           Lottery=F, it=0, prop.dead.by.it=0, id.obs=1:nrow(Pnae))

Pnae.random<-plotea.null.landscape.D50(random_vs_real_D50.Pnae)
random_vs_real_D50.Pnae_Febrero_2023<-random_vs_real_D50.Pnae
#################################################################################################################
proportion.removed<-seq(0.01,0.95,,30)
#############################################################################################################################
##############################
### Remove from lower to higher centrality for different dispersal levels
############################################################
ppp.random_d50.500.Pnae<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Pnae, m.pool=0.01,
                                      M.xy.J=Pnae, free.cores=2, JJ_in=8,
                                      D50=500, m.max=1, loss.by="random",by_column_n=8,
                                      id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                      Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Gr.d50.500.Pnae<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Pnae, m.pool=0.01,
                                    M.xy.J=Pnae, free.cores=2, JJ_in=8,
                                    D50=500, m.max=1, loss.by="by.column",by_column_n=4,
                                    id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                    Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)


#ppp.R.Bt.d50.500.Pnae<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=28, Meta.pool=Meta.comm.Pnae, m.pool=0.01,
#                                    M.xy.J=Pnae, free.cores=2, JJ_in=8,
#                                    D50=500, m.max=1, loss.by="by.column",by_column_n=6,
#                                    id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
#                                    Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Cc.d50.500.Pnae<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Pnae, m.pool=0.01,
                                    M.xy.J=Pnae, free.cores=2, JJ_in=8,
                                    D50=500, m.max=1, loss.by="by.column",by_column_n=5,
                                    id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                    Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

#ppp.R.Module.d50.500.Pnae<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=28, Meta.pool=Meta.comm.Pnae, m.pool=0.01,
#                                        M.xy.J=Pnae, free.cores=2, JJ_in=8,
#                                        D50=500, m.max=1, loss.by="by.column",by_column_n=7,
#                                        id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
#                                        Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)


ppp.R.500.Pnae<-list("random" = ppp.random_d50.500.Pnae, "Degree"=ppp.R.Gr.d50.500.Pnae, "Closennes" = ppp.R.Cc.d50.500.Pnae)

#par(mfcol=c(2,4))
#plotea(a = ppp.R.500.Pnae[[1]], b = ppp.R.500.Pnae)
#ppp->ppp.random_d50.500
###################################################################################3



ppp.random_d50.1000.Pnae<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Pnae, m.pool=0.01,
                                       M.xy.J=Pnae, free.cores=2, JJ_in=8,
                                       D50=1000, m.max=1, loss.by="random",by_column_n=8,
                                       id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                       Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Gr.d50.1000.Pnae<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Pnae, m.pool=0.01,
                                     M.xy.J=Pnae, free.cores=2, JJ_in=8,
                                     D50=1000, m.max=1, loss.by="by.column",by_column_n=4,
                                     id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                     Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)


#ppp.R.Bt.d50.1000<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=28, Meta.pool=Meta.comm.Pnae, m.pool=0.01,
#                                     M.xy.J=Pnae, free.cores=2, JJ_in=8,
#                                     D50=1000, m.max=1, loss.by="by.column",by_column_n=6,
#                                     id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
#                                     Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Cc.d50.1000.Pnae<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Pnae, m.pool=0.01,
                                     M.xy.J=Pnae, free.cores=2, JJ_in=8,
                                     D50=1000, m.max=1, loss.by="by.column",by_column_n=5,
                                     id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                     Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

#ppp.R.Module.d50.1000<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=28, Meta.pool=Meta.comm.Pnae, m.pool=0.01,
#                                         M.xy.J=Pnae, free.cores=2, JJ_in=8,
#                                         D50=1000, m.max=1, loss.by="by.column",by_column_n=7,
#                                         id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
#                                         Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)


ppp.R.1000.Pnae<-list("random" = ppp.random_d50.1000.Pnae, "Degree"=ppp.R.Gr.d50.1000.Pnae, "Closennes" = ppp.R.Cc.d50.1000.Pnae)
#plotea(a = ppp.R.1000.Pnae[[1]], b = ppp.R.1000)
#ppp->ppp.random_d50.500



############

ppp.random_d50.100.Pnae<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Pnae, m.pool=0.01,
                                      M.xy.J=Pnae, free.cores=2, JJ_in=8,
                                      D50=100, m.max=1, loss.by="random",by_column_n=8,
                                      id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                      Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Gr.d50.100.Pnae<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Pnae, m.pool=0.01,
                                    M.xy.J=Pnae, free.cores=2, JJ_in=8,
                                    D50=100, m.max=1, loss.by="by.column",by_column_n=4,
                                    id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                    Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)


#ppp.R.Bt.d50.100.Pnae<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=28, Meta.pool=Meta.comm.Pnae, m.pool=0.01,
#                                    M.xy.J=Pnae, free.cores=2, JJ_in=8,
#                                    D50=100, m.max=1, loss.by="by.column",by_column_n=6,
#                                    id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
#                                    Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Cc.d50.100.Pnae<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Pnae, m.pool=0.01,
                                    M.xy.J=Pnae, free.cores=2, JJ_in=8,
                                    D50=100, m.max=1, loss.by="by.column",by_column_n=5,
                                    id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                    Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

#ppp.R.Module.d50.100<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=28, Meta.pool=Meta.comm.Pnae, m.pool=0.01,
#                                        M.xy.J=Pnae, free.cores=2, JJ_in=8,
#                                        D50=100, m.max=1, loss.by="by.column",by_column_n=7,
#                                        id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
#                                        Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)


ppp.R.100.Pnae<-list("random" = ppp.random_d50.100.Pnae, "Degree"=ppp.R.Gr.d50.100.Pnae, "Closennes" = ppp.R.Cc.d50.100.Pnae)
#plotea(a = ppp.R.100[[1]], b = ppp.R.100)
#ppp->ppp.random_d50.500
###################################################################################3



ppp.random_d50.2000.Pnae<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Pnae, m.pool=0.01,   # Random loss as refference
                                       M.xy.J=Pnae, free.cores=2, JJ_in=8,
                                       D50=2000, m.max=1, loss.by="random",by_column_n=8,
                                       id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                       Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Gr.d50.2000.Pnae<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Pnae, m.pool=0.01,      # Degree
                                     M.xy.J=Pnae, free.cores=2, JJ_in=8,
                                     D50=2000, m.max=1, loss.by="by.column",by_column_n=4,
                                     id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                     Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)


#ppp.R.Bt.d50.2000<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=28, Meta.pool=Meta.comm.Pnae, m.pool=0.01,      # Betweenness
#                                     M.xy.J=Pnae, free.cores=2, JJ_in=8,
#                                     D50=2000, m.max=1, loss.by="by.column",by_column_n=6,
#                                     id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
#                                     Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Cc.d50.2000.Pnae<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Pnae, m.pool=0.01,      # Closeness
                                     M.xy.J=Pnae, free.cores=2, JJ_in=8,
                                     D50=2000, m.max=1, loss.by="by.column",by_column_n=5,
                                     id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                     Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

#ppp.R.Module.d50.2000<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=28, Meta.pool=Meta.comm.Pnae, m.pool=0.01,   # Modules
#                                         M.xy.J=Pnae, free.cores=2, JJ_in=8,
#                                         D50=2000, m.max=1, loss.by="by.column",by_column_n=7,
#                                         id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
#                                         Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ii.max.diff<-which(f.hill.delta.pnae(seq(1,3000,,10000))==max(f.hill.delta.pnae(seq(1,3000,,10000))))
d50.max.S.ratio<-seq(1,3000,,10000)[ii.max.diff]
ii.max.diff.b<-which(f.hill.deltab.pnae(seq(1,3000,,10000))==min(f.hill.deltab.pnae(seq(1,3000,,10000))))
d50.max.B.ratio<-seq(1,3000,,10000)[ii.max.diff.b]

d50.max.diver<-(d50.max.S.ratio+d50.max.B.ratio)/2

ppp.R.Cc.d50.d50.max.diver.Pnae<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Pnae, m.pool=0.01,      # Closeness
                                          M.xy.J=Pnae, free.cores=2, JJ_in=8,
                                          D50=d50.max.diver, m.max=1, loss.by="by.column",by_column_n=5,
                                          id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                          Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.2000.Pnae<-list("random" = ppp.random_d50.2000.Pnae, "Degree"=ppp.R.Gr.d50.2000.Pnae, "Closennes" = ppp.R.Cc.d50.2000.Pnae)
#plotea(a = ppp.R.2000[[1]], b = ppp.R.2000)
#ppp->ppp.random_d50.500
#### PLot Lines
par(mfcol=c(2,5), lwd=1.3, cex=0.7, mar=c(2.5,2.5,.5,.5))

############################################################################################################
# D50 with the maximum effect on biodiversity
###################################################################################3
ii.max.diff<-which(f.hill.delta.pnae(seq(1,3000,,10000))==max(f.hill.delta.pnae(seq(1,3000,,10000))))
d50.max.S.ratio<-seq(1,3000,,10000)[ii.max.diff]
ii.max.diff.b<-which(f.hill.deltab.pnae(seq(1,3000,,10000))==min(f.hill.deltab.pnae(seq(1,3000,,10000))))
d50.max.B.ratio<-seq(1,3000,,10000)[ii.max.diff.b]

d50.max.diver<-(d50.max.S.ratio+d50.max.B.ratio)/2


ppp.random_d50.d50.max.diver.Pnae<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Pnae, m.pool=0.01,   # Random loss as refference
                                            M.xy.J=Pnae, free.cores=2, JJ_in=8,
                                            D50=d50.max.diver, m.max=1, loss.by="random",by_column_n=8,
                                            id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                            Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Gr.d50.d50.max.diver.Pnae<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Pnae, m.pool=0.01,      # Degree
                                          M.xy.J=Pnae, free.cores=2, JJ_in=8,
                                          D50=d50.max.diver, m.max=1, loss.by="by.column",by_column_n=4,
                                          id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                          Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Cc.d50.d50.max.diver.Pnae<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Pnae, m.pool=0.01,      # Closeness
                                          M.xy.J=Pnae, free.cores=2, JJ_in=8,
                                          D50=d50.max.diver, m.max=1, loss.by="by.column",by_column_n=5,
                                          id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                          Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.Cc.d50.d50.max.diver.Pnae<-H2020_Random.loss(proportion.removed=proportion.removed,n.random=112, Meta.pool=Meta.comm.Pnae, m.pool=0.01,      # Closeness
                                                   M.xy.J=Pnae, free.cores=2, JJ_in=8,
                                                   D50=d50.max.diver, m.max=1, loss.by="by.column",by_column_n=5,
                                                   id.fixed=0, D50.fixed=NULL, m.max.fixed=NULL, comm.fixed=NULL,
                                                   Lottery=FALSE, it=NULL, prop.dead.by.it=NULL)

ppp.R.d50.max.diver.Pnae<-list("random" = ppp.random_d50.d50.max.diver.Pnae, "Degree"=ppp.R.Gr.d50.d50.max.diver.Pnae, "Closennes" = ppp.R.Cc.d50.d50.max.diver.Pnae)

###########################################################################################################3
#### Figure 3 plot previous results
par(mfcol=c(2,5), lwd=1.3, cex=0.7, mar=c(2.5,2.5,.5,.5))
plotea(a = ppp.R.d50.max.diver.Pnae[[1]], b =ppp.R.d50.max.diver.Pnae)
plotea(a = ppp.R.100.Pnae[[1]], b = ppp.R.100.Pnae)
plotea(a = ppp.R.500.Pnae[[1]], b = ppp.R.500.Pnae)
plotea(a = ppp.R.1000.Pnae[[1]], b = ppp.R.1000.Pnae)
plotea(a = ppp.R.2000.Pnae[[1]], b = ppp.R.2000.Pnae)


land.lost<-plotea.gg(a = ppp.R.d50.max.diver.Pnae[[1]], b = ppp.R.d50.max.diver.Pnae)
grid.arrange(S.o.Pnae, B.o.Pnae, Pnae.random[[1]], Pnae.random[[2]], land.lost[[1]],land.lost[[2]],ncol=6, widths=c(1.5,1.5,.6,.6, .6,.6))
plotea.gg(a = ppp.R.Cc.d50.d50.max.diver.Pnae[[1]], b =ppp.R.Cc.d50.d50.max.diver.Pnae, titulo = "Pnae")
######################################################################################################

