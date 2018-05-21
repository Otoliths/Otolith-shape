

    ############################################################################################################################
    ##  Functions used in the tutorial "Investigating morphospace occupation in multi-scale ecological and evolutionary data  ##
	##                                  using regression trees: Case studies and perspectives"                                ##
    ##  Accompanies Vignon (2016) Evolutionary Biology                                                                        ##
	##  matthias.vignon@univ-pau.fr or mvignon@st-pee.inra.fr                                                                 ##
    ############################################################################################################################

###################################################
# Extract variable importance from any "mvpart" object
VarImp<-function(	b,			# a "mvpart" object.
					pretty=T	# Coloured plot when TRUE.
				){

	obj<-b
	splits<-which(obj$frame[,1]!="<leaf>")
	#nodes<-as.numeric(row.names(obj$frame)[splits]) 
	obj$frame[splits,]->f
	#ff<-f[rev(order(f$complexity)),]
	fff<-f[rev(order(obj$frame[splits,6])),c(1,6)]
	fff$var<-factor(fff$var)
	Lev<-nlevels(fff$var)
	Imp<-data.frame(Variables=levels(fff$var),Importance=rep(0,length=nlevels(fff$var)))
	for(i in 1:nlevels(fff$var)){
		Imp$Importance[i]<-sum(fff$complexity[which(fff$var==levels(fff$var)[i])])
		}
	Results<-Imp[rev(order(Imp$Importance)),]
	Absolute<-Results
	Results$Importance<-100/sum(Results$Importance)*Results$Importance

	if(pretty==T) {
		jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
		co<-jet.colors(round(max(Results[2]),0)*2)
		RRR<-matrix(0,nrow=200,ncol=Lev)#new matrix for barplot beside = F
		colnames(RRR)<-Results$Variables
		for (i in 1:Lev){
			j<-1
			while (j<=round(Results$Importance[i],0)*2){
				RRR[j,i]<-1
				j<-j+1
			}
		}
		RRR<-RRR/2
		nbsplit<-max(b$cptable[,2])
		barplot(RRR,col=co,border=co,names.arg=colnames(RRR),xlab="Variables selected in the CART model",ylab="Relative importance (%)")
		#title(main=paste("Relative importance for model",italic(substitute(b)),"using",max(b$cptable[,2]max(b$cptable[,2])),"splits"))
		title(main=substitute(paste("Relative importance for model ",italic(b)," using ",nbsplit," splits" )))
		par(new=TRUE);barplot(round(Results[,2],0),col="#6495ED00",xlab="",ylab="",axes=F)
	} else {
		#windows()
		barplot(Results[,2],names.arg=Results[,1],xlab="Variables selected in the CART model",ylab="Relative importance (%)",col="grey")
		title(main=paste("Relative importance for model",substitute(b),"using",max(b$cptable[,2]),"splits"))
	}
	#if(sum(b$frame$ncompete)>0) cat(paste("Please, note that",sum(b$frame$ncompete),"competing values, corresponding to",length(which(b$frame$ncompete!=0)),"splits were omitted"),"\n")
	#if(sum(b$frame$nsurrogate)>0) cat(paste("Please, note that",sum(b$frame$nsurrogate),"surrogate values, corresponding to",length(which(b$frame$nsurrogate !=0)),"splits were omitted"),"\n")
	#if(sum(b$frame$nsurrogate)<1&sum(b$frame$ncompete)<1) cat(paste("Note: No surrogate or competing value in this model","\n"))
	
	cat("Variance explained by model",substitute(b),"= ")
	cat(sum(Absolute$Importance)*100,"%",'\n')
	ans<-list(Absolute=Absolute,Relative=Results,R2=sum(Absolute$Importance))
}


###################################################
# Plot groups belonging as obtained using "mvpart" function. 2D
PP<-function(	model,		# a "mvpart" object.
				coord		# a dataframe containing 2D coordinates in columns(e.g as obtained by PCA, MDS etc.). X in the first column and Y in the second.
			)
{
    if(require("ade4")=="False"){
		print("trying to install mvpart")
		install.packages("ade4")
		if(require("ade4")){
		  print("ade4 installed and loaded")
		} else {
		  stop("could not install ade4. Please try manually")
		}
	}
  
	gr<-model$where
	aa<-1;gr2<-rep(1,length(gr))      # Renumber clusters sequentially.
	for(i in 1:length(levels(as.factor(gr)))) {
		gr2[which(gr==levels(as.factor(gr))[i])]<-i
	  }
	gr2<-as.factor(gr2)
	colnames(coord)[1]<-"X"
	colnames(coord)[2]<-"Y"
	plot(coord$X,coord$Y,xlab="PC1",ylab="PC2",col=gr2)
	grid()
	coul = rainbow(nlevels(gr2))
	for (i in 1:nlevels(gr2)){
		X<-coord$X[which(gr2==levels(gr2)[i])]
		Y<-coord$Y[which(gr2==levels(gr2)[i])]
		hpts <- chull(X,Y)
		hpts <- c(hpts, hpts[1])
		lines(X[hpts], Y[hpts], col = coul[i])
	}
	d<-data.frame(X=coord$X,Y=coord$Y)
	s.class(d[,1:2],fac=gr2, add.p = TRUE,col=coul,axesell=FALSE,cell=0)

}

###################################################
# Plot Voronoi tessellation for groups obtained using "mvpart" function. 2D
PV<-function(	model,			# a "mvpart" object.
				XY,				# a dataframe containing 2D coordinates in columns(e.g as obtained by PCA, MDS etc.). X in the first column and Y in the second.
				hull="TRUE"		# should convexhull be plotted (does not affect area measures that always considers hull).
			){
  if(require("deldir")=="FALSE"){
    print("trying to install deldir")
    install.packages("deldir")
    if(require("deldir")){
      print("deldir installed and loaded")
    } else {
      stop("could not install deldir. Please try manually")
    }
  }
    if(require("spatstat")=="FALSE"){
    print("trying to install spatstat")
    install.packages("spatstat")
    if(require("spatstat")){
      print("spatstat installed and loaded")
    } else {
      stop("could not install spatstat. Please try manually")
    }
  }
	if (dim(as.matrix(XY))[2]!=2){
    stop("XY must constain two colums (2D coordinates)")
    }
	x<-XY[,1]
	y<-XY[,2]
	gr<-model$where
	aa<-1;gr2<-rep(1,length(gr))      # Renumber clusters sequentially.
	for(i in 1:length(levels(as.factor(gr)))) {
		gr2[which(gr==levels(as.factor(gr))[i])]<-i
	  }
	G<-as.factor(gr2)# to be used as colors
	co = rainbow(nlevels(G))
	coul<-rep(0,length=length(gr2))
	for(i in 1:nlevels(G)){
		coul[which(G==levels(G)[i])]<-co[i]
	}
	
	
	if (hull=="TRUE"){
		o<-owin(xrange=range(2*x), yrange=range(2*y)) # range trop grand, permet de faire ultérieurement graphique avec les points du convexhull
		tess <- dirichlet(ppp(x, y,window=o))
		index<-as.numeric(tileindex(x, y, tess)) # il ne doit pas y avoir de valeurs vides
		o<-convexhull.xy(x,y)
		tess <- dirichlet(ppp(x, y,window=o))
		plot(ppp(x, y,window=o),col="white",main="Voronoi morphospace for groups")
		ar<-numeric(length=tess$n)
		t<-tile.areas(tess)
		for(i in 1:tess$n){
			polygon(tess$tiles[[index[i]]]$bdry[[1]],col=coul[i])
			points(x[i],y[i],col=coul[i])
			ar[i]<-t[index[i]]
		}
	} else {
		d <- deldir(x,y)
		td <- tile.list(d)
		plot(XY,axes=F,xlab="",ylab="",main="Voronoi morphospace for groups",asp=1)
		plot(td,fillcol=coul,close=F,showpoints=F,add=T)
		
		o<-owin(xrange=range(2*x), yrange=range(2*y))
		tess <- dirichlet(ppp(x, y,window=o))
		index<-as.numeric(tileindex(x, y, tess))
		o<-convexhull.xy(x,y)
		tess <- dirichlet(ppp(x, y,window=o))
		ar<-numeric(length=tess$n)
		t<-tile.areas(tess)
		for(i in 1:tess$n){
			ar[i]<-t[index[i]]
		}
		
	}
	
	Area<-numeric(length=nlevels(as.factor(G)))
	for(i in 1:nlevels(as.factor(G))){
		Area[i]<-sum(ar[which(G==levels(G)[i])])
	}
	Areaweight<-100*Area/sum(Area)
	id<-data.frame(Group=as.numeric(levels(as.factor(gr))),N=as.numeric(table(model$where)),Area=Area,Areaperc=Areaweight)
	
	ans<-list(Area=Area,Index=id)
}

###################################################
# Plot groups belonging as obtained using "mvpart" function. interactive 3D
P3D<-function(	model,	# a "mvpart" object.
				Coord	# a dataframe containing 3D coordinates in columns(e.g as obtained by PCA, MDS etc.). By default, X in the first column, Y in the second and Z in the third.
				)
{
    if(is(model)!="rpart"){
    stop("Model is not an object of class 'rpart'")}
    if(ncol(Coord)!=3){
    stop("Coord must constain 3D coordinates")}
    if(require("rgl")=="False"){
		print("trying to install rgl")
		install.packages("rgl")
		if(require("rgl")){
		  print("rgl installed and loaded")
		} else {
		  stop("could not install rgl. Please try manually")
		}
	}
    if(require("geometry")=="False"){
		print("trying to install geometry")
		install.packages("geometry")
		if(require("geometry")){
		  print("geometry installed and loaded")
		} else {
		  stop("could not install geometry. Please try manually")
		}
	}
model$where<-as.factor(model$where)
N<-nlevels(model$where)
coul = rainbow(N)
Vol<-numeric(N)
colnames(Coord)[1]<-"x"
colnames(Coord)[2]<-"y"
colnames(Coord)[3]<-"z"
for (i in 1:N){ # pour groupe N°i
  X<-Coord$x[which(model$where==levels(model$where)[i])]
  Y<-Coord$y[which(model$where==levels(model$where)[i])]
  Z<-Coord$z[which(model$where==levels(model$where)[i])]
  d<-cbind(X,Y,Z)
  Vol[i]<-convhulln(d,options="FA")$vol
  del<-delaunayn(d)
  s<-t(surf.tri(d,delaunayn(d)))
  rgl.points(d,col=coul[i])
  rgl.triangles(d[s,1], d[s,2], d[s,3],col=coul[i],alpha=0.6,clear=F)
  }
ans<-list(Volume=Vol)
}


###################################################
# % Variance explained by group differences (from model) for distance data
#
# Fit of the model is defined by relative error (RE, accessible in model$cptable); the total impurity of the leaves
# divided by the impurity of the root node (the undivided data). However, RE gives an over-optimistic estimate of how
# accurately a tree will predict for new data, and predictive accuracy is better estimated from the cross-validated
# relative error (CVRE [0-1]). Alternatively, one may be interested in the observed (REAL) variance explained by the model.
Perc<-function(	P,		# a 'dist' object.
				model	# a 'mvpart' object that contains group belonging (model$where).
			  ){
    if(require("vegan")=="False"){
		print("trying to install vegan")
		install.packages("vegan")
		if(require("vegan")){
		  print("vegan installed and loaded")
		} else {
		  stop("could not install vegan. Please try manually")
		}
	}
P<-as.dist(P)
Groups<-as.factor(model$where)
betadisper(P, Groups)->a
WithinSS<-rep(0,length=nlevels(Groups))
#WithinVar<-rep(0,length=nlevels(Groups))
for (i in 1:nlevels(Groups)){
WithinSS[i]<-sum(a$distances[which(Groups==levels(Groups)[i])])}

TotWinthinSS<-sum(WithinSS)
#TotWinthinVar<-TotWinthinSS/(dim(as.matrix(P))[1]-4)

G<-as.factor(rep(1,length=dim(as.matrix(P))[1]))
betadisper(P,G)->b
TotSS<-sum(b$distances)
#TotVar<-TotSS/(dim(as.matrix(P))[1]-1)
BetweenSS<-TotSS-TotWinthinSS
#BetweenVar<-BetweenSS/3
#pie(c(TotSS,BetweenSS))
cat("% variance explained by groups differences :")
100/TotSS*BetweenSS # R2 % variance explained
}

###################################################
# % Variance explained by asymmetry
# Special case of the above-mentioned 'Perc' function for pair objects
Asym<-function(Dist,Groups){
    if(require("vegan")=="False"){
		print("trying to install vegan")
		install.packages("vegan")
		if(require("vegan")){
		  print("vegan installed and loaded")
		} else {
		  stop("could not install vegan. Please try manually")
		}
	}
	if(class(Dist)!="dist"){
		Dist<-as.dist(Dist)
		#Dist<-as.dist(Dist, diag = FALSE, upper = T)
	}
	Groups<-as.factor(Groups)
	table(Groups)->a
	as.factor(a)->b
	if (any(levels(b)!=2)){
	stop("Not two shapes per individual, asymmetry impossible to quantify!")
	}
	betadisper(Dist, Groups)->a
	WithinSS<-rep(0,length=nlevels(Groups))
	#WithinVar<-rep(0,length=nlevels(Groups))
	for (i in 1:nlevels(Groups)){
	WithinSS[i]<-sum(a$distances[which(Groups==levels(Groups)[i])])}

	TotWinthinSS<-sum(WithinSS)
	#TotWinthinVar<-TotWinthinSS/(dim(as.matrix(P))[1]-4)

	G<-as.factor(rep(1,length=dim(as.matrix(Dist))[1]))
	betadisper(Dist,G)->b
	TotSS<-sum(b$distances)
	#TotVar<-TotSS/(dim(as.matrix(P))[1]-1)
	BetweenSS<-TotSS-TotWinthinSS
	#BetweenVar<-BetweenSS/3
	#pie(c(TotSS,BetweenSS))
	cat("% total variance explained by asymmetry: ",100-(100/TotSS*BetweenSS),'\n') # 1-R2 = % variance unexplained 
	asym<-1-(BetweenSS/TotSS)
	ans<-list(Side=asym)
}

###################################################
# plot convex hull
Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
} 


###################################################
# plot R2, apparent and relative(CV-based) Error for a mvpart object
rsqrpart<-function (x) 
    {
      if (!inherits(x, "rpart")) 
        stop("Not legitimate rpart")
      p.rpart <- x$cptable
      xstd <- p.rpart[, 5]
      xerror <- p.rpart[, 4]
      rel.error <- p.rpart[, 3]
      nsplit <- p.rpart[, 2]
      method <- x$method
     # if (!method == "anova") 
     #  cat("May not be applicable for this method\n")
      plot(nsplit, 1 - rel.error, xlab = "Number of Splits", ylab = "R-square",ylim = c(0, 1), type = "l")
      grid()
      A<-par("usr"); rect(A[1],A[3],A[2],A[4],col="#6495ED25")
      points(nsplit, 1 - rel.error,pch=21,col="black",bg="white")
      par(new = TRUE)
      plot(nsplit, 1 - xerror, xlab = "", ylab = "",ylim = c(0, 1), type = "l",lty=2)
      points(nsplit, 1 - xerror,pch=21,col="black",bg="grey")
      legend(2,0.2, c("Apparent", "X Relative"), lty = 1:2, col="black",pch=21,cex=1,box.col="black",pt.bg=c("white","grey"),bg="white",bty="0")
      ylim <- c(min(xerror - xstd) - 0.1, max(xerror + xstd) + 
        0.1)
      plot(nsplit, xerror, xlab = "Number of Splits", ylab = "X Relative Error", ylim = ylim, type = "l")
      grid()
      A<-par("usr"); rect(A[1],A[3],A[2],A[4],col="#6495ED25")
      points(nsplit, xerror,pch=21,col="black",bg="grey")
      segments(nsplit, xerror - xstd, nsplit, xerror + xstd)
      invisible()
}

###################################################
# Interaction plot. "res" must be obtained using the mvpart function from the mvpart package
MCC<-function(res,weight="TRUE",standard="TRUE",reord="TRUE",Fig=TRUE){
	if(missing(res)){
       stop("Tree is missing")
    }
	if(require("reshape2")=="FALSE"){
		print("trying to install reshape2")
		install.packages("reshape2")
		if(require("reshape2")){
		  print("reshape2 installed and loaded")
		} else {
		  stop("could not install reshape2. Please try manually")
		}
    }
    if(require("ggplot2")=="FALSE"){
		print("trying to install ggplot2")
		install.packages("ggplot2")
		if(require("ggplot2")){
		  print("ggplot2 installed and loaded")
		} else {
		  stop("could not install ggplot2. Please try manually")
		}
    } 
	if(weight!="TRUE"){ty<-"count"
		}else{ty<-"weight"}
	leaf<-as.numeric(rownames(res$frame[which(res$frame[,1]!="<leaf>"),]))
	r<-rownames(res$frame[which(res$frame[,1]!="<leaf>"),])
	s<-res$frame[which(res$frame[,1]!="<leaf>"),4]
	ss<-s[rev(order(s))]
	rr<-as.numeric(r[rev(order(s))])		# donne ordre de 'pas-leaf' par ordre décroissant de deviance)
	terms<-attr(res$terms,"term.labels")
	mat<-matrix(0,ncol=length(terms),nrow=length(terms))
	colnames(mat)<-rownames(mat)<-terms

	for (i in 2:length(leaf)){
		pat<-path(res,rr[i],print.it=F)
		for (j in 2:length(pat)){						# élément j de pat (chenmin d'un noeud). Avance facteur par facteur.
			for (k in 1:length(terms)){
				if(is.element(terms[k], pat[c(1:(j-1))])){	
					if(weight=="TRUE") {
						mat[which(rownames(mat)==pat[j]),k]<-mat[which(rownames(mat)==pat[j]),k]+1*ss[i] #  frequency si weigth
					} else{mat[which(rownames(mat)==pat[j]),k]<-mat[which(rownames(mat)==pat[j]),k]+1}
				}	
			}			# lecture de mat[] : Est-ce que le facteur de cette ligne est précédé par lélément de cette colonne
		}				#                    Est-ce que l'élément de cette colonne précède le facteur de cette ligne
	}
	#heatmap(mat,Rowv=NA,Colv=NA)
	if(reord!="TRUE"){
		reorder <- function(mat){
			hc <- hclust(as.dist(mat))
			mat <-mat[hc$order, hc$order]
		}
		mat <- reorder(mat)
	}
	diag(mat)<-NA
	if(standard=="TRUE"&weight=="TRUE") {mat<-mat/max(mat,na.rm = TRUE)}
	if(Fig==TRUE){
		if(weight!="TRUE"){
			res<-melt(mat)
			colnames(res)<-c("Variable1","Variable1","value")
			p<-ggplot(res, aes(res[,1], res[,2],group=Variable1),environment = environment()) +
			geom_tile(aes(fill = value)) + 
			geom_text(aes(fill = res$value, label = round(res$value, 1)),na.rm=T) +
			xlab("Primary predictor") +
			ylab("Interacting predictor") +
			# scale_fill_gradient(low = "darkblue", high = "red",name="count",na.value="black")
			# scale_fill_gradient(low = "darkblue", high = "yellow",name="count",na.value="black")
			# scale_fill_gradient(low = "darkgreen", high = "yellow",name="count",na.value="black")
			scale_fill_gradient(low = "white", high = "red",name="count",na.value="black")	
		} else {
			res<-melt(mat)
			#res$value<-res$value/max(na.omit(res$value))
			colnames(res)<-c("Variable1","Variable1","value")
			p<-ggplot(res, aes(res[,1], res[,2],group=Variable1),environment = environment()) +
			geom_tile(aes(fill = value)) + 
			geom_text(aes(fill = res$value, label = round(res$value, 1)),na.rm=T) +
			xlab("Primary predictor") +
			ylab("Interacting predictor") +
			# scale_fill_gradient(low = "darkblue", high = "red",name="count",na.value="black")
			# scale_fill_gradient(low = "darkblue", high = "yellow",name="count",na.value="black")
			# scale_fill_gradient(low = "darkgreen", high = "yellow",name="count",na.value="black")
			scale_fill_gradient(low = "white", high = "red",name="Relative weight",na.value="black")
		}
	plot(p)
	}
	# The x-axis is the sequence number of the primary predictor and the y-axis the sequence number of the potential interacting predictor. 
	# The intensity expresses the frequency/count when the potential interacting predictor precedes the primary predictor in a tree.
	ans<-list(interact=mat,type=ty,standard=standard)
}
	

plotBootMCC<-function(mat,standard="TRUE",weight="TRUE",reord="TRUE"){

	if(reord!="TRUE"){
		reorder <- function(mat){
			hc <- hclust(as.dist(mat))
			mat <-mat[hc$order, hc$order]
		}
		mat <- reorder(mat)
	}
	diag(mat)<-NA
	if(standard=="TRUE"&weight=="TRUE") {mat<-mat/max(mat,na.rm = TRUE)}
		if(weight!="TRUE"){
			res<-melt(mat)
			colnames(res)<-c("Variable1","Variable1","value")
			p<-ggplot(res, aes(res[,1], res[,2],group=Variable1),environment = environment()) +
			geom_tile(aes(fill = value)) + 
			geom_text(aes(fill = res$value, label = round(res$value, 1)),na.rm=T) +
			xlab("Primary predictor") +
			ylab("Interacting predictor") +
			scale_fill_gradient(low = "white", high = "red",name="count",na.value="black")	
		} else {
			res<-melt(mat)
			colnames(res)<-c("Variable1","Variable1","value")
			p<-ggplot(res, aes(res[,1], res[,2],group=Variable1),environment = environment()) +
			geom_tile(aes(fill = value)) + 
			geom_text(aes(fill = res$value, label = round(res$value, 1)),na.rm=T) +
			xlab("Primary predictor") +
			ylab("Interacting predictor") +
			scale_fill_gradient(low = "white", high = "red",name="Relative weight",na.value="black")
		}
	plot(p)
}

###########################################################################
## ALL THE FUNCTIONS BELOW ARE SUBFUNCTIONS/ROUTINES (not to be changed) ##
###########################################################################

# The "sedit", "substring.location", "substring2", "replace.substring.wild", "numeric.string" and "all.digits" functions (below) were extracted from the Hmisc package (Version 3.15-0), so that script will work without installing the complete package.

sedit <- function(text, from, to, test=NULL, wild.literal=FALSE)
{
  to <- rep(to, length=length(from))
  for(i in 1:length(text)) {
    s <- text[i]
    if(length(s))
      for(j in 1:length(from)) {
        old <- from[j]
        front <- back <- FALSE
        if(!wild.literal) {
          if(substring(old,1,1)=='^') {
            front <- TRUE;
            old <- substring(old,2)
          }

          if(substring(old,nchar(old))=='$') { 
            back <- TRUE; old <- substring(old, 1, nchar(old)-1)
          }
        }

        new <- to[j]

        lold <- nchar(old)
        if(lold > nchar(s))
          next

        ex.old <- substring(old, 1:lold, 1:lold)
        if(!wild.literal && any(ex.old=='*')) 
          s <- replace.substring.wild(s, old, new, test=test, front=front, back=back)
        else {
          l.s <- nchar(s)
          is <- 1:(l.s-lold+1)
          if(front)
            is <- 1

          ie <- is + lold - 1
          if(back)
            ie <- l.s

          ss <- substring(s, is, ie)
          k <- ss==old
          if(!any(k))
            next

          k <- is[k]
          substring2(s, k, k+lold-1) <- new
        }
      }

    text[i] <- s
  }

  text
}


substring.location <- function(text, string, restrict)
{
  if(length(text)>1)
    stop('only works with a single character string')
  
  l.text <- nchar(text)
  l.string <- nchar(string)
  if(l.string > l.text)
    return(list(first=0,last=0))
  
  if(l.string==l.text)
    return(if(text==string)
             list(first=1,last=l.text)
           else 
             list(first=0,last=0))

  is <- 1:(l.text-l.string+1)
  ss <- substring(text, is, is+l.string-1)
  k <- ss==string
  if(!any(k))
    return(list(first=0,last=0))
  
  k <- is[k]
  if(!missing(restrict))
    k <- k[k>=restrict[1] & k<=restrict[2]]
  
  if(length(k)==0)
    return(list(first=0,last=0))
  
  list(first=k, last=k+l.string-1)
}


## if(version$major < 5)  14Sep00
substring2 <- function(text, first, last=100000L)
  base::substring(text, first, last)

'substring2<-' <- function(text, first, last=100000, value)
{
  if(is.character(first)) {
    if(!missing(last))
      stop('wrong # arguments')
    
    return(sedit(text, first, value))  ## value was setto 25May01
  }

  lf <- length(first)

  if(length(text)==1 && lf > 1) {
    if(missing(last))
      last <- nchar(text)

    last <- rep(last, length=lf)
    for(i in 1:lf) {
      text <- paste(if(first[i]>1) 
                      substring(text, 1, first[i]-1),
                    value,
                    substring(text, last[i]+1), sep='')

      if(i < lf) {
        j <- (i+1):lf
        w <- nchar(value) - (last[i]-first[i]+1)
        first[j] <- first[j] + w  
        last[j] <- last[j] +  w
      }
    }

    return(text)
  }
  text <- paste(ifelse(first>1,substring(text, 1, first-1),''), value,
                substring(text, last+1), sep='')
  text
}


replace.substring.wild <- function(text, old, new, test=NULL, 
                                   front=FALSE, back=FALSE)
{
  if(length(text)>1)
    stop('only works with a single character string')

  if(missing(front) && missing(back)) {
    if(substring(old,1,1)=='^') {
      front <- TRUE;
      old <- substring(old,2)
    }

    if(substring(old, nchar(old))=='$') {
      back <- TRUE
      old <- substring(old, 1, nchar(old)-1)
    }
  }
  if((front || back) && old!='*') 
    stop('front and back (^ and $) only work when the rest of old is *')

  star.old <- substring.location(old,'*')
  if(length(star.old$first)>1)
    stop('does not handle > 1 * in old')
  
  if(sum(star.old$first)==0)
    stop('no * in old')
  
  star.new <- substring.location(new,'*')
  if(length(star.new$first)>1)
    stop('cannot have > 1 * in new')

  if(old=='*' && (front | back)) {
    if(front && back)
      stop('may not specify both front and back (or ^ and $) with old=*')
    
    if(length(test)==0)
      stop('must specify test= with old=^* or *$')
    
    et <- nchar(text)
    if(front) {
      st <- rep(1, et);
      en <- et:1
    } else {
      st <- 1:et;
      en <- rep(et,et)
    }

    qual <- test(substring(text, st, en))
    if(!any(qual))
      return(text)
    
    st <- (st[qual])[1]
    en <- (en[qual])[1]
    text.before <- if(st==1)''
                   else substring(text, 1, st-1)
    
    text.after  <- if(en==et)''
                   else substring(text, en+1, et)
    
    text.star   <- substring(text, st, en)
    new.before.star <-
      if(star.new$first>1) 
        substring(new, 1, star.new$first-1)
      else ''

    new.after.star <- if(star.new$last==length(new))''
                      else substring(new, star.new$last+1)

    return(paste(text.before, new.before.star, text.star, new.after.star,
                 text.after, sep=''))
  }

  old.before.star <- if(star.old$first==1)''
                     else substring(old, 1, star.old$first-1)
  
  old.after.star  <- if(star.old$last==nchar(old))''
                     else substring(old, star.old$first+1)

  if(old.before.star=='')
    loc.before <- list(first=0, last=0)
  else {
    loc.before <- substring.location(text, old.before.star)
    loc.before <- list(first=loc.before$first[1], last=loc.before$last[1])
  }

  if(sum(loc.before$first+loc.before$last)==0)
    return(text)

  loc.after <- if(old.after.star=='') list(first=0, last=0)
               else {
                 la <- substring.location(text, old.after.star, 
                                          restrict=c(loc.before$last+1,1e10))
                 lastpos <- length(la$first)
                 la <- list(first=la$first[lastpos], last=la$last[lastpos])
                 if(la$first+la$last==0)
                   return(text)

                 la
               }

  loc.star <- list(first=loc.before$last+1, 
                   last=if(loc.after$first==0) nchar(text)
                        else loc.after$first-1)
  
  star.text <- substring(text, loc.star$first, loc.star$last)
  if(length(test) && !test(star.text))
    return(text)

  if(star.new$first==0)
    return(paste(if(loc.before$first>1)substring(text,1,loc.before$first-1),
                 new, sep=''))

  new.before.star <- if(star.new$first==1)''
                     else substring(new, 1, star.new$first-1)
  new.after.star  <- if(star.new$last==nchar(new)) ''
                     else substring(new, star.new$first+1)

  paste(if(loc.before$first>1)substring(text,1,loc.before$first-1),
        new.before.star,
        substring(text,loc.star$first,loc.star$last),
        new.after.star,
        if(loc.after$last<nchar(text) && loc.after$last>0) 
          substring(text,loc.after$last+1),
        sep='')
}


## Some functions useful as test= arguments to replace.substring.wild, sedit
numeric.string <- function(string)
{
  ##.Options$warn <- -1  6Aug00
  oldopt <- options(warn=-1)
  on.exit(options(oldopt))
  !is.na(as.numeric(string))
}


all.digits <- function(string)
{
  k <- length(string)
  result <- logical(k)
  for(i in 1:k) {
    st <- string[i]
    ls <- nchar(st)
    ex <- substring(st, 1:ls, 1:ls)
    result[i] <- all(match(ex,c('0','1','2','3','4','5','6','7','8','9'),nomatch=0)>0)
  }
  
  result
}




# All the below functions belong to the "MRTsummaryplot.R" file and were extracted from the "MVPARTwrap_0.1-9.2.tar" tarball package.
# This allows using plot functions without installing the complete "MVPARTwrap" package that was archived due to its dependency to "mvpart".
# The codes below were maintained by Marie-Helene Ouellette with contributions from Pierre Legendre.
# In any case, "mvpart" needs to be installed to compute the analysis.


# ---------------------------------------------------------- #
# Get determinant species : this has to be done at each node #
# ---------------------------------------------------------- #

MRT<-function(obj,percent=10,species=NULL,LABELS=FALSE,...)
# This function recalculates the complexity of each node, giving the species' contribution to the R2 at each node, thus the output is the table 1 in Dea'th (2002). It contains the total species variance partitionned by species, by the  tree, and by the splits of the tree.

{
	
	#require(Hmisc)

# obj is the mvpart object
# percent : percentage considered discriminant species
# typeplot : a vector with states tree the usual tree output is needed, and triordi if the triplot ordination representation of the tree is of interest
# ... : further options for mvpart. See rpart.option, rpart and mvpart functions for more details
# species : A vector of species names used to build the table of partitioned variance. If it is NULL, the colnames of obj$y will be used.

	# Check if obj is of proper class
	if(class(obj)!='rpart') stop('obj is not of class rpart')
	
	# splits contains the line numbers of table 'frame' of the mvpart object that are nodes, not leafs
	splits<-which(obj$frame[,1]!="<leaf>")
	
	# node numbers NOT in order of their contribution (their is no need yet)
	nodes<-as.numeric(row.names(obj$frame)[splits])  

	# Initilizing the matrix to calculate the contributions to R2 of each species, at each node
	SSE<-mat.or.vec(length(nodes)*2+1,ncol(obj$y))
	MOYs<-mat.or.vec(length(nodes)*2+1,ncol(obj$y))
	moy<-apply(obj$y,MARGIN=2,FUN=mean)	
	MOY<-rep(moy,times=nrow(obj$y))
	MOY<-matrix(MOY,nrow(obj$y),length(moy),byrow=TRUE)
	SSE[1,]<-colSums((obj$y-MOY)^2)
	# List of the objects in each node
	LWHERE<-list()
	RWHERE<-list()
	
	compt_SSE<-2 # Already got the SST... so we start at the second slot
	
	# Loop that computes the SSEs for all nodes (2 per nodes, 1 per child), and create LWHERE and RWHERE
	
	for(i in 1:length(nodes)) 
	{

		# Initialisation of :
		
		# left child a leaf, than this will change to nodes[i]*2
		Lleaf=0
		# rigth child a leaf, than this will change to nodes[i]*2+1
		Rleaf=0
		# left child a node, than this will change to nodes[i]*2
		Lnode=0
		# rigth child a node , than this will change to nodes[i]*2+1
		Rnode=0
		
		# What are the two children ? Are they leafs or nodes ?
		# Left child ?
		ifelse(any(nodes==nodes[i]*2),Lnode<-nodes[i]*2,Lleaf<-nodes[i]*2) # is it a node or a leaf ?
		# Rigth child ?
		ifelse(any(nodes==nodes[i]*2+1),Rnode<-nodes[i]*2+1,Rleaf<-nodes[i]*2+1) # is it a node or a leaf ?
		
		
		# For each child seperatly (the bipartition) , 
		
		# -- # Left child # -- #
		
		if(Lleaf==0)
		{
			
			# Initialisation
			node=0
			leaf=0
			ifelse(any(nodes==Lnode*2),node<-Lnode*2,leaf<-Lnode*2)
			ifelse(any(nodes==Lnode*2+1),node<-c(node,Lnode*2+1),leaf<-c(leaf,Lnode*2+1))
			
			# We know at this point that 'node' is not empty and Lnode was a node (at least the right one)

			for(j in 1:length(nodes)) # Since left child is a node : check younger relatives for leafs
			{
				# Get this node's children if it is one decendant of the left node
				
				if(any(node==nodes[j])) # if the node we're workin on is in the node vector (left child of interest) : look at it's children and classify them
				{
					# Check if these children are nodes or leafs and update consequently 'node' and 'leaf' ('leaf' is really the important one, but we need 'node' for 'leaf')
					ifelse(any(nodes==nodes[j]*2),node<-c(node,nodes[j]*2),leaf<-c(leaf,nodes[j]*2))
					ifelse(any(nodes==nodes[j]*2+1),node<-c(node,nodes[j]*2+1),leaf<-c(leaf,nodes[j]*2+1))
				}
			}
			
			
	
		# Than create the Lwhere : which object are in the left node
		
		Lwhere<-0
		for(k in 1:length(leaf)) Lwhere<-c(Lwhere,which(obj$where==which(as.numeric(row.names(obj$frame))==leaf[k])))
		Lwhere<-Lwhere[-1]
		
		}
		
		if(Lnode==0)
		{
			
			Lwhere<-which(obj$where==which(as.numeric(row.names(obj$frame))==nodes[i]*2))
			
		}
		
		
		
		# -- # Rigth child # -- #
	
		if(Rleaf==0)
		{
			
			# Initialisation
			node=0
			leaf=0
			ifelse(any(nodes==Rnode*2),node<-Rnode*2,leaf<-Rnode*2)
			ifelse(any(nodes==Rnode*2+1),node<-c(node,Rnode*2+1),leaf<-c(leaf,Rnode*2+1))
			
			# We know at this point that 'node' is not empty and Rnode was a node (at least the right one)
			# Since right child is a node : check younger relatives for leafs
			for(j in 1:length(nodes)) 
			{
				# Get this node's children if it is one decendant of the rigth node
				
				if(any(node==nodes[j])) # if the node we're workin on is in the node vector (rigth child of interest) : look at it's children and classify them
				{
					# Check if these children are nodes or leafs and update consequently 'node' and 'leaf' ('leaf' is really the important one, but we need 'node' for 'leaf')
					ifelse(any(nodes==nodes[j]*2),node<-c(node,nodes[j]*2),leaf<-c(leaf,nodes[j]*2))
					ifelse(any(nodes==nodes[j]*2+1),node<-c(node,nodes[j]*2+1),leaf<-c(leaf,nodes[j]*2+1))
					
				}
			}
			
		
			# Than create the Rwhere : the objects in the rigth child
			Rwhere<-0
			for(k in 1:length(leaf)) Rwhere<-c(Rwhere,which(obj$where==which(as.numeric(row.names(obj$frame))==leaf[k])))
			Rwhere<-Rwhere[-1]
			
		}
	
		if(Rnode==0)
		{
			
			Rwhere<-which(obj$where==which(as.numeric(row.names(obj$frame))==nodes[i]*2+1))
			
		}
		
		
		
		
		# We now have all elements needed to calculate the SSEs and means
		if(length(Lwhere)==1)
		{
			moy<-obj$y[Lwhere,]
			MOY<-moy
			SSE[compt_SSE,]<-(moy-MOY)^2
			MOYs[compt_SSE,]<-moy
			compt_SSE<-compt_SSE+1
		}
		if(length(Lwhere)>1)
		{
			moy<-apply(obj$y[Lwhere,],MARGIN=2,FUN=mean)
			MOY<-rep(moy,times=nrow(obj$y[Lwhere,]))
			MOY<-matrix(MOY,nrow(obj$y[Lwhere,]),length(moy),byrow=TRUE)
			SSE[compt_SSE,]<-colSums((obj$y[Lwhere,]-MOY)^2)
			MOYs[compt_SSE,]<-moy
			compt_SSE<-compt_SSE+1
		}
		
		
		if(length(Rwhere)==1)
		{
			moy<-obj$y[Rwhere,]
			MOY<-moy
			SSE[compt_SSE,]<-(moy-MOY)^2
			MOYs[compt_SSE,]<-moy
			compt_SSE<-compt_SSE+1
		}
		if(length(Rwhere)>1)
		{
			moy<-apply(obj$y[Rwhere,],MARGIN=2,FUN=mean)
			MOY<-rep(moy,times=nrow(obj$y[Rwhere,]))
			MOY<-matrix(MOY,nrow(obj$y[Rwhere,]),length(moy),byrow=TRUE)
			SSE[compt_SSE,]<-colSums((obj$y[Rwhere,]-MOY)^2)
			MOYs[compt_SSE,]<-moy
			compt_SSE<-compt_SSE+1
		}
		
		# Keep the L and R where for futur use (i.e. objects present at each node/child)
		LWHERE[[i]]<-Lwhere
		RWHERE[[i]]<-Rwhere
		
		
	}	
	
	# Now that we have all SSEs... let's calculate the contribution to R2 of each species, for each node.
	
	# Build the node/leaf for each level matrix
	
	Tree<-mat.or.vec(length(nodes),length(as.numeric(row.names(obj$frame))))
	rownames(Tree)<-as.character(nodes)
	colnames(Tree)<-row.names(obj$frame)
	
	for(i in 1:length(nodes))
	{
		# For every node, which ones are nodes and leafs #
		
		# Repeat previous line (nodes r nodes, leafs r leaf)
		if(i!=1) Tree[i,]=Tree[i-1,]
		
		# Update : This node is a node
		Tree[i,which(as.numeric(row.names(obj$frame))==nodes[i])]<-1
		
		# Update : It's leafs are leafs
		Tree[i,which(as.numeric(row.names(obj$frame))==nodes[i]*2)]<-2
		Tree[i,which(as.numeric(row.names(obj$frame))==nodes[i]*2+1)]<-2
		
		# 
		
	}
	
	# Calculation of contribution #
	
	SST<-sum(SSE[1,])
	
	R2<-mat.or.vec(length(nodes),ncol(obj$y))
	
	for(i in 1:length(nodes))
	{
		
		## leafs before ##
		if(i!=1)
		{
			bef_pos<-which(Tree[i-1,]==2)
			leafs<-1
			for(j in 1:length(bef_pos))
			{
				
				leaf<-as.numeric(colnames(Tree)[bef_pos[j]])
			
				if(leaf%%2==0)
				{
					# if it is a left node
					# where is it in SSE ?
					leafs<-cbind(leafs,which(nodes==leaf/2)*2)
				}
			
				if(leaf%%2==1)
				{
					# if it is a right node
					# where is it in SSE ?
					leafs<-cbind(leafs,which(nodes==(leaf-1)/2)*2+1)
				}
			}
			
			if(length(leafs)!=2)
			{
				
				bef<-apply(SSE[leafs[-1],],FUN=sum,MARGIN=2) # sum of leaf
			}
			if(length(leafs)==2)
			{
				
				bef<-SSE[leafs[-1],] # sum of leaf
			}
		}
		if(i==1) bef<-SSE[1,]
		
		
		## leafs now ##
		now_pos<-which(Tree[i,]==2)
		leafs<-1
		for(j in 1:length(now_pos))
		{
			leaf<-as.numeric(colnames(Tree)[now_pos[j]])
			
			if(leaf%%2==0)
			{
				# if it is a left node
				# where is it in SSE ?
				leafs<-cbind(leafs,which(nodes==leaf/2)*2)
			}
			
			if(leaf%%2==1)
			{
				# if it is a left node
				# where is it in SSE ?
				leafs<-cbind(leafs,which(nodes==(leaf-1)/2)*2+1)
			}
		}
		
		if(length(leafs)!=2)
		{
			
			now<-apply(SSE[leafs[-1],],FUN=sum,MARGIN=2) # sum of leaf
		}
		if(length(leafs)==2)
		{
			
			now<-SSE[leafs[-1],] # sum of leaf
		}
		
		## R2 for this node ##
		R2[i,]<-(bef-now)/SST
	}
	
	# Now get the discriminant species #
	
	
	sumR2<-apply(R2, 1,sum)

	mat_sum<-rep(sumR2, ncol(R2))
	mat_sum<-matrix(mat_sum,nrow(R2),ncol(R2))
	pourct<-(R2/mat_sum)*100
	colnames(pourct)<-colnames(obj$y)
	
	
	
	# Give nodes in real order of R2
	
	#R2_perct<-R2*100
	#R2_pernode<-apply(R2_perct,MARGIN=1,FUN=sum)
	#ii<-order(R2_pernode)
	#nodes<-as.numeric(row.names(obj$frame)[splits][ii])
	

	# Return results needed for complete summary
	
	
	MOYs<-MOYs[-1,]
	
	# Build table as table 1 in Dea'th (2002) article
	
	# Species names
	
	ifelse(length(species)==0,col1names<-colnames(obj$y),col1names<-species)

	coli<-t(R2)
	
	col_total_tree<-rowSums(coli)
	col_total_species<-SSE[1,]/sum(SSE[1,])
	
	#TABLE1
	TABLE1<-cbind(coli,col_total_tree,col_total_species)
	
	rownames(TABLE1)<-col1names
	
	# Add sum of columns
	total_col<-colSums(TABLE1)
	names(total_col)<-'SUMS'
	
	TABLE1<-rbind(TABLE1,total_col)
	
	ii<-ORDER_NODES(obj,R2,splits)
	if(LABELS==FALSE) mat_labels_NOTINORDER<-mat_labels_fct(obj)
	if(!LABELS==FALSE) mat_labels_NOTINORDER<-LABELS
	mat_labels<-mat_labels_NOTINORDER[ii,]
	if(length(splits)>1) colnames(TABLE1)[1:(ncol(TABLE1)-2)]=sedit(mat_labels[,2],c('<','>','=',' '),c('','','',''))
	
	if(length(splits)==1) colnames(TABLE1)[1:(ncol(TABLE1)-2)]=sedit(mat_labels[2],c('<','>','=',' '),c('','','',''))
	
	TABLE1<-TABLE1*100
	
	# nŽcessaire pour plot :
	# typeplot,Cex,widthtree,heighttree,widthtriordi,heighttriordi,test.F=test.F,silent=silent
	
	res<-list(nodes,pourct,R2,obj,percent,MOYs,RWHERE,LWHERE,TABLE1,LABELS,mat_labels)
	
	
	names(res)<-c('nodes','pourct','R2','obj','percent','MOYs','RWHERE','LWHERE','TABLE1','LABELS','mat_labels')
	class(res)<-'MRT'
	return(res)
	

	
}

### ------------------------------------ ###
###              plotmvpart              ###
### Gives usual tree, TRIplot and BIplot ###
### ------------------------------------ ###

plot.MRT<-function(x,NodeMarking=TRUE,typeplot=c('tree'),Cex=0.5,widthtree=7, heighttree=9,R2A=FALSE,X=NULL,T=NULL,...)
#ADD test.F=FALSE,silent=TRUE,'triordi',widthtriordi=7, heighttriordi=9,
{
	if(R2A & is.null(X)) stop ('If R2A is TRUE, X must be provided (non NULL)')
	if(R2A & is.null(T)) stop ('If R2A is TRUE, T must be provided (non NULL)')
	# If we want the tree plot
	if(any(typeplot=='tree'))
	{
		plot_tree(x$obj,NodeMarking,x$R2,Cex,widthtree,heighttree,x$LABELS,R2A=R2A,X=X,T=T,...)
	}

	# If we want the ordination plot
#	if(any(typeplot=='triordi'))
#	{
#		
#		plot_triplot(x$obj,x$LWHERE,x$RWHERE,x$R2,Cex,widthtriordi,heighttriordi,test.F=test.F,silent=silent,x$LABELS)
#	}
	
}

plot_tree<-function(obj,NodeMarking,R2,Cex,widthtree,heighttree,LABELS,R2A=FALSE,X=NULL,T=NULL,...)
{

	#require(Hmisc)

	# points in node
	a=3
	# points between nodes
	b=2
	
	# nodes designated by numbers following the frame imposed by mvpart
	nodes<-as.numeric(rownames(obj$frame))
	
	maxs=max(nodes)
	
	# RElates to how long the x axis should be (don't remember exaclty what this is)
	# Beleive this is how many maximum splits there is from the root node (leaf counting as a node)
	s<-floor(log(maxs)/log(2))+1
	# Splits : which one are not leafs, in others words their position in the table obj$frame #
	splits<-which(obj$frame[,1]!="<leaf>")

	# length of x axis is (a+b)s+b   3 |2|  3  |2|  3
	#par(mar=c(5,5,4,2),cex.axis=1,col.axis=184)
	
	# Getting the real order of the nodes (correction for children explaining more variance than parents)
	ii<-ORDER_NODES(obj,R2,splits)
	if(LABELS==FALSE) mat_labels_NOTINORDER<-mat_labels_fct(obj)
	if(!LABELS==FALSE) mat_labels_NOTINORDER<-LABELS
	R2_node<-rowSums(R2)
	
	
	mat_labels<-mat_labels_NOTINORDER[ii,]
	
	
	# Cumulative complexity (cumulative R square of the tree from top to bottom)
	#cum_R<-cumsum(R2_node[ii[order(length(ii):1)]])[order(length(ii):1)][order(ii)]
	
	#cum_R<-cumsum(R2_node[order(-ii)])[order(-c(1:length(R2_node)))][ii]
	
	#cum_R<-cumsum(R2_node[ii])
	
	cum_R<-c(0,cumsum(R2_node[ii]))
	# If we want to show to up to the maximum value of R2
	#quartz(width=widthtree, height=heighttree)
	par(cex=Cex)#,mar=c(5,5,0,2))
	plot(c(0,(a+b)*2^(s-2)+b),c(max(cum_R*100),0),xlim=c(0,(a+b)*2^(s-2)+b),ylim=c(max(cum_R*100),0),xlab="",ylab=expression(R^2),type='n',xaxt='n', yaxt='n',bty='n',plt=c(0.9,0.9,0.9,0.9))

		
	
	# If we want to show the whole axis no matter what ??? not shure about this
	#else
	#plot(c(0,(a+b)*2^(s-2)+b),c(100,0),xlab="",ylab=expression(R^2),type='n',xaxt='n', yaxt='n')
	
	#axis(side=2, at = c(0,c(0,cum_R*100)), labels = as.character(round(c(0,cumsum(R2_node[ii[order(length(ii):1)]][c(length(R2_node):1)]),1)*100,digits=2))[(length(R2_node)+2):1],pos=-0.75,las=1)
	
	axis(side=2, at = cum_R*100, labels = as.character(round(cum_R*100,digits=2)),pos=-0.5,las=1)
	
	# Buils mat_labels and mat_done for pretty labels attibution #
	# See if thisis done before #
	# call matlabels #
	
	#print(length(R2_node))
	
	
	# Get rid of data.matrix once and for all
	if(length(R2_node)==1) mat_labels<-as.data.frame(t(mat_labels))
	mat_labels[,2]<-substring(mat_labels[,2],substring.location(mat_labels[1,2], ').')$last+1)
	mat_labels[,3]<-substring(mat_labels[,3],substring.location(mat_labels[1,3], ').')$last+1)
	
	# Build vertical lines at the R2 level , and build horiontal lines #
	
	compt<-1
	complexity<-obj$frame$complexity
	
	# What if we change the complexity vector right away ?
	
	for(i in 1:length(nodes))
	{
		if(obj$frame$var[i]!='<leaf>') 
		{
			complexity[i]<-R2_node[compt]
			compt<-compt+1
		}	
	}
	
	
	compt<-1
	for(i in 1:length(nodes))
	{	
		# If the node is not a leaf
		if(obj$frame$var[i]!='<leaf>')
		{
			
			# This node is at which level of the hierarchy ?
			si<-floor(log(nodes[i])/log(2))+1
			# At this level, which position from the left does it occupy ?
			posi<-nodes[i]-2^(si-1)+1
			# Change complexity for the right one
			#complexity[i]<-R2_node[compt]
			
			
			# Horizontale lines, R2 printed, size of node if leaf
			
			if(si==(s-1))
			{
				min<-(posi-1)*(a+b)+b
				max<-posi*(a+b)
				lines(c(min,max),c(cum_R[which(ii==compt)]*100,cum_R[ which(ii==compt)]*100))
				#points(-1,obj$frame$complexity[i],pch=23)
				
				# R under the node
				middle<-mean(c(min,max))
				
				aaa<-as.numeric(R2_node[compt])*100
				aaa<-format(aaa,digits=2,nsmall=2)
				
				pretty_middle<-paste('(',nodes[i],')\n',aaa,'%',sep='')
				text(middle,cum_R[which(ii==compt)]*100,pos=1,labels=pretty_middle,cex=1)
					
				# Here we add the labels on the node ?
				caract_min<-mat_labels[which(as.numeric(mat_labels[,1])==nodes[i]),2] #which(as.numeric(splits[,1])==nodes[i])
				caract_max<-mat_labels[which(as.numeric(mat_labels[,1])==nodes[i]),3]
				text(min,cum_R[which(ii==compt)]*100-0.2,caract_min,pos=2,cex=1,offset=0.5)
				text(max,cum_R[which(ii==compt)]*100-0.2,caract_max,pos=4,cex=1,offset=0.5)
				compt=compt+1
			}

			if(si==(s-2))
			{
				min<-(2*posi-1)*(a+b)-(1/2)*a
				max<-2*posi*(a+b)-(1/2)*a
				lines(c(min,max),c(cum_R[which(ii==compt)]*100,cum_R[which(ii==compt)]*100))
				# points(-1,obj$frame$complexity[i],pch=23)
				
				middle<-mean(c(min,max))
				
				aaa<-as.numeric(R2_node[compt]*100)
				aaa<-format(aaa,digits=2,nsmall=2)
				
				pretty_middle<-paste('(',nodes[i],')\n',aaa,'%',sep='')
				text(middle,cum_R[which(ii==compt)]*100,pos=1,labels=pretty_middle,cex=1)
				#text(-1,obj$frame$complexity[i],pos=4,labels=format(obj$frame$complexity[i],digits=3),cex=0.5)
				caract_min<-mat_labels[which(as.numeric(mat_labels[,1])==nodes[i]),2]
				caract_max<-mat_labels[which(as.numeric(mat_labels[,1])==nodes[i]),3]
				text(min,cum_R[which(ii==compt)]*100-0.2,caract_min,pos=2,cex=1,offset=0.5)
				text(max,cum_R[which(ii==compt)]*100-0.2,caract_max,pos=4,cex=1,offset=0.5)
				compt=compt+1
			}

			if(si<(s-2))
			{
				# DiffÃ©rence de 2^(s-si-3)
				min<-(2^(s-si-2)+(posi-1)*2^(s-si-1)-2^(s-si-3))*(a+b)+(1/2)*b
				max<-(2^(s-si-2)+(posi-1)*2^(s-si-1)+2^(s-si-3))*(a+b)+(1/2)*b
				# Vertical line
				lines(c(min,max),c(cum_R[which(ii==compt)]*100,cum_R[which(ii==compt)]*100))
				# points(-1,obj$frame$complexity[i],pch=23)
				
				middle<-mean(c(min,max))
				
				aaa<-as.numeric(R2_node[compt])*100
				aaa<-format(aaa,digits=2,nsmall=2)
				
				pretty_middle<-paste('(',nodes[i],')\n',aaa,'%',sep='')
				text(middle,cum_R[which(ii==compt)]*100,pos=1,labels=pretty_middle,cex=1)
				#text(-1,obj$frame$complexity[i],pos=4,labels=format(obj$frame$complexity[i],digits=3),cex=0.5)
				caract_min<-mat_labels[which(as.numeric(mat_labels[,1])==nodes[i]),2]
				caract_max<-mat_labels[which(as.numeric(mat_labels[,1])==nodes[i]),3]
				text(min,cum_R[which(ii==compt)]*100-0.2,caract_min,pos=2,cex=1,offset=0.5)
				text(max,cum_R[which(ii==compt)]*100-0.2,caract_max,pos=4,cex=1,offset=0.5)
				compt=compt+1
			}
			
			
			# And the verticals #
			
			# Find the level of R2 of the children, if it's not a leaf
			
			# Number of the child to the right and the child to the left
			child_right<-which(as.numeric(row.names(obj$frame))==as.numeric(row.names(obj$frame)[i])*2)
			child_left<-which(as.numeric(row.names(obj$frame))==as.numeric(row.names(obj$frame)[i])*2+1)
			
			R_child_right<-complexity[child_right]
			R_child_left<-complexity[child_left]
			
			# These are the complexity of the children, if they are not leafs.
			posR_child_right<-match(R_child_right,R2_node[ii],nomatch=999)
			posR_child_left<-match(R_child_left,R2_node[ii],nomatch=999)
			
		
			
			if(posR_child_right==999) R_cum_child_right=max(cum_R)
			if(posR_child_right!=999) R_cum_child_right<-cum_R[posR_child_right]
			if(posR_child_left==999) R_cum_child_left=max(cum_R)
			if(posR_child_left!=999) R_cum_child_left<-cum_R[posR_child_left]
			
			
			if(NodeMarking)
			{
				if(obj$frame[child_right,1]=='<leaf>') text(min,max(cum_R)*100,paste('n=',obj$frame[child_right,2],'\n(#',child_right,')',sep=''),cex=0.75,pos=1)
			
				if(obj$frame[child_left,1]=='<leaf>') text(max,max(cum_R)*100,paste('n=',obj$frame[child_left,2],'\n(#',child_left,')',sep=''),cex=0.75,pos=1)
			}
			
			if(!NodeMarking)
			{
				if(obj$frame[child_right,1]=='<leaf>') text(min,max(cum_R)*100,paste('n=',obj$frame[child_right,2],sep=''),cex=0.75,pos=1)
			
				if(obj$frame[child_left,1]=='<leaf>') text(max,max(cum_R)*100,paste('n=',obj$frame[child_left,2],sep=''),cex=0.75,pos=1)
			}
			
			
			
			# Draw both vertical lines
			# For right child
		
			lines(c(min,min),c(cum_R[which(ii==(compt-1))]*100,R_cum_child_right*100))
			# For left child
			lines(c(max,max),c(cum_R[which(ii==(compt-1))]*100,R_cum_child_left*100))
			
		}
		
		# If the node is a leaf #
		# if(obj$frame$var[i]=='<leaf>')
		# {
			# Add number of objects per node
			# This leaf is at which level of the hierarchy ?
		#	si<-floor(log(nodes[i])/log(2))+1
			# At this level, which position from the left does it occupy ?
		#	posi<-nodes[i]-2^(si-1)+1
			
		#	text(,0,paste('n= '.arbre$frame[i,2]),cex=1,offset=0.5,pos=1)
		#}
	}
	
	
	# Add adjusted R squared and other interesting values #
	
	len <- dim(obj$cptable)[1]
	exp1<- expression(R^2)
	bbb<-signif(1-obj$cptable[len, 3],digits=3)
	if(R2A)
	{
		R2adj<-R2aGDF(MRT(obj),T=T,X=X,tau_const=0.6,...)
		foot0 <- paste('R2a: ',round(R2adj*100,2),'%')
		foot1 <- paste('R2 : ',bbb*100,'%')
		foot2 <- paste("  Error : ", signif(obj$cptable[len, 3], 
                digits=3))
		foot <- paste(foot0,foot1,foot2, "  CV Error : ", signif(obj$cptable[len, 
                  4], digits=3), "  SE : ", signif(obj$cptable[len, 
                  5], digits=3))
         mtext(foot, side = 1, line = 3.5, cex = 0.85)
    }
    
    if(!R2A)
	{
		foot1 <- paste('R2 : ',bbb*100,'%')
		foot2 <- paste("  Error : ", signif(obj$cptable[len, 3], 
                digits=3))
		foot <- paste(foot1,foot2, "  CV Error : ", signif(obj$cptable[len, 
                  4], digits=3), "  SE : ", signif(obj$cptable[len, 
                  5], digits=3))
         mtext(foot, side = 1, line = 3.5, cex = 0.85)
    }
            
   
	
	
}


### -------------------------------- ###
###         mat_lables_fct           ###
### -------------------------------- ###

# Le nŽcessaire pour construire 'mat_labels'
mat_labels_fct<-function(obj)
{
	split<-which(obj$frame$var!='<leaf>')
	
	Labels<-labels(obj,pretty=pretty)[-1]
	
	# Les splits en ordre ...
	splits<-as.numeric(row.names(obj$frame)[split])
	
	# Ajout des colonnes d'enfants de gauche et droite
	D_child<-2*splits+1
	G_child<-2*splits
	
	mat_splits<-cbind(splits,D_child,G_child)
	
	# Get the imaginary children out
	for(i in 1:nrow(mat_splits))
	{
		for (j in 2:3)
		{
			if(!any(splits==mat_splits[i,j]))
			{
				mat_splits[i,j]=0
			}
		}
	}


	# Count the number of relatives on the left
	G_relatives<-mat.or.vec(nrow(mat_splits),nrow(mat_splits))
	mat_splits_relatives<-cbind(mat_splits,G_relatives)
	
	
	
	for(i in 1:nrow(mat_splits)) # Pour tout les noeuds
	{
		temp<-mat_splits_relatives[i,3] # the first left child
		for(j in 1:nrow(mat_splits)) # Pour tout les noeuds qui pourraient tre parents
		{
			if(any(mat_splits_relatives[j,1]==temp)) # On retrouve ce noeuds dans les parents de gauche
			{
				mat_splits_relatives[j,i+3]<-1
				temp<-c(temp,mat_splits_relatives[j,c(2,3)])
				# retire les 0
				temp<-unique(temp)
				
			}
		}
	}
	

	# How many much left relatives ? :o)
	
	if(nrow(mat_splits_relatives)==1) 
	{
		nrelatives=0
	}
	else 
	{
		nrelatives<-apply(mat_splits_relatives[,-c(1:3)],FUN=sum,MARGIN=2)
	}
	# Caculating the positions of the labels
	
	# Initialisation
	pos_labels<-mat.or.vec(length(Labels),2)
	
	compt<-1
	for(i in 1:nrow(mat_splits)) # Pour tout les noeuds
	{
	
		# First empty space
		zero<-which(pos_labels[,1]==0)[1]
		
		# Left child : the first empty space
		pos_labels[zero,1]<-compt
		compt<-compt+1
		pos_labels[zero,2]<-1
		
		
		# Right child : pos Left child + 2xnb relatives on left + 1
		pos_labels[zero+2*nrelatives[i]+1,1]<-compt
		compt<-compt+1
		pos_labels[zero+2*nrelatives[i]+1,2]<-2
		
	}
	
	# Buiding labels.. the real ones
	Labels_ordre<-Labels[order(pos_labels[,1])]
	
	# matrix()
	Labels<-matrix(Labels_ordre,nrow(mat_splits),2,byrow=TRUE)
	Labels<-cbind(splits, Labels)
	return(Labels)
}

# ------------------------------------------------------------ #
# this function makes the summary output for the MRT analsyis. #
# ------------------------------------------------------------ #
	
summary.MRT<-function(object,IndvalPART=TRUE,IndvalNODE=TRUE,...)
{
	# IndvalPART : add indval species for final partition
	# IndvalNODE : add indval species for each node
	
	if(object$LABELS==FALSE) mat_labels<-mat_labels_fct(object$obj)
	if(!object$LABELS==FALSE) mat_labels<-object$LABELS
	

	cat('Portion (%) of deviance explained by species for every particular node','\n','\n')
	
	for(i in 1:length(object$nodes))
	{
		cat('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
		# Print node number, it's complexity, it's discriminant species, it's objects
		
		# Node number + complexity
		cat('                    --- Node',object$nodes[i],'---\n','Complexity(R2)',sum(object$R2[i,])*100,'\n',mat_labels[i,-1],'\n\n')
		
		# Discriminant species #
		cat('~ Discriminant species :','\n')
		
		# Create table of results #
		table<-rbind(object$pourct[i,object$pourct[i,]>object$percent],object$MOYs[i*2-1,object$pourct[i,]>object$percent],object$MOYs[i*2,object$pourct[i,]>object$percent])
		rownames(table)<-c('% of expl. deviance','Mean on the left','Mean on the right')
	
		print(table)
		
		# Build table of indval if 
		if(IndvalNODE)
		{
			# INDVAL species #
			
			LWHEREnode<-object$LWHERE[[i]]
			RWHEREnode<-object$RWHERE[[i]]
			Ynode<-object$obj$y[c(LWHEREnode,RWHEREnode),]
		
			if (any(apply(Ynode>0,2,sum)==0)) {
    		cat('eliminating null columns\n')
    		Ynode <- Ynode[,apply(Ynode>0,2,sum)>0]
			}
			
			clustnode<-c(mat.or.vec(length(LWHEREnode),1)+1,mat.or.vec(length(RWHEREnode),1)+2)

			INDVALnode<-indval(Ynode,clustering=clustnode,numitr=1000)

			cat('\n','~ INDVAL species for this node: : left is 1, right is 2','\n',sep='')
			summary(INDVALnode, p=0.05, type='short', ...)
		}
		
		cat('\n')
		
		cat('\n')
	}
	
	# Objects in each leaf
	cat('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	cat('               --- Final partition ','---\n','\n')
	for(i in 1:length(object$obj$frame$var))
	{
		if(object$obj$frame$var[i]=='<leaf>')
		{
			cat('Sites in the leaf #',i,'\n')
			print(row.names(object$obj$y)[object$obj$where==i])
			cat('\n','\n')
		
		}
	}
	
	# Build table of indval if 
	if(IndvalPART)
	{
		# INDVAL species #
		cat('~ INDVAL species of the final partition:','\n')

		INDVALpart<-indval(object$obj$y,clustering=object$obj$where,numitr=1000)
		summary(INDVALpart, p=0.05, type='short',...)
		
		# Correspondance between $where and numbers given by indval
		k=1
		cat('\n')
		cat('~ Corresponding MRT leaf number for Indval cluster number:','\n')
		ordered_where<-sort(unique(object$obj$where))
		for(j in 1:length(object$obj$frame$var))
		{ 
			
			if(object$obj$frame$var[j]=='<leaf>')
			{
				cat('MRT leaf #',ordered_where[k],'is Indval cluster no.',k,'\n')
				k=k+1
		
			}
		}
	}
}

# This function takes Lwhere (list of left objects in every node) and Rwhere (list of rigth object in every node), and the labels from mat_labels, to create what is described in the function.


code.tree <- function(obj,LWHERE,RWHERE,LABELS)
# This function codes the tree.
# Each coding veriable in matrix 'code' represents one of the splits of the MRT.
# The objects on the left and right of the split have weights corresponding to 
# the inverse of the number of observations in each group.
# The objects not concerned with a split have a weight of zero.
# As a consequence, each column sums to 0.
#
#	   Marie-HŽlne Ouellette, August 2009
#      Modified from Pierre Legendre, August 2009
#      
{
	
	n <- nrow(obj$y)

	ifelse(length(LWHERE)==1,nvar<-1,nvar <- nrow(LABELS))
	
	code <- matrix(0,n,nvar)

	
	# For all splits, create code

	for(i in 1:length(LWHERE))
	{
		vec1 <- LWHERE[[i]]
		
		n1 <- length(vec1)
		vec2 <- RWHERE[[i]]
		
		n2 <- length(vec2)
		code[vec1,i] <- 1/n1
		
		code[vec2,i] <- -1/n2
		
	}

	# Cut the crap out (< > =) of LABELS
	if(length(LWHERE)==1)LABELS<-as.data.frame(t(LABELS))
	colnames(code) = sedit(LABELS[,2],c('<','>','=',' '),c('','','',''))
	return(code)
}

# ----------------------- #
# This function takes a coded tree from code.tree and builds the TRiplot as suggested by PL
# ----------------------- #

#plot_triplot<-function(obj,LWHERE,RWHERE,R2,Cex,widthtriordi,heighttriordi,test.F=test.F,silent=silent,LABELS,...)
#{
#	require(rdaTest)
#	group<-obj$where
#	splits<-which(obj$frame[,1]!="<leaf>")
#	ii<-ORDER_NODES(obj,R2,splits)
#	if(LABELS==FALSE) LABELS_NOTINORDER<-mat_labels_fct(obj)
#	if(!LABELS==FALSE) LABELS_NOTINORDER<-LABELS
#	LABELS<-LABELS_NOTINORDER[ii,]
	
#	#dummy = model.matrix(~ group)       # group4 was eliminated by model.matrix
#	#grouplast = dummy[,1] - apply(dummy[,-1],1,sum)
#	#dummy2 = cbind(grouplast,dummy[,-1])   # Put back group4 into the matrix
	
#	coded.tree = code.tree(obj,LWHERE,RWHERE,LABELS)
#	apply(coded.tree,2,sum)  # All columns sum to 0
	
#	# Change row names to site membership for plot
#	rownames(obj$y)<-obj$where
#	rownames(coded.tree)<-obj$where
#	rda.res3 = rdaTest(obj$y, coded.tree,test.F=test.F,silent=silent)
	
#	# Plot without objects
#	plot(rda.res3, graph.type="Z",cex=Cex,width=widthtriordi, height=heighttriordi,scaling=1)

	# Bimultivariate redundancy statistic (canonical R-square):
	# R-square = 0.7856865;   adjusted R-square =  0.7484146
#}

# ----------------------- #
# This function takes a coded tree from code.tree and builds the BIplot as suggested by PL
# ----------------------- #

# This function returns ii, the order of the nodes in terms of R2 of splits
ORDER_NODES<-function(obj,R2,splits)
{
	ii<-1
	if(length(splits)>1)
	
	{
		# The order shall be given by the children matrix (correction for when children explain more variance than parent)
	
		parent<-1 # the
		nodes2<-1
		pos_nodes<-1
		pos_parent<-1
		
		R2_node<-rowSums(R2)
		nodes_name<-as.numeric(row.names(obj$frame))[splits]
	
		for(i in 1:(length(splits)-1))
		{
			# For this node, which children are nodes also (update nodes with children) ?
			if(any(as.numeric(row.names(obj$frame))[splits]==parent*2))
			{
				nodes2<-c(nodes2,parent*2)
				pos_nodes<-c(pos_nodes,which(nodes_name==parent*2))
			}
		
			if(any(as.numeric(row.names(obj$frame))[splits]==(parent*2+1)))
			{
				nodes2<-c(nodes2,parent*2+1)
				pos_nodes<-c(pos_nodes,which(nodes_name==(parent*2+1)))
			}
		
			# Get this parent (node) out of nodes and pos_nodes
			pos_nodes<-pos_nodes[-which(nodes2==parent)]
			pos_nodes
			nodes2<-nodes2[-which(nodes2==parent)]
			nodes2

		
			# Get position in R2 and splits of the nodes
		
			# Of these, which has the smallest R2 ? This is our next parent
			parent<-nodes_name[which(R2_node==max(R2_node[pos_nodes]))]
			parent
			pos_parent<-which(nodes_name==parent)
			pos_parent
		
			# Stock it's position as it is the next node
		
			ii<-c(ii,which(R2_node==max(R2_node[pos_nodes])))
	
		}
	}
return(ii)
}


   #########################################################
   #### Following functions needed for interaction plot  ###
   #########################################################

"lab" <-function(object, digits=4, minlength=1, pretty,
                  collapse=TRUE, ...) {
    if (missing(minlength) && !missing(pretty)) {
    if (is.null(pretty)) minlength <-1
    else if (is.logical(pretty)) {
        if (pretty) minlength <- 4
        else        minlength <- 0
        }
    else minlength <- 0
    }

    ff <- object$frame
    n  <- nrow(ff)
    if (n==1) return("root")  #special case of no splits

    is.leaf <- (ff$var == "<leaf>")
    whichrow <- !is.leaf
    vnames <- ff$var[whichrow]  #the variable names for the primary splits

    index <- cumsum(c(1, ff$ncompete + ff$nsurrogate + 1*(!is.leaf)))
    irow  <- index[c(whichrow, FALSE)]     #we only care about the primary split
    ncat  <- object$splits[irow, 2]

    # Now to work: first create labels for the left and right splits,
    #  but not for leaves of course
    #
    lsplit <- rsplit <- vector(mode='character', length= length(irow))

    if (any(ncat <2)) {  # any continuous vars ?
    jrow <- irow[ncat <2]
    cutpoint <- formatg(object$splits[jrow,4], digits)
    temp1 <- (ifelse(ncat<0, "< ", ">="))[ncat <2]
    temp2 <- (ifelse(ncat<0, ">=", "< "))[ncat <2]
    lsplit[ncat<2] <- paste(temp1, cutpoint, sep='')
    rsplit[ncat<2] <- paste(temp2, cutpoint, sep='')
    }

    if (any(ncat >1)) { # any categorical variables ?
    xlevels <- attr(object, 'xlevels')
    #
    # jrow will be the row numbers of factors within lsplit and rsplit
    # crow the row number in "csplit"
    # and cindex the index on the "xlevels" list
    #
    jrow <- (seq(along=ncat))[ncat>1]
    crow <- object$splits[irow[ncat>1],4]    #row number in csplit
    cindex <- (match(vnames, names(xlevels)))[ncat >1]

    # Now, abbreviate the levels
    if (minlength ==1) {
        if (any(ncat>52))
        warning(paste("More than 52 levels in a predicting factor,",
                  "truncated for printout"))
        xlevels <- lapply(xlevels,
                   function(z) {
                   k <- length(z)
                   k <- pmin(1:k, 52)
                   c(letters, LETTERS)[k]
                   })
        }
    else if (minlength >1)
        xlevels <- lapply(xlevels, abbreviate, minlength, ...)

    # Now tuck in the labels
    # I'll let some other clever person vectorize this
    for (i in 1:(length(jrow))) {
        j <- jrow[i]
        splits <- object$csplit[crow[i],]
        # splits will contain 1=left, 2=right, 3= neither
        ltemp <- (1:length(splits))[splits== 1]
        rtemp <- (1:length(splits))[splits== 3]
        if (minlength==1) {
        lsplit[j] <- paste((xlevels[[cindex[i]]])[ltemp], collapse='')
        rsplit[j] <- paste((xlevels[[cindex[i]]])[rtemp], collapse='')
        }
        else {
        lsplit[j] <-paste((xlevels[[cindex[i]]])[ltemp], collapse=',')
        rsplit[j] <-paste((xlevels[[cindex[i]]])[rtemp], collapse=',')
        }
        }
    }

    if (!collapse) {  #called by no routines that I know of
    ltemp <- rtemp <- rep("<leaf>", n)
    ltemp[whichrow] <- lsplit
    rtemp[whichrow] <- rsplit
    return(cbind(ltemp, rtemp))
    }

    lsplit <- paste(ifelse(ncat<2, "", "="), lsplit, sep='')
    rsplit <- paste(ifelse(ncat<2, "", "="), rsplit, sep='')

    #
    # Now match them up to node numbers
    #   The output will have one label per row of object$frame, each
    #   corresponding the the line segement joining this node to its parent
    #
    varname <- (as.character(vnames))
    node <- as.numeric(row.names(ff))
    parent <- match(node %/% 2, node[whichrow])
    odd <- (as.logical(node %%2))

    labels <- vector('character', length=n)
    labels[odd] <- varname[parent[odd]]
    labels[!odd]<- varname[parent[!odd]]
    #labels[1] <- "root"
    labels<-labels[-1]
    labels
    }

"path" <-
function(tree, nodes, pretty = 0, print.it = TRUE)
{
        if(!inherits(tree, "rpart"))
                stop("Not legitimate tree")
        splits <- lab(tree, pretty = pretty)
        frame <- tree$frame
        n <- row.names(frame)
        node <- as.numeric(n)
        which <- descendants(node)      #ancestors are columns
        path <- list()
        if(missing(nodes)) {
                xy <- rpartco(tree)
                while(length(i <- identify(xy, n = 1, plot = FALSE)) > 0) {
                        path[[n[i]]] <- path.i <- splits[which[, i]]
                        if(print.it) {
                                cat("\n", "node number:", n[i], "\n")
                                cat(paste("  ", path.i), sep = "\n")
                        }
                }
        }
        else {

                if(length(nodes <- node.match(nodes, node)) == 0)
                        return(invisible())
                for(i in nodes)
                       { path[[n[i]]] <- path.i <- splits[which[, i]]
            if(print.it) {
                                cat("\n", "node number:", n[i], "\n")
                                cat(paste("  ", path.i), sep = "\n")
                                }
                       }
        }
        invisible(path)
	  p<-path.i
}


"rpartco" <- function(tree, parms)
    {
     uniform <- parms$uniform
     nspace <- parms$nspace
     nbranch <- parms$nbranch
     minbranch <- parms$minbranch

    frame <- tree$frame
    method <- tree$method
    node <- as.numeric(row.names(frame))
    depth <- tree.depth(node)
    is.leaf <- (frame$var == '<leaf>')

    if(uniform) y <- (1 + max(depth) -depth) / max(depth,4)
    else {                    #make y- (parent y) = change in deviance
    y <- dev <- frame$dev
        temp <- split(seq(node), depth)     #depth 0 nodes, then 1, then ...
        parent <- match(floor(node/2), node)
        sibling <- match(ifelse(node %% 2, node - 1, node + 1), node)

    # assign the depths
        for(i in temp[-1]) {
        temp2 <- dev[parent[i]] - (dev[i] + dev[sibling[i]])
            y[i] <- y[parent[i]] - temp2
        }
    #
    # For some problems, classification & loss matrices in particular
    #   the gain from a split may be 0.  This is ugly on the plot.
    # Hence the "fudge" factor of  .3* the average step
    #
    fudge <-  minbranch * diff(range(y)) / max(depth)
        for(i in temp[-1]) {
        temp2 <- dev[parent[i]] - (dev[i] + dev[sibling[i]])
        haskids <- !(is.leaf[i] & is.leaf[sibling[i]])
        y[i] <- y[parent[i]] - ifelse(temp2<=fudge & haskids, fudge, temp2)
        }
    y <- y / (max(y))
        }

    # Now compute the x coordinates, by spacing out the leaves and then
    #   filling in
    x   <-  double(length(node))         #allocate, then fill it in below
    x[is.leaf] <- seq(sum(is.leaf))      # leaves at 1, 2, 3, ....
    left.child <- match(node * 2, node)
    right.child <- match(node * 2 + 1, node)

    # temp is a list of non-is.leaf, by depth
    temp <- split(seq(node)[!is.leaf], depth[!is.leaf])
    for(i in rev(temp))
            x[i] <- 0.5 * (x[left.child[i]] + x[right.child[i]])

    if (nspace < 0) return(list(x=x, y=y))

    #
    # Now we get fancy, and try to do overlapping
    #
    #  The basic algorithm is, at each node:
    #      1: get the left & right edges, by depth, for the left and
    #           right sons, of the x-coordinate spacing.
    #      2: find the minimal free spacing.  If this is >0, slide the
    #           right hand son over to the left
    #      3: report the left & right extents of the new tree up to the
    #           parent
    #   A way to visualize steps 1 and 2 is to imagine, for a given node,
    #      that the left son, with all its descendants, is drawn on a
    #      slab of wood.  The left & right edges, per level, give the
    #      width of this board.  (The board is not a rectangle, it has
    #      'stair step' edges). Do the same for the right son.  Now
    #      insert some spacers, one per level, and slide right hand
    #      board over until they touch.  Glue the boards and spacer
    #      together at that point.
    #
    #  If a node has children, its 'space' is considered to extend left
    #    and right by the amount "nspace", which accounts for space
    #    used by the arcs from this node to its children.  For
    #    horseshoe connections nspace usually is 1.
    #
    #  To make it global for a recursive function, the x coordinate list
    #    is written into frame 0.
    #
    compress <- function(me, depth) {
        lson <- me +1
    x <- x
    if (is.leaf[lson]) left <- list(left=x[lson], right=x[lson],
                        depth=depth+1, sons=lson)
        else               left <- compress(me+1, depth+1)

        rson <- me + 1 + length(left$sons)        #index of right son
    if (is.leaf[rson]) right<- list(left=x[rson], right=x[rson],
                        depth=depth+1, sons=rson)
    else               right<- compress(rson, depth+1)

    maxd <- max(left$depth, right$depth) - depth
        mind <- min(left$depth, right$depth) - depth

    # Find the smallest distance between the two subtrees
    #   But only over depths that they have in common
    # 1 is a minimum distance allowed
    slide <- min(right$left[1:mind] - left$right[1:mind]) -1
    if (slide >0) { # slide the right hand node to the left
        x[right$sons] <- x[right$sons] - slide;
        x[me] <- (x[right$sons[1]] + x[left$sons[1]])/2
#       assign("x", x)
            x <<- x
        }
    else slide <- 0

    # report back
        if (left$depth > right$depth) {
        templ <- left$left
            tempr <- left$right
            tempr[1:mind] <- pmax(tempr[1:mind], right$right -slide)
        }
        else {
        templ <- right$left  - slide
        tempr <- right$right - slide
        templ[1:mind] <- pmin(templ[1:mind], left$left)
        }

    list(left = c(x[me]- nspace*(x[me] -x[lson]), templ),
         right= c(x[me]- nspace*(x[me] -x[rson]), tempr),
         depth= maxd+ depth, sons=c(me, left$sons, right$sons))
    }
    compress(1, 1)
    list(x = x, y = y)
}

"node.match" <-
function(nodes, nodelist, leaves, print.it = TRUE)
{
    node.index <- match(nodes, nodelist, nomatch = 0)
    bad <- nodes[node.index == 0]
    if(length(bad) > 0 & print.it)
        warning(paste("supplied nodes", paste(bad, collapse = ","),
                      "are not in this tree"))
    good <- nodes[node.index > 0]
    if(!missing(leaves) && any(leaves <- leaves[node.index])) {
        warning(paste("supplied nodes",
                      paste(good[leaves], collapse = ","), "are leaves"))
        node.index[node.index > 0][!leaves]
    }
    else node.index[node.index > 0]
}


