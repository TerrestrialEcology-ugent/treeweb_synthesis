##########################

# helper functions for R script

# on TREEWEB synthesis

##############

##### compute bayesian R2 for single function
bayes_R2 <- function(obs,pred,summary = TRUE){
  e <- -1 * sweep(pred,2,obs)
  var_ypred <- apply(pred,1,var)
  var_e <- apply(e,1,var)  
  if(summary){
    return(quantile(var_ypred / (var_ypred + var_e),probs=c(0.1,0.5,0.9)))
  }
  else{
    return(var_ypred / (var_ypred + var_e))    
  }
}
#apply this across the functions
#sapply(1:K,function(k) bayes_R2(dat_std[,(k+1)],ypred[,,k]))


#### function to generate positive definite matrix (covariance)
Posdef <- function (n, ev = runif(n, 0, 10)) 
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

##### function to help initialization
init_fn <- function(){
  list(L_omega = chol(Posdef(K)))
}

######## function to plot ternary graphs
make_ind_gg <- function(predicted, fun = "Predation",fragm="Low", pal=viridis(10)){
  subs <- subset(predicted,Function==fun & Fragm==fragm)
  out <-   ggtern(subs,aes(fsyl,qrob,qrub))+
    theme_bw()+
    theme_nomask() +
    #facet_wrap(~Function) +
    
    geom_tri_tern(bins=4,fun=mean,aes(value=Med_real,fill=..stat..)) +
    stat_tri_tern(bins=4,fun=mean,geom="text",aes(value=Med_real,
                                                  label=sprintf("%.2f",..stat..)),
                  color="darkorange3",size=3,centroid=TRUE) +
    annotate(geom = "text",x=0.5,y=0.9,z=-0.3,label=fun) +
    scale_fill_gradientn(colours = pal) +
    theme(legend.position = "none")
  return(out)
}

###### function to make custom correlation graphs
##### modifying code from corrr::network_plot

network_df <- function(rdf,
                       min_cor = .30,
                       legend = TRUE,
                       colours = c("indianred2", "white", "skyblue1"),
                       repel = TRUE,
                       curved = TRUE,
                       colors) {
  
  require(dplyr)
  require(ggplot2)
  
  if (min_cor < 0 || min_cor > 1) {
    stop ("min_cor must be a value ranging from zero to one.")
  }
  
  if (!missing(colors))
    colours <- colors
  
  rdf %>%
    as_matrix(diagonal = 1) -> rdf
  distance <- sign(rdf) * (1 - abs(rdf))
  
  # fixed position of variables in 4 different columns
  points <- data.frame(x = c(rep(-0.5,6),rep(-0.17,7),rep(0.17,5),rep(0.5,6)),y=c(seq(-0.5,0.5,length=6),
                                                                                  seq(-0.5,0.5,length=7),
                                                                                  seq(-0.5,0.5,length=5),
                                                                                  seq(-0.5,0.5,length=6)),id=nice_name)
  #points$id <- rownames(points)
  
  # Create a proximity matrix of the paths to be plotted.
  proximity <- abs(rdf)
  proximity[upper.tri(proximity)] <- NA
  diag(proximity) <- NA
  proximity[proximity < min_cor] <- NA
  
  # Produce a data frame of data needed for plotting the paths.
  n_paths <- sum(!is.na(proximity))
  paths <- matrix(nrow = n_paths, ncol = 6) %>% data.frame()
  colnames(paths) <- c("x", "y", "xend", "yend", "proximity", "sign")
  path <- 1
  for(row in 1:nrow(proximity)) {
    for(col in 1:ncol(proximity)) {
      path_proximity <- proximity[row, col]
      if (!is.na(path_proximity)) {
        path_sign <- sign(distance[row, col])
        x    <- points$x[row]
        y    <- points$y[row]
        xend <- points$x[col]
        yend <- points$y[col]
        paths[path, ] <- c(x, y, xend, yend, path_proximity, path_sign)
        path <- path + 1
      }
    }
  }
  
  plot_ <- list(
    # For plotting paths
    if (curved) geom_curve(data = paths,
                           aes(x = x, y = y, xend = xend, yend = yend,
                               alpha = proximity, size = proximity,
                               colour = proximity*sign)), 
    if (!curved) geom_segment(data = paths,
                              aes(x = x, y = y, xend = xend, yend = yend,
                                  alpha = proximity, size = proximity,
                                  colour = proximity*sign)), 
    scale_alpha(limits = c(0, 1)),
    scale_size(limits = c(0, 1)),
    scale_colour_gradientn(limits = c(-1, 1), colors = colours),
    # Plot the points
    #geom_point(data = points,
    #           aes(x, y),
    #           size = 3, shape = 19, colour = "white"),
    # Plot variable labels
    if (repel) ggrepel::geom_text_repel(data = points,
                                        aes(x, y, label = id),
                                        fontface = 'bold', size = 5,
                                        segment.size = 0.0,
                                        segment.color = "white"),
    if (!repel) geom_text(data = points,
                          aes(x, y, label = id),
                          fontface = 'bold', size = 5,color="grey90"),
    # expand the axes to add space for curves
    expand_limits(x = c(min(points$x) - .1,
                        max(points$x) + .1),
                  y = c(min(points$y) - .1,
                        max(points$y) + .1)
    ),
    # Theme and legends
    theme_void(),
    guides(size = "none", alpha = "none"),
    if (legend)  labs(colour = NULL),
    if (!legend) theme(legend.position = "none")
  )
  
  ggplot() + plot_ + theme(plot.background = element_rect(fill="grey70"))
  
}
