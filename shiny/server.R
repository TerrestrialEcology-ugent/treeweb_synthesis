############################################
### A shiny app to explore desirability  ###
############################################

#Below is the server code

library(shiny)
library(plyr)
library(ggplot2)
#library(ggtern)
#library(viridis)


#setwd("Documents/PostDoc_Ghent/Synthesis_stuff/Shiny/multi_desire/")
#helper functions
revv <- function(vec){
  return(- vec + (max(vec) + min(vec)))
}

desirab <- function(mat,direction,importance){
  for(i in 1:ncol(mat)){
    if(direction[i] == -1){
      mat[,i] <- revv(mat[,i])
    }
  }
  out <- apply(mat,1,function(x) sum(importance * x) / sum(importance))
  return(out)
}

dd_mat <- function(ypred,direction,importance,probs=c(0.25,0.5,0.75),type){
  tmp <- apply(ypred,1,function(mat) desirab(mat,direction,importance))
  out <- adply(tmp,1,quantile,probs,na.rm=TRUE)
  names(out)[2:4] <- c("LCI","Med","UCI")
  if(type=="div"){
    out <- cbind(out,X_div)
    out$divF <- factor(out$specrich,labels=c(1,2,3))
    out$fragmF <- factor(out$fragm,labels=c("High","Medium","Low"))
  }
  if(type=="comp"){
    out$speccomb <- rep(comp,each=3)
    out$fragm <- rep(c(-2.42,0,1.69),7)
  }
  if(type=="kir"){
    out <- cbind(out,X_kir)
  }
  return(out)
}

#load data
#the model prediction
pred_div <- readRDS("data/prediction_diversity.rds")
#pred_comp <- readRDS("data/prediction_comp.rds")
#pred_kir <- readRDS("data/prediction_kirwan.rds")
#the newdata matrix
X_div <- read.table("data/X_div.csv",head=TRUE)
#X_comp <- read.table("data/X_comp.csv",head=TRUE)
#X_kir <- read.table("data/X_kir.csv",head=TRUE)

#load the standard direction
dirr <- read.table("data/dirr.csv",head=TRUE)
dirr <- dirr[,1]

comp <- c("qrob_qrub","fsyl","fsyl_qrub","all","fsyl_qrob","qrub","qrob")

# Define server 
shinyServer(function(input, output) {
  #define the importance scores based on the choice of the user (default, all equal or manual)
  importance <- reactive({
    if(input$type == "default"){
      tmp <- c(7.5,6.5,6.5,6.5,6,6,9.5,7,4.66666666666667,6.5,6.75,
               9.33333333333333,9.33333333333333,4,7.25,4.5,2.33333333333333,
               3,3,4,4,0.666666666666667,0.666666666666667,0.666666666666667)
    }
    
  
  if(input$type == "equal"){
    tmp <- rep(1,24)
  }
  if(input$type == "manual"){
    tmp <- c(input$C_stock,input$BS,input$pH,input$CN,input$P,input$Decomp,input$Biomass,input$Cover,
             input$Veg_div,input$GLI,input$LAI,input$P_germ,input$Seed_biom,input$Arth_div,input$Herbivory,
             input$Predation,input$Fit_spider,input$Spider_diet,input$Bird_smi,input$Breed_succ,input$Bird_div,input$Egg_vol,
             input$Egg_bact,input$Egg_IgY)
  }
  if(input$type == "drop"){
    tmp <- c(7.5,6.5,6.5,6.5,6,6,9.5,7,4.66666666666667,6.5,6.75,
             9.33333333333333,9.33333333333333,4,7.25,4.5,2.33333333333333,
             3,3,4,4,0.666666666666667,0.666666666666667,0.666666666666667)
      
    tmp[!c(input$C_stock2,input$BS2,input$pH2,input$CN2,input$P2,input$Decomp2,input$Biomass2,input$Cover2,
           input$Veg_div2,input$GLI2,input$LAI2,input$P_germ2,input$Seed_biom2,input$Arth_div2,input$Herbivory2,
           input$Predation2,input$Fit_spider2,input$Spider_diet2,input$Bird_smi2,input$Breed_succ2,input$Bird_div2,
           input$Egg_vol2,input$Egg_bact2,input$Egg_IgY2)] <- 0
    
   
  }
    return(tmp)
  })


  
  #compute desirability for each model
  d_div <- reactive({dd_mat(pred_div,dirr,importance(),type="div")})
  #d_comp <- reactive({dd_mat(pred_comp,dirr,importance(),type="comp")})
  #d_kir <- reactive({dd_mat(pred_kir,dirr,importance(),type="kir")})
  #the plots
  gg_div <- reactive({ggplot(d_div(),aes(x=specrich,y=Med,group=fragmF)) +
    geom_ribbon(aes(ymin=LCI,ymax=UCI,fill=fragmF),alpha=0.2) +
    geom_path(aes(color=fragmF)) +
    scale_x_continuous(breaks = c(-0.96,0.45,1.88),labels = 1:3) +
    labs(fill="Fragmentation intensity",color="Fragmentation intensity",x="Tree species richness",y="Desirability score")})
  
  
  # gg_kir <- reactive({ggtern(d_kir(),aes(fsyl,qrob,qrub))+
  #   theme_bw()+
  #   geom_tri_tern(bins=4,fun=mean,aes(value=Med,fill=..stat..)) +
  #   stat_tri_tern(bins=4,fun=mean,geom="text",aes(value=Med,
  #                                                 label=sprintf("%.2f",..stat..)),
  #                 color="lightgrey",centroid=TRUE) +
  #   scale_fill_viridis(option="C") +
  #   theme(legend.position = "none")})
  # 
  # gg_comp <- reactive({ggplot(d_comp(),aes(x=fragm,y=Med,group=speccomb)) +
  #   geom_ribbon(aes(ymin=LCI,ymax=UCI,fill=speccomb),alpha=0.1) +
  #   geom_path(aes(color=speccomb),size=2) +
  #   labs(fill="Species\ncomposition",color="Species\ncomposition",x="Fragmentation intensity (from high to low)",y="Desirability score")})
  # 
  output$div_f_plot <- renderPlot({
    
    gg_div()
    
  })
  
  
  #output$kir_f_plot <- renderPlot({
    
  #  gg_kir()
    
  #})
  
  
  #output$comp_f_plot <- renderPlot({
  
  #  gg_comp()
    
  #})
  
  output$dlPlot <- downloadHandler(
    filename = function() paste0("desire_div_",Sys.Date(),".png"),
    content = function(file) {
      tmp <- grid.arrange(gg_div(),ncol=1,bottom=paste("Weights :",paste(round(importance()[1:12],2),collapse = " "),"\n",paste(round(importance()[13:24],2),collapse = " ")))
      ggplot2::ggsave(file,plot=tmp,width=15,height=10,units="cm")
    },
    contentType = "image/png"
  )
  
  output$debug <- renderText(importance())
  
})
