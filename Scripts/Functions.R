# Function - Density Plot ---------------------------------------------------
densityPlot <- function( data=data, title=title, legendTitle = "Sample" ,phenoData = pheno, fill = "phenotype", x.label = element_blank(), y.label = element_blank() ) {
  require(ggplot2)
  data.plot <- reshape2::melt(data)
  data.plot$Var3 <- phenoData[,fill][match( data.plot$Var2, phenoData$sampleNames)]
  p <- ggplot2::ggplot(
    data = data.plot
    ,  aes(x = value
           , color = Var3
           , fill =Var2
    )
  ) +
    geom_density(
      alpha=0.05
    ) +
    labs(
      title=title
      , x= x.label
      , y = y.label
      , color = legendTitle
      , fill = legendTitle
    ) +
    theme_classic() + 
    guides(fill=FALSE)
  return(p)
}


# Function - Box Plot -------------------------------------------------------
ggBoxPlot <- function( data=data, title=title, legendTitle = "Sample", phenoData = pheno, fill = "phenotype", x.label = element_blank(), y.label = element_blank() ) {
  require(ggplot2)
  require(reshape2)
  data.plot <- reshape2::melt(data)
  
  
  data.plot$Var3 <- phenoData[,fill][match( data.plot$Var2, phenoData$sampleNames)]
  
  p <- ggplot2::ggplot(
    data = data.plot
    ,  aes(x = Var2
           , y = value
           , fill = Var3
    )
  ) +
    geom_boxplot(
      alpha=0.5
    ) +
    labs(
      title=title
      , x= x.label
      , y = y.label
      , fill = legendTitle
    ) +
    #scale_fill_manual(values=fill.scheme) +
    theme_classic()  +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1)
      #, plot.title = element_text(hjust = 0.5)
    )
  
  p
  return(p)
}

# Function - PCA ------------------------------------------------------------
pcaPlot <- function (data=data, title=title, shape = colnames(data), color = colnames(data), shape.legend = "Cell Type", color.legend = "Cell Type", shape.manual = F, shape.value=NULL  ) {
  require(ggplot2)
  PCA_raw <- prcomp(t(data), scale. = FALSE)
  percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
  sd_ratio <- sqrt(percentVar[2] / percentVar[1])
  
  plot.data <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2]
                          , Phenotype =colnames( data )
  )
  p<-ggplot(plot.data, 
            aes(
              x = PC1,
              y = PC2)) +
    
    geom_point(
      size = 2
      , aes(
        shape = shape
        , colour = color
      )
    ) +
    ggtitle(title) +
    xlab(paste0("PC1, VarExp: ", percentVar[1], "%")
    ) +
    ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
    labs(
      colour = color.legend
      , shape = shape.legend
      , size = color.legend
    )+
    theme(
      panel.grid.major = element_blank()
      , panel.grid.minor = element_blank()
      , panel.background = element_blank()
      , axis.line = element_line(colour = "black")
    ) +
    scale_shape_manual(values=shape.value) + 
    geom_vline(xintercept = 0, color = "black", linetype=2)  + 
    geom_hline(yintercept = 0, color = "black", linetype=2)
  
  coord_fixed(ratio = sd_ratio) 
  return(p)
}
