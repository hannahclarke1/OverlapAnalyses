Overlap_func <- function(list1, list2, universe, nadj=1, trunc, trunc_col=NULL, trunc_row=NULL, plot, plotmax=0.5,list1_title="Default", split=NULL, name_trunc=NULL, tidy_col=NULL, tidy_row=NULL, save_res = NULL, save_fig = NULL, hm_title = NULL, width_fig = NULL, height_fig=NULL) {
  
  myPal <- colorRamp2(seq(0, plotmax, length.out = 100), viridis(100))
  
  overlap <- newGOM(list1, list2, genome.size=universe)
  pval <- getMatrix(overlap, name="pval")
  
  intersection <- getMatrix(overlap,name="intersection")
  lengthsA <- sapply(getGsetA(overlap), length)
  lengthsB <- sapply(getGsetB(overlap), length) 
  
  overlap_coeff <- outer(rownames(intersection), colnames(intersection), Vectorize(function(a, b) {
    inter <- intersection[a, b]
    min_size <- min(lengthsA[a], lengthsB[b])
    inter / min_size
  }))
  
  rownames(overlap_coeff) <- rownames(intersection)
  colnames(overlap_coeff) <- colnames(intersection)
  
  if (nadj==1) {
    adj <- matrix(p.adjust(as.vector(as.matrix(pval)),method="BH"), ncol=ncol(pval), nrow=nrow(pval))
    rownames(adj) <- rownames(pval)
    colnames(adj) <- colnames(pval)
  } else if (nadj==3) {
    adj1 <- matrix(p.adjust(as.vector(as.matrix(pval[,1:2])),method="BH"), ncol=2, nrow=nrow(pval))
    adj2 <- matrix(p.adjust(as.vector(as.matrix(pval[,3:4])),method="BH"), ncol=2, nrow=nrow(pval))
    adj3 <- matrix(p.adjust(as.vector(as.matrix(pval[,5:8])),method="BH"), ncol=4, nrow=nrow(pval))
    
    rownames(adj1) <- rownames(pval)
    rownames(adj2) <- rownames(pval)
    rownames(adj3) <- rownames(pval)
    
    colnames(adj1) <- colnames(pval)[1:2]
    colnames(adj2) <- colnames(pval)[3:4]
    colnames(adj3) <- colnames(pval)[5:8]
    
    adj <- cbind(adj1,adj2,adj3)
    rm(adj1,adj2,adj3)
  }

  if (trunc==TRUE) {
    if (is.null(trunc_row)) {
      stop("trunc_row cannot be NULL when trunc==TRUE")} 
    if (is.null(trunc_col)) {
      stop("trunc_col cannot be NULL when trunc==TRUE")}
    
    if (nrow(overlap_coeff)<ncol(overlap_coeff)) {
    overlap_coeff_trunc <- overlap_coeff[, !colnames(overlap_coeff) %in% trunc_col]
    overlap_coeff_trunc <- overlap_coeff_trunc[!rownames(overlap_coeff_trunc) %in% trunc_row, ]
    adj_trunc <- adj[, !colnames(adj) %in% trunc_col]
    adj_trunc <- adj_trunc[!rownames(adj_trunc) %in% trunc_row, ] } else {
      
      overlap_coeff_trunc <- overlap_coeff[!rownames(overlap_coeff) %in% trunc_row, ]
      overlap_coeff_trunc <- overlap_coeff_trunc[, !colnames(overlap_coeff_trunc) %in% trunc_col]
      adj_trunc <- adj[!rownames(adj) %in% trunc_row, ]
      adj_trunc <- adj_trunc[, !colnames(adj_trunc) %in% trunc_col]
      
    }
    
  
  if (is.null(dim(overlap_coeff_trunc)) == TRUE) {
    overlap_coeff_trunc <- matrix(overlap_coeff_trunc,nrow = 1,dimnames = list(name_trunc,names(overlap_coeff_trunc)))
    adj_trunc <- matrix(adj_trunc,nrow = 1,dimnames = list(name_trunc,names(adj_trunc)))
  } else {
    rownames(overlap_coeff_trunc) == rownames(adj_trunc)
    colnames(overlap_coeff_trunc) == colnames(adj_trunc)}
  }
  
  if (is.null(tidy_col)==FALSE && trunc==TRUE) { 
    colnames(overlap_coeff_trunc) <- tidy_col
  } else if (is.null(tidy_col)==FALSE && trunc==FALSE) {
    colnames(overlap_coeff) <- tidy_col
  }
  
  if (is.null(tidy_row)==FALSE && trunc==TRUE) { 
    rownames(overlap_coeff_trunc) <- tidy_row
  } else if (is.null(tidy_row)==FALSE && trunc==FALSE) {
    rownames(overlap_coeff) <- tidy_row
  }
      
  
  if (plot==TRUE) {
    if (trunc==TRUE) {
      if (is.null(split)==FALSE) {
        split_c=split
      } else {split_c=rep("",ncol(overlap_coeff_trunc))}
      
        if (max(overlap_coeff_trunc)>plotmax) {
          stop("Error: Maximum value larger than given scale!")} 
      
    test <- Heatmap(overlap_coeff_trunc,
                    name = "Szymkiewicz–Simpson",
                    col = myPal,
                    cluster_rows = FALSE,
                    cluster_columns = FALSE,
                    rect_gp = gpar(col = "grey", lwd = 0.2),
                    row_title= list1_title,
                    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                    show_column_names = TRUE,
                    column_split = split_c,
                    column_names_gp = gpar(fontsize = 10),
                    show_row_names = TRUE,
                    row_names_gp = gpar(fontsize = 10),
                    row_names_side = "left",
                    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                    column_title_side = "top",
                    column_title_rot =  0,
                    cell_fun = function(j, i, x, y, w, h, fill) {
                      
                      stars <- if (adj_trunc[i, j] < 0.0001) {
                        grid.text("****", x, y, gp = gpar(fontface = "bold", fontsize=14, col="white"))
                      } else if (adj_trunc[i, j] < 0.001) {
                        grid.text("***", x, y, gp = gpar(fontface = "bold", fontsize=14, col="white"))
                      } else if (adj_trunc[i, j] < 0.01) {
                        grid.text("**", x, y, gp = gpar(fontface = "bold", fontsize=14, col="white"))
                      } else if (adj_trunc[i, j] < 0.05) {
                        grid.text("*", x, y, gp = gpar(fontface = "bold", fontsize=14, col="white"))
                      } else {
                        ""
                      }
                    })
    } else if (trunc==FALSE) {
      if (is.null(split)==FALSE) {
        split_c=split
      } else {split_c=rep("",ncol(overlap_coeff))}
      
      if (max(overlap_coeff)>plotmax) {
        stop("Error: Maximum value larger than given scale!")} 
      
      test <- Heatmap(overlap_coeff,
                      name = "Szymkiewicz–Simpson",
                      col = myPal,
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      rect_gp = gpar(col = "grey", lwd = 0.2),
                      row_title= list1_title,
                      row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                      show_column_names = TRUE,
                      column_split = split_c,
                      column_names_gp = gpar(fontsize = 10),
                      show_row_names = TRUE,
                      row_names_gp = gpar(fontsize = 10),
                      row_names_side = "left",
                      column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                      column_title_side = "top",
                      column_title_rot =  0,
                      cell_fun = function(j, i, x, y, w, h, fill) {
                        
                        stars <- if (adj[i, j] < 0.0001) {
                          grid.text("****", x, y, gp = gpar(fontface = "bold", fontsize=14, col="white"))
                        } else if (adj[i, j] < 0.001) {
                          grid.text("***", x, y, gp = gpar(fontface = "bold", fontsize=14, col="white"))
                        } else if (adj[i, j] < 0.01) {
                          grid.text("**", x, y, gp = gpar(fontface = "bold", fontsize=14, col="white"))
                        } else if (adj[i, j] < 0.05) {
                          grid.text("*", x, y, gp = gpar(fontface = "bold", fontsize=14, col="white"))
                        } else {
                          ""
                        }
                      })
      
      
    }
    
    draw(test)
    svglite(save_fig, width = width_fig, height = height_fig)
    draw(test)
   dev.off()
    
  }
  
  if (trunc==TRUE) {
    res <- list(overlap_coeff_trunc=overlap_coeff_trunc,
                       adj_trunc=adj_trunc)
  } else if (trunc==FALSE) {
    res <- list(overlap_coeff=overlap_coeff,
                       adj=adj)
  }
  
  if (is.null(save_res)==FALSE) {
    as.data.frame(res)
    write.csv(res, file=save_res)
  }
  
  
  if (trunc==TRUE) {
    return(res)
  } else if (trunc==FALSE) {
    return(res)
  }
}




