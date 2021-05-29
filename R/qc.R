.get_mtRNA_genes = function(species = c("mmusculus", "hsapiens")){
  MT_RNA = c('mt-Nd1', 'mt-Nd2', 'mt-Nd3', 'mt-Nd4','mt-Nd4l', 'mt-Nd5', 'mt-Nd6', 'mt-Co1', 'mt-Co2', 'mt-Co3', 'mt-Atp6', 'mt-Atp8')

  species = match.arg(species)
  if(species == "mmusculus")
    MT_RNA = c(MT_RNA, 'mt-Cytb')
  else if(species == "hsapiens")
    MT_RNA = c(toupper(MT_RNA), 'MT-CYB')
  else{
    out = sprintf("Invalid species.\n")
    stop(out)
  }
  return(MT_RNA)
}

get_mtRNA_stats <- function(ace, by = NULL, groups_use = NULL, features_use = NULL, assay = "counts", species = c("mmusculus", "hsapiens"), metric = c("pct", "ratio", "counts")){

  require(stats)
  species = match.arg(species)
  metric = match.arg(metric)

  features_use = ACTIONet:::.preprocess_annotation_features(ace, features_use = features_use)
  mask = features_use %in% .get_mtRNA_genes(species)
  mat = assays(ace)[[assay]]
  cs_mat = ACTIONet::fastColSums(mat)
  mm = mat[mask, , drop = F]
  cs_mm = ACTIONet::fastColSums(mm)

  if(!is.null(by)){

    IDX = ACTIONet:::.get_attr_or_split_idx(ace, by, groups_use)

    if(metric == "pct"){
      frac.list = lapply(IDX, function(idx){
        m = cs_mm[idx]/cs_mat[idx]
      })
    } else if(metric == "ratio"){
      frac.list = lapply(IDX, function(idx){
        m = cs_mm[idx]/(cs_mat[idx] - cs_mm[idx])
      })

    } else{
      frac.list = lapply(IDX, function(idx){
        m = cs_mm[idx]
      })
    }
    return(frac.list)
  } else {

    if(metric == "pct"){
        frac = cs_mm/cs_mat
    } else if(metric == "ratio"){
        frac = cs_mm/(cs_mat - cs_mm)
    } else{
        frac = cs_mm
    }
    return(frac)

  }

}

plot.mtRNA.dist.by.attr <- function(
  ace,
  by,
  groups_use = NULL,
  features_use = NULL,
  assay = "counts",
  log_scale = FALSE,
  species = c("mmusculus", "hsapiens"),
  metric = c("pct", "ratio", "counts"),
  to_return = c("plot", "data"),
  palette = NULL,
  x_label = NULL,
  y_label = NULL,
  plot_title = NULL
){

  require(stats)
  to_return <- match.arg(to_return)

  frac.list = get_mtRNA_stats(
    ace = ace,
    by = by,
    groups_use = groups_use,
    features_use = features_use,
    assay = assay,
    species = species,
    metric = metric
  )

  df = lapply(1:length(frac.list), function(l) data.frame(attr = names(frac.list)[l], frac = frac.list[[l]]))
  df = do.call(rbind, df)

  if(log_scale){
    df$frac[df$frac == 0] = 1
    df$frac = log10(df$frac)
    y_label = ifelse(is.null(y_label), ylab("frac (log10)"), y_label)
  } else{
    y_label = ifelse(is.null(y_label), ylab("frac"), y_label)
  }

  if(to_return == "data")
    return(df)
  else{
    p = .plot_gg_violin(
      df,
      x=df$attr,
      y=df$frac,
      fill=df$attr,
      x_label = x_label,
      y_label = y_label,
      plot_title = plot_title,
      palette = palette
    )
    return(p)
  }
}

plot.counts.by.attr <- function(
  ace,
  attr,
  nonzero = FALSE,
  log_scale = FALSE,
  assay = "counts",
  to_return = c("plot", "data"),
  palette = NULL,
  x_label = NULL,
  y_label = NULL,
  plot_title = NULL
){

  to_return <- match.arg(to_return)

  require(ggplot2)

  IDX = ACTIONet:::.get_attr_or_split_idx(ace, attr)
  mat = assays(ace)[[assay]]

  if(nonzero == TRUE)
    mat = mat > 0

  cs_mat = ACTIONet::fastColSums(mat)
  sums.list = lapply(IDX, function(idx){
      cs_mat[idx]
  })

  df = lapply(1:length(sums.list), function(l) data.frame(attr = names(sums.list)[l], sums = sums.list[[l]]))
  df = do.call(rbind, df)

  if(log_scale){
    df$sums[df$sums == 0] = 1
    df$sums = log10(df$sums)
    y_label = ifelse(is.null(y_label), ylab("sums (log10)"), y_label)
    } else{
    y_label = ifelse(is.null(y_label), ylab("sums"), y_label)
  }

  if(to_return == "data")
    return(df)
  else{
    p = .plot_gg_violin(
      df,
      x=df$attr,
      y=df$sums,
      fill=df$attr,
      x_label = x_label,
      y_label = y_label,
      plot_title = plot_title,
      palette = palette
    )
    return(p)
  }

}

.plot_gg_violin <- function(df, x, y, fill, x_label = NULL, y_label = NULL, plot_title = NULL, palette = NULL){

  require(ggplot2)
  if(is.null(x_label))
    x_label = element_blank()
  if(is.null(y_label))
    y_label = element_blank()
  if(is.null(palette)){
    palette = ACTIONet::CPal_default[seq.int(unique(x))]
  }

  all_labels = unique(x) %>% sort

  p <- ggplot(df, aes(x=x, y=y, fill=fill)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), scale = "width") +
    scale_fill_manual(values = palette) +
    scale_x_discrete(
      labels = all_labels,
      position = "bottom"
      ) +
    labs(x = x_label, y = y_label) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      panel.border = element_rect(color = "black", fill = NA),
      legend.position = "none"
    )

  if(!is.null(plot_title))
    p <- p + ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5))

  return(p)
}
