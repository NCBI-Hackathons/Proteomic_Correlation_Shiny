#' Plot the intensity profiles across several fractions
#' 
#' @description 
#' Basic plot for visualizing the intensity of individual proteins across 
#' multiple fractions.
#' 
#' @details 
#' Expects a data.frame in the long format typical of ggplot2.
#' 
#' @param long.df data.frame fit for ggplot2
#' @param x expression indicating the column name to be used for the x-axis
#' values, e.g. "fraction"
#' @param y expression indicating the column name to be used for the x-axis
#' values, one of:
#' c("raw.intensity", "MS.MS.count", "peptides.count",
#' "unique.peptides.only", "razor.and.unique.peptides") 
#' @param what indicate how the points that will be connected by lines should be
#' grouped together, e.g. c("gene_symbol", "expt_id") will draw a
#' separate line for each unique combination of gene_symbol and expt_id
#' @param color.by expression indicating the factor that will be used to assign
#'  different colors, e.g. "gene.symbol"
#' @param split.by expression indicating the factor that will be used to split
#'  the single plot into several ones, e.g. "date"
#' @param split.by.col factor that will be used to split the plot column-wise.
#' The default '.' is equivalent to 'NA'
#' @param color.palette name of RColorBrewer palette used for plotting; 
#' recommended palettes: "Set1", "Set2", "Set3", "Pastel2", "Paired", "Dark2",
#'  "Accent"; 
#' see \url{http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/figure/
#' unnamed-chunk-14-1.png} for more palettes
#' @param point.shape either a \emph{number} defining the shape for the points for all
#' groups or a \emph{string} which will be interpreted as a factor that the
#'  point shape should depend on; see
#'  \url{http://www.cookbook-r.com/Graphs/Shapes_and_line_types/} for shape types
#'  and the corresponding integers.
#' @param point.size number defining the size for the points for all
#' groups
#' @param line.type either a \emph{number} defining the line type for all
#' groups or a \emph{string} which will be interpreted as a factor that the 
#' line type should depend on; see \url{http://www.cookbook-r.com/Graphs/Shapes_and_line_types/}
#'  for line types and the corresponding integers
#' @param line.size number defining the width for the line for all
#' groups
#' @param y.lab title for the y axis; default: "value"
#' @param title title for the plot
#' @param line.smooth Boolean indicating whether a smoothened line should be displayed.
#' 
#' @return A line plot visualizing the value of interest for the gene of
#' interest across all fractions.
#' 
#' @seealso \code{\link{MQ_to_ggplot2}} 
#' 
#' @export
plot_profile <- function(long.df, x = "fraction", y = "value",
                         what = c("gene_symbol", "expt_id"), color.by = "gene_symbol",
                         split.by = NULL, split.by.col = ".", 
                         color.palette = NULL, point.shape = 16L, point.size=1L,
                         line.type = 1L, line.size = 1L, y.lab = "value",
                         title = NULL, line.smooth = FALSE){
  
  if(!( is.numeric(point.size) & is.numeric(line.size) ) ){
    stop("point.size and line.size must be numeric.")
  }
  
 # check_columns(c(what, color.by, split.by, y.lab, x, y), long.df, "long.df", "plot_profile")
  
  if(  suppressWarnings(is.na(as.numeric(line.type)) )){
    check_columns(line.type, long.df, "long.df", "plot_profile")
  }
  if( suppressWarnings(is.na(as.numeric(point.shape)) )){
    check_columns(point.shape, long.df, "long.df", "plot_profile")
  }
  ###############
  # make unique combinations of the indicated factors
  # using interaction()
  p <- ggplot(data=long.df, aes_string(x = x, y = y, 
                                       group=paste0( "interaction(",paste0( what, collapse = "," ), ")" ), 
                                       colour = color.by)) + 
    theme_bw(base_size = 16) + ylab(y.lab) 
  
  ## make points
  if( suppressWarnings(is.na(as.numeric(point.shape))) ){
    p <- p + geom_point(aes_string(shape = point.shape), size = point.size)  
  } else {
    p <- p + geom_point(shape = is.numeric(point.shape), size = point.size) 
  }
  
  ## add line
  if (line.smooth) {
   
    dat.smooth <- smooth_for_plot(long.df, what = what, x.values = x, y.values = y)
    
    if( suppressWarnings(is.na(as.numeric(line.type))) ){
      
      p <- p + geom_line(data=dat.smooth, 
                         aes_string(linetype = line.type),
                         size = line.size)
      } else {
      p <- p + geom_line(data=dat.smooth,
                         linetype = as.numeric(line.type), size = line.size)
      }
    
  } else {
    if( suppressWarnings(is.na(as.numeric(line.type)) )){
      p <- p + geom_line(data=long.df,
                         aes_string(linetype = line.type), size = line.size)
    } else {
      #print(str(long.df))
      #print(levels(interaction(long.df[what])))
      p <- p + geom_line(linetype = as.numeric(line.type), size = line.size)
    }
  }
  
  if (!is.null(split.by)) {
    p <- p + facet_grid(as.formula(paste(paste(split.by, "~"), split.by.col)))
  }
  
  if (!is.null(color.palette)) {
    n.values <- length(unique(long.df[,color.by]))
    colors.customized <- color.palette(n.values)
    p <- p + scale_colour_manual(values = colors.customized)
  }
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  return(p)
}

#' Smoothening the data for plotting
#' 
#' This calculates X-spline values relative to control points.
#' 
#' @param in.df Data.frame 
#' @param what indicate how the points that will be connected by lines should be
#' grouped together, e.g. c("gene_symbol", "expt_id") will draw a
#' separate line for each unique combination of gene_symbol and expt_id
#' @param x.values string indicating the column containing the values to be plotted
#' on the x.axis
#' @param y.values string indicating the column containing the values to be plotted on the y.axis
#' 
#' @details The data.frame will be sorted first by factor2, then factor1, then the x.values.
#' The smoothened values will be calculated using ddplyr(in.df, factor2 ~ factor1).
#' @return Data frame with smoothened values.
#' @seealso \code{\link{xspline}}, \code{\link{plot_profile}}
smooth_for_plot <- function(in.df, what = c("gene_symbol" , "expt_id"), x.values = "fraction", y.values = "value"){
  
  check_columns(c(x.values, y.values), in.df, "in.df", "smooth_for_plot()")
  
  if(!all( dim(in.df) > 0) ){
    stop(paste("The data.frame supplied to smooth_for_plot()
               does not have sufficient data. Its dimensions are:", dim(in.df) ))
  }
  
  out.df <- as.data.frame(in.df)
  # make sure that the x.values are in correct, ascending order
  out.df <- out.df[ order( out.df[x.values]), ]
  
  # calculate splines per group which is defined by the factors remaining after
  # ignoring x and y values
  grouping <- names(out.df)[! ( grepl(x.values, names(out.df)) | grepl(y.values, names(out.df)) ) ]
  
  plot.new() # this is a quirky xspline requirement...
  
  out.df <- plyr::ddply(out.df, grouping, function(x) data.frame( xspline(x[,c(x.values, y.values)], 
                                                                               shape=-0.5, draw=F) ) )
  
  names(out.df)[names(out.df) == "x"] <- x.values
  names(out.df)[names(out.df) == "y"] <- y.values
  
  if( !all( dim(out.df) > 0) ){
    stop("Something went wrong with the smoothening. Possibly because wrong factors were indicated?")
  }
  
  if( dim(out.df)[2] != dim(in.df)[2] ){
    stop(paste("The smoothening changed the columns. Original columns from in.df:", names(in.df),
         "New columns after smoothening:", names(out.df)))
  }
  
  return(out.df)
}


#' Plotting "metadata"
#'
#' @description 
#' For plotting the sum of certain values per fraction,
#' e.g., the sum of all unique peptides identified per fraction.
#' 
#' @param long.df data.frame or data.table fit for ggplot2, 
#' generated with MQ_to_cumulativePlot
#' @param x.interval number indicating the interval for labeling the fractions
#' @param color.palette name of RColorBrewer palette used for plotting; 
#' see http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/figure/
#' unnamed-chunk-14-1.png for palette names
#' @param title title for the plot
#' 
#' @return A bar plot visualizing the sum of the value of interest across all
#' fractions.
#' 
#' @seealso \code{\link{MQ_to_cumulativePlot}}
#' 
#' @export
plot_cumulativeValues <- function(long.df, x.interval = 5,
                                  color.palette = NULL,
                                  title = NULL, split.by=NULL,
                                  y.lab = NULL){
  
  l <- length(long.df$fraction)
  Y <-  "value" # names(long.df)[!grepl("fraction", names(long.df))]
  
  if( !("fraction" %in% names(long.df)) )(stop("The input file must have a variable names 'fraction'." ))
  #if( length(names(long.df)) > 3 )(stop("The input file must not have more than 2 variables." ))
  if( !is.numeric(long.df[[Y]]) )(stop("The input file must have numeric entries, at least in the column that does not contain the fractions."))
  
  p <- ggplot( data = long.df,
               aes_string(x="fraction", y = Y,
                          fill = sprintf("factor(%s)", "fraction"))) +
    geom_bar(stat="identity") +
    theme_bw(base_size = 16) + 
    scale_x_continuous(breaks=seq(x.interval, l, x.interval) ) +
    ylab(paste("sum of", Y)) +
    guides(fill=FALSE) + ylab(paste("sum of", y.lab))
  
  if (!is.null(color.palette)) {
    colors.customized <- color.palette(l)
    p <- p + scale_fill_manual(values = colors.customized)
  }else{
    p <- p + scale_fill_manual(values = rep("cornflowerblue", l)) 
  }
  
  if (!is.null(split.by)) {
    p <- p + facet_grid(as.formula(paste(split.by, "~ .")))
  }
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  return(p)
}
