#' Construct a moving block bootstrap resampling scheme
#'
#' A spatial moving block bootstrap resampling scheme is constructed, that can then
#' be used in resampling-based inference (e.g. using \code{\link[mvabund]{anova.manyglm}}).
#' Blocks are constructed as discs or tiles of observations of a user-specified size and shape.
#'
#' @param x the easting of spatial coordinate for the site.
#' @param y the northing of the spatial coordinate for the site.
#' @param block_L the size of the spatial blocks to be resampled.
#' @param nBoot the number of resamples required (defaults to 499).
#' @param Grid_space the resolution of the lattice used to sample centres of spatial blocks.
#' Defaults to a third of the resolution of \code{block_Ls}.
#' @param lookuptables.folderpath (optional) path to a directory where lookup tables can be found
#' that assign observations to spatial blocks. Such tables would considerably speed up the process.
#' @param shape shape of the spatial blocks to resample. Default is \code{disc}, but \code{square} 
#' tiles can also be specified.
#' @param ... additional arguments passed to \code{resample_blocks_by_area}.
#'
#' @details
#' A spatial moving block bootstrap resampling scheme is constructed, that can then
#' be used in resampling-based inference. A matrix of IDs is returned, with different resamples 
#' in different rows of the matrix, which can be used as input to functions designed for 
#' resampling-based inference (such as \code{\link[mvabund]{anova.manyglm}}).
#'
#' Blocks are constructed as discs or tiles of observations, whose size is controlled by \code{blockl_L},
#' these blocks are kept together in resampling. The centre of each block is chosen 
#' at random along a lattice specified via \code{Grid_space}. 
#' 
#' The most computationally intensive part of this process is working out which observations 
#' belong in which blocks. If repeated analyses are to be undertaken in the same spatial domain 
#' using the same sites, it is best to run this process only once and save the result as a set of
#' lookup tables via the function \bold{Eve to advise what function to use}.
#' Then the path to the directory containing these tables is specified via \code{lookuptables.folderpath}
#' and the whole thing runs heaps faster.
#' 
#' @return A matrix of IDs for each moving block bootstrap resample, with different
#' resamples in rows and observations in columns.
#' 
#' @author Eve Slavich <eve.slavich@@unsw.edu.au>
#' 
#' @seealso \code{\link[mvabund]{anova.manyglm}}, \code{\link{lm}}
#' @examples
#' \donttest{
#' data(Myrtaceae)
#' # fit a lm:
#' library(mvabund)
#' Myrtaceae$logrich=log(Myrtaceae$richness+1)
#' mft_richAdd = manylm(logrich~soil+poly(TMP_MAX,degree=2)+
#'                      poly(TMP_MIN,degree=2)+poly(RAIN_ANN,degree=2),
#'                                         data=Myrtaceae)
#'                                         
#' # construct a boot ID matrix: 
#' BootID = BlockBootID(x = Myrtaceae$X, y = Myrtaceae$Y, block_L = 20,
#'             nBoot = 199, Grid_space = 5)
#' anova(mft_richAdd,resamp="case",bootID=BootID)
#' }
#' 
#' @export

BlockBootID = function (x ,
                           y,
                           block_L,
                           nBoot = 499,
                           Grid_space,
                           lookuptables.folderpath = NA,
                           shape = "disc",
                           ...) {
  
  #set the spacing between sampling points
  if(missing(Grid_space)) {
    Grid_space = block_L/3
  }

    if (block_L > 0) {
      lookup_table=NA
      lookup.coords=NA
      if (is.na(lookuptables.folderpath) ==FALSE){ #check if a lookup table has been created to speed things up, if it has load it, and check the lookup table is for data with same x,y, coordinates
        load_file_if_exists(paste0(lookuptables.folderpath,"lookup_table","_L",block_L,"_grid_space_",Grid_space,"_",shape,".RData"))
        
        if(is.na(lookup.coords)==FALSE){
          if( identical ( lookup.coords$lookup.x , x ) == FALSE |  identical ( lookup.coords$lookup.y , y ) == FALSE ){
            lookup_table=NA
          }
        }
      }
      
      #Create the new sample (BlockBootID)
      new_sample = resample_blocks_by_area(
        x = x,
        y = y,
        nBoot = nBoot,
        lookup_table = lookup_table,
        block_L = block_L ,
        Grid_space = Grid_space,
        area_or_sites = "sites",
        shape = shape,
        lookuptables.folderpath = lookuptables.folderpath,
        ...
      )
      
    }
    #If block_L =0, do an iid bootstrap
    if (block_L == 0) {
      new_sample = list()
      for (i in 1:nBoot) {
        new_sample [[i]] = sample(1:length(x), size = length(x), replace = T)
      }
    }
  new_sample = t(as.data.frame(new_sample))
  rownames(new_sample) = paste0("boot.",1:nBoot)
  ;
  new_sample
}

create_moving_block_lookup = function(x, y, block_L, Grid_space, Grid, shape,...){
  
  if(missing(Grid)) {
    Grid = create_Grid(Grid_space = Grid_space, x = x , y = y)
  }
  sites = data.frame( x, y , ind = 1:length(x))
  ind = list()
  l = 1
  if(shape == "disc"){
    for ( i in 1:dim(Grid)[1]){
      ref_point_x = Grid[i,"x"]
      ref_point_y = Grid[i,"y"]
      
      # a square around the reference point of area (2BlockL)^2 
      x1 = Grid [ i , "x" ] - block_L
      x2 = Grid [ i , "x" ] + block_L
      y1 = Grid [ i , "y" ] - block_L
      y2 = Grid [ i , "y" ] + block_L
      #subset to get only the sites within a disc of the reference point
      sites_in_square = sites[sites$x>x1 & sites$x<=x2 & sites$y>y1 & sites$y<=y2,c("x","y","ind")]
      sites_in_square$site_dist = sqrt ( ( sites_in_square$x - ref_point_x)^2  + (sites_in_square$y - ref_point_y)^2 )
      sites_in_block = sites_in_square[ sites_in_square$site_dist<block_L,]
      if (length (sites_in_block$ind) > 0){
        ind[[l]]=sites_in_block$ind
        l=l+1}
    }
  }
  if(shape == "square"){
    for ( i in 1:dim(Grid)[1]){
      x1 = Grid [ i , "x" ] - block_L/2
      x2 = Grid [ i , "x" ] + block_L/2
      y1 = Grid [ i , "y" ] - block_L/2
      y2 = Grid [ i , "y" ] + block_L/2
      sites_in_block = sites[ x>x1 & x<=x2 & y>y1 & y<=y2,]
      if (length (sites_in_block$ind) > 0){
        ind[[l]]=sites_in_block$ind
        l=l+1}
    }
  }
  return(ind)
}


