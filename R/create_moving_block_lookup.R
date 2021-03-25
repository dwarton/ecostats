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
      sites_in_block = sites %>%
        filter (x>x1, x<x2, y>y1, y<y2)%>%
        mutate (site_dist = sqrt ( ( x - ref_point_x)^2  + (y - ref_point_y)^2 ))%>%
        subset (site_dist<block_L)%>%
        dplyr::select (ind)
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
      sites_in_block = sites %>%
        filter (x>x1, x<x2, y>y1, y<y2)%>%
        dplyr::select (ind)
      if (length (sites_in_block$ind) > 0){
        ind[[l]]=sites_in_block$ind
        l=l+1}
    }
  }
  return(ind)
}



