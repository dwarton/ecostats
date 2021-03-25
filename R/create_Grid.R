
create_Grid = function(Grid_space, x, y, shape,...){
  #Grid_space is the spacing of the Grid points which form the corners/ centres of the sampling units (squares/ circles)
  #x , y are coordinates of the sites
  
  #order x,y by x then y
  o=order(x,partial=y)
  
  x_range = max(x)-min(x)
  y_range = max(y)-min(y)
  nx = ceiling ( x_range/Grid_space ) #number of vertical grid lines
  ny = ceiling ( y_range/Grid_space ) #number of horizontal grid lines
  trimx =  ceiling ( x_range/Grid_space ) - x_range/Grid_space 
  trimy =  ceiling ( y_range/Grid_space ) - y_range/Grid_space 
  x0 = min ( x ) - trimx/2 * Grid_space  #where the vertical grid lines start
  y0 = min ( y ) - trimy/2 * Grid_space  #where the horizontal grid lines start
  
  xM = x0 + nx*Grid_space
  yM = y0 + ny*Grid_space
  
  
  #the x and y Grid coordinates
  if (shape=="disc"){
    fineGrid_x = x0 + 1:(nx-1)  * Grid_space
    fineGrid_y = y0 + 1:(ny-1)  * Grid_space
  }
  if(shape == "square"){
    seqx = c ( 0, 1:(nx-1))
    seqy = c ( 0, 1:(ny-1))
    fineGrid_x = x0 + seqx  * Grid_space
    fineGrid_y = y0 + seqy  * Grid_space
  }
  
  
  #the x and y Grid coordinates to resample
  
  Grid = merge ( fineGrid_x, fineGrid_y, all.x = TRUE, all.y = TRUE )
  return(Grid)
}