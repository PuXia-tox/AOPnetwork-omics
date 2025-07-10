annotation_compass <- function(label,
                               position = c('N','NE','E','SE','S','SW','W','NW'),
                               padding = grid::unit(c(0.5,0.5),"line"), ...){
  position <- match.arg(position)
  x <- switch (position,
               N = 0.5,
               NE = 1,
               E = 1,
               SE = 1,
               S = 0.5, 
               SW = 0,
               W = 0, 
               NW = 0
  )
  y <- switch (position,
               N = 1,
               NE = 1,
               E = 0.5,
               SE = 0,
               S = 0, 
               SW = 0,
               W = 0.5, 
               NW = 1
  )
  hjust <- switch (position,
                   N = 0.5,
                   NE = 1,
                   E = 1,
                   SE = 1,
                   S = 0.5, 
                   SW = 0,
                   W = 0, 
                   NW = 0
  )
  vjust <- switch (position,
                   N = 1,
                   NE = 1,
                   E = 0.5,
                   SE = 0,
                   S = 0, 
                   SW = 0,
                   W = 0.5, 
                   NW = 1
  )
  f1 <- switch (position,
                N = 0,
                NE = -1,
                E = -1,
                SE = -1,
                S = 0, 
                SW = 1,
                W = 1, 
                NW = 1
  )
  f2 <- switch (position,
                N = -1,
                NE = -1,
                E = 0,
                SE = 1,
                S = 1, 
                SW = 1,
                W = 0, 
                NW = -1
  )
  annotation_custom(grid::textGrob(label,
                                   x=grid::unit(x,"npc") + f1*padding[1] , 
                                   y=grid::unit(y,"npc") + f2*padding[2],
                                   hjust=hjust,vjust=vjust, ...))
}