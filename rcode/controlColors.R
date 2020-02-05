# extract_colors_for_func_units-----
# Extract colors used for units
# get basic plot, that is already plotted...
funcUnitColors <- as_tibble(ggplot_build(plt)$data[[1]])
funcUnitColors <- funcUnitColors[ ,"colour"]
funcUnitColors$visMoveType <- tail(vertices$visMovType,-4)
funcUnitColors<- unique(funcUnitColors)
# funcUnitColors
# # A tibble: 4 x 2
# colour  visMoveType
# <chr>   <chr>      
# 1 #00BFC4 Vis        
# 2 #F8766D Mov        
# 3 #C77CFF VisMov     
# 4 #7CAE00 Other  

# other-----

df <- tibble(val=seq(1,256,by=1)   )
df$valAsChar <- sprintf("%03d",df$val)

p <- ggplot(df, aes(x = val, y = val, color = valAsChar, fill=valAsChar)) +
  geom_point(shape = 21, size = 4)

p
ZZ<-ggplot_build(p)

dat<-as_tibble(ZZ$data[[1]])
dat$val<-df$val
dat$valAsChar<-df$valAsChar



#p2 <- ggplot(dat[seq(1, 256, by = 1), ], aes(x = val, y = val, color = colour, fill = colour, stroke = 1)) +
  p2 <- ggplot(dat[seq(1, 256, by = 10), ], aes(x = 1, y = val, color = colour)) +
  geom_point(shape = 19, size = 2) +
    scale_fill_identity(guide = "none")
  scale_color_identity(
    "colorValue",
    labels = dat$colour,
    breaks = dat$colour,
    guide = "legend"
  )

p2
