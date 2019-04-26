## Spatial models: R club 26/04 2019 Kristian H. Haugen
# Packages ----------------------------------------------------------------
library(pacman)
p_load(data.table, sp, spdep, ggplot2)

# Loading data ----------------------------------------------------------------
data <- fread('https://raw.githubusercontent.com/KHHaugen/RURED/master/grunnstruktur%20med%20geokoder.csv')


# Prepping data ----------------------------------------------------------------
# Removing schools located outside the country
data$lon[data$lat < 57.5] <- NA
data$lat[data$lat < 57.5] <- NA

data$lat[(data$lat < 66 | data$lat > 72) & data$lon > 15] <- NA
data$lon[is.na(data$lat)] <- NA

data <- data[!is.na(data$lat), ]

# Removing Språksenteret and schools with reported 0 students
data <- data[data$Orgnr != 813993422,]
data <- data[data$grpst2 != 0,]

## The data is panel data, and we need it to be cross sectional. We therefore aggregate the data.
data <- aggregate(data[,c('grpoeng', 'grpst2', 'utdanning_komb', 'folkemengde', 'eleverN', 'lon', 'lat')], list(data$Orgnr, data$Navn), mean, na.rm=T)

data <- data[!is.na(data$utdanning_komb),]
data <- data[!is.na(data$grpst2),]
data <- data[!is.na(data$grpoeng),]
data$folkemengde[data$folkemengde == 0] <- NA
data <- data[!is.na(data$folkemengde),]

# Letting R know it's spatial data and what system is used ----------------------------------------------------------------
# Defining coordinates and CRS
sp_point <- cbind(data$lon, data$lat)
proj <- CRS("+init=epsg:4326")


# We are working with points, and a problem arise when two points are exactly the same. 
# This is an issue here, and one way to handle this is to set an offset to the location of duplicates.

sp_point[,2][duplicated(sp_point)] <- sp_point[,2][duplicated(sp_point)] + 0.00001

data.sp <- SpatialPointsDataFrame(coords=sp_point,data,proj4string=proj)
data.sp <- spTransform(data.sp, CRS("+init=epsg:25833"))
map_crd <- coordinates(data.sp)
rownames(map_crd) <- data.sp$id

# Defining neighbors ----------------------------------------------------------------
k <- 3
W_knn1 <- knn2nb(knearneigh(map_crd, k=k))
W_knn1_mat <- nb2listw(W_knn1, style = 'W', zero.policy = TRUE)

# Displaying connectivity ----------------------------------------------------------------
plot(W_knn1_mat,coords=map_crd,pch=19, cex=0.1, col="gray")

# Testing for global spatial correlation ----------------------------------------------------------------
moran.test(data$grpoeng, listw = W_knn1_mat)

# Showing distribution of GPA
data.sp$percentile<-(rank(data.sp$grpoeng)/length(data.sp$grpoeng))*100
#mid <- mean(data.sp$grpoeng)

GPA <- ggplot(as.data.frame(data.sp), aes(map_crd[,1], map_crd[,2])) +
  geom_point(aes(colour = percentile),
             size = 0.5) +
  scale_color_gradient(low="red", high="blue") +
  labs(title = 'GPA of schools', x = 'lon', y = 'lat') +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

GPA

# Models ----------------------------------------------------------------
mod1 <- lm(grpoeng ~ grpst2 + I(log(folkemengde)) + utdanning_komb, data = data) # OLS
mod2 <- errorsarlm(grpoeng ~ grpst2 + I(log(folkemengde)) + utdanning_komb, data = data, listw = W_knn1_mat) # simultaneous autoregressive model
mod3 <- lagsarlm(grpoeng ~ grpst2 + I(log(folkemengde)) + utdanning_komb, data = data, listw = W_knn1_mat) # Spatially lagged dependent variable

summary(mod1)
summary(mod2)
summary(mod3)

# Geostatistical models ----------------------------------------------------------------
# Geostatistical models are spatial models where the connectivity is based on distance, and not on neighbors necessarily.
dist <- unlist(nbdists(W_knn1, map_crd)) # Unlisting to use the information
distance <- dnearneigh(map_crd, d1=0, d2=max(dist)) # Calculating distances. d1 and d2 creates a band of distances to be used. Here the band is set so that all points have a weight, but it can be set so that for example only points in a set distance will be used.
distance.neigh.list <- nbdists(distance, coordinates(map_crd))
W_inv.dist <- lapply(distance.neigh.list, function(x) 1/x) # Weighting on the inverse of the distance. The points nearer to will have a larger weight than points further away

W_inv.distance <- nb2listw(distance, glist=W_inv.dist, style="W") 

mod4 <-  errorsarlm(grpoeng ~ grpst2 + I(log(folkemengde)) + utdanning_komb, data = data, listw=W_inv.distance, zero.policy=T) # Geostatistical error correcting model
mod5 <- lagsarlm(grpoeng ~ grpst2 + I(log(folkemengde)) + utdanning_komb, data = data.sp, listw=W_inv.distance, zero.policy=T) # Geostatistical spatially lagged dependent

summary(mod4)
summary(mod5)