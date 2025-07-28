library(rcarbon)
library(dplyr)
library(sf)
library(spatstat)


Data<- read.csv(here("Data", "/After_Combine.csv"))

# Excluding outliers
Data_clean <- Data %>% filter(is.na(Mismatch_Flag)) %>% filter(!is.na(Longitude) & !is.na(Latitude)) 

# Loading map
# http://www.gisdeveloper.co.kr/download/admin_shp/ctprvn_20230729.zip
map_korea <- st_read(here("ctprvn_20230729", "CTPRVN.shp"))

map_korea$CTPRVN_CD <- iconv(map_korea$CTPRVN_CD,
                             from = 'CP949',
                             to = 'UTF-8',
                             sub = NA)

# 2. Filtering for Han River basin
target_region_codes <- c('11', '41', '51')
subset_map <- map_korea %>% filter(CTPRVN_CD %in% target_region_codes)

# 3. Sf to spatial
subset_map_sp <- as(subset_map, "Spatial")

class(subset_map_sp)

Han_River_Basin<- as.owin(subset_map_sp)
plot(Han_River_Basin)

#Converting
convert_to_EPSG5179 <- function(data, lon_col, lat_col) {
  crs_wgs84 <- st_crs(4326)
  crs_korea2000 <- st_crs(5179)
  
  input_sf <- st_as_sf(data.frame(longitude = data[[lon_col]], latitude = data[[lat_col]]),
                       coords = c("longitude", "latitude"), crs = crs_wgs84)
  
  output_sf <- st_transform(input_sf, crs = crs_korea2000)
  output_df <- st_coordinates(output_sf)
  
  cbind(data, output_df)
}


dir.create(file.path("im"), showWarnings = FALSE)
mapData <- convert_to_EPSG5179(Data_clean, lon_col = "Longitude", lat_col = "Latitude")

mapData.cal <- calibrate(x=mapData$BP, errors=mapData$error, normalised=TRUE)

stkde <- stkde(x=mapData.cal, coords=mapData[,c("X","Y")],
                win=Han_River_Basin, sbw=40000, cellres=5000, focalyears=seq(3500, 1100, -200),
                tbw=50, bins=NA, backsight=200, outdir="im")

png("HanRiverBasinfocalchange3500-1100.png",width=9000,height=6000,res=500)
par(mfrow=c(3,4), mar=c(2,2,3,1))
plot(stkde, 3500, type="change",cex.main=3, cex.axis=1.8,font=2, main="(a) 3500-3300 cal BP")
plot(stkde, 3300, type="change",cex.main=3, cex.axis=1.8,font=2, main="(b) 3300-3100 cal BP")
plot(stkde, 3100, type="change",cex.main=3, cex.axis=1.8,font=2,main="(c) 3100-2900 cal BP")
plot(stkde, 2900, type="change",cex.main=3, cex.axis=1.8,font=2, main="(d) 2900-2700 cal BP")
plot(stkde, 2700, type="change",cex.main=3, cex.axis=1.8,font=2, main="(e) 2700-2500 cal BP")
plot(stkde, 2500, type="change",cex.main=3, cex.axis=1.8,font=2, main="(f) 2500-2300 cal BP")
plot(stkde, 2300, type="change",cex.main=3, cex.axis=1.8,font=2, main="(g) 2300-2100 cal BP")
plot(stkde, 2100, type="change",cex.main=3, cex.axis=1.8,font=2, main="(h) 2100-1900 cal BP")
plot(stkde, 1900, type="change",cex.main=3, cex.axis=1.8,font=2, main="(i) 1900-1700 cal BP")
plot(stkde, 1700, type="change",cex.main=3, cex.axis=1.8,font=2, main="(j) 1700-1500 cal BP")
plot(stkde, 1500, type="change",cex.main=3, cex.axis=1.8,font=2, main="(k) 1500-1300 cal BP")
plot(stkde, 1300, type="change",cex.main=3, cex.axis=1.8,font=2, main="(l) 1300-1100 cal BP")
dev.off()


png("HanRiverBasinfocalintensity3500-1100.png",width=9000,height=6000,res=500)
par(mfrow=c(3,4), mar=c(2,2,3,1))
plot(stkde, 3500, type="focal",cex.main=3, cex.axis=1.8,main="(a) 3500-3300 cal BP")
plot(stkde, 3300, type="focal",cex.main=3, cex.axis=1.8,main="(b) 3300-3100 cal BP")
plot(stkde, 3100, type="focal",cex.main=3, cex.axis=1.8,main="(c) 3100-2900 cal BP")
plot(stkde, 2900, type="focal",cex.main=3, cex.axis=1.8, main="(d) 2900-2700 cal BP")
plot(stkde, 2700, type="focal",cex.main=3, cex.axis=1.8, main="(e) 2700-2500 calBP")
plot(stkde, 2500, type="focal",cex.main=3, cex.axis=1.8, main="(f) 2500-2300 cal BP")
plot(stkde, 2300, type="focal",cex.main=3, cex.axis=1.8, main="(g) 2300-2100 cal BP")
plot(stkde, 2100, type="focal",cex.main=3, cex.axis=1.8, main="(h) 2100-1900 cal BP")
plot(stkde, 1900, type="focal",cex.main=3, cex.axis=1.8, main="(i) 1900-1700 cal BP")
plot(stkde, 1700, type="focal",cex.main=3, cex.axis=1.8, main="(j) 1700-1500 cal BP")
plot(stkde, 1500, type="focal",cex.main=3, cex.axis=1.8, main="(k) 1500-1300 cal BP")
plot(stkde, 1300, type="focal",cex.main=3, cex.axis=1.8,main="(l) 1300-1100 cal BP")
dev.off()

