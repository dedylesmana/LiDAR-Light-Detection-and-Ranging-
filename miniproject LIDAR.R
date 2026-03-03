
############################################################
# MINI PROJECT LIDAR
# Plot 06 
############################################################

# ===============================
# 0. INSTALL & LOAD PACKAGE
# ===============================
install.packages(c("lidR","terra","sf","rgl","dplyr","ggplot2"))

library(lidR)
library(terra)
library(sf)
library(rgl)
library(dplyr)
library(ggplot2)
library(viridis)

###############################################
# 1. LOAD DATA & Visualisasi Las Data (TUGAS 1)
###############################################
las <- readLAS("plot_06.las")
las
summary(las)

if (is.empty(las)) stop("Data LAS tidak valid atau tidak berisi point cloud. 
                        Silakan periksa kembali file dan path")
las_kosong <- las[0]
is.empty(las_kosong)

# Visualisasi LAS Data
plot(las,
     color = "Z",
     size = 3,
     axis = TRUE,
     xlab = "X (meter)",
     ylab = "Y (meter)",
     zlab = "Tinggi (meter)")

##################################################################################
# TUGAS 2 - DETEKSI POHON menggunakan Canopy Height Model (PIT-FREE & Point Cloud)
##################################################################################
# NORMALISASI HEIGHT
dtm <- rasterize_terrain(las, res = 1, algorithm = knnidw())
las_norm <- normalize_height(las, dtm)
las_norm <- filter_poi(las_norm, Z >= 0)

# CHM PIT-FREE + DETEKSI POHON DARI CHM
chm <- rasterize_canopy(las, res = 0.5,
                        algorithm = pitfree(subcircle = 0.2))
chm_smooth <- focal(chm, w = 3, fun = mean, na.rm = TRUE)
ttops_chm <- locate_trees(chm_smooth, lmf(ws = 5))

# DETEKSI LANGSUNG DARI POINT CLOUD (LMF)
ttops_pc <- locate_trees(las, lmf(ws = 5))

# ================================
# VISUALISASI
# ================================
# CHM dari pitfree
chm_pf <- rasterize_canopy(
  las_norm,
  res = 0.5,
  algorithm = pitfree(subcircle = 0.15)
)

# Deteksi puncak pohon dari CHM (LMF)
ttops <- locate_trees(chm, lmf(ws = 5))

# DETEKSI POHON DARI POINT CLOUD 
# Deteksi puncak pohon langsung dari point cloud ternormalisasi
ttops_pc <- locate_trees(las_norm, lmf(ws = 5, h = 10))

# VISUAL COMPARISON (2 PANEL)
par(mfrow = c(1,2),
    mar = c(4,4,4,6),   
    oma = c(0,0,2,0))  

# ---- Panel A: CHM Pit-free + treetops dari CHM
plot(chm_pf,
     col = height.colors(60),
     main = "A) Deteksi Pohon pada CHM (Pit-free)",
     axes = TRUE)

plot(st_geometry(ttops_chm),
     add = TRUE,
     col = "red",
     pch = 16,
     cex = 0.8)

mtext(paste0("Jumlah treetop: ", nrow(ttops_chm)), side = 1, line = 3, cex = 0.9)

# ---- Panel B: Point cloud (Top-view) + treetops dari point cloud
plot(chm_pf,
     col = height.colors(60),
     main = "B) Deteksi Pohon dari Point Cloud (LMF)",
     axes = TRUE)

plot(st_geometry(ttops_pc),
     add = TRUE,
     col = "red",
     pch = 16,
     cex = 0.8)

mtext(paste0("Jumlah treetop: ", nrow(ttops_pc)), side = 1, line = 3, cex = 0.9)

mtext("Perbandingan Deteksi Pohon: CHM Pit-free vs Point Cloud",
      outer = TRUE, cex = 1.2, font = 2)

############################################################
# TUGAS 3 - SEGMENTASI POHON
############################################################
# a) SEGMENTASI DALPONTE 2016
algo_dalponte <- dalponte2016(chm, ttops)
las_dalponte <- segment_trees(las_norm, algo_dalponte)
n_dalponte <- length(unique(las_dalponte$treeID))
cat("Jumlah pohon Dalponte:", n_dalponte, "\n") #jumlah pohon hasil segmentasi

# b) SEGMENTASI LI (2012)
algo_li <- li2012()
las_li  <- segment_trees(las_norm, algo_li)
n_li <- length(unique(las_li$treeID))
cat("Jumlah pohon Li:", n_li, "\n") #jumlah pohon hasil segmentasi

df_dalponte <- las_dalponte@data
df_li       <- las_li@data

df_dalponte <- df_dalponte %>% filter(!is.na(treeID) & treeID > 0)
df_li       <- df_li %>% filter(!is.na(treeID) & treeID > 0)

hitung_area <- function(df){
  df %>%
    group_by(treeID) %>%
    summarise(
      area_m2 = {
        pts <- cbind(X, Y)
        if(nrow(pts) < 3) return(NA_real_)
        ch  <- chull(pts)
        poly <- pts[c(ch, ch[1]), ]
        abs(sum(poly[-1,1]*poly[-nrow(poly),2] -
                  poly[-nrow(poly),1]*poly[-1,2]))/2
      }
    )
}

crown_dalponte <- hitung_area(df_dalponte)
crown_li       <- hitung_area(df_li)

# EKSTRAKSI METRICS PER POHON
metrics_dalponte <- df_dalponte %>%
  group_by(treeID) %>%
  summarise(
    z_max  = max(Z),
    z_mean = mean(Z),
    z_sd   = sd(Z),
    n_pts  = n()
  )

metrics_li <- df_li %>%
  group_by(treeID) %>%
  summarise(
    z_max  = max(Z),
    z_mean = mean(Z),
    z_sd   = sd(Z),
    n_pts  = n()
  )

# GABUNGKAN AREA + METRICS
hasil_dalponte <- left_join(metrics_dalponte, crown_dalponte, by="treeID")
hasil_li       <- left_join(metrics_li, crown_li, by="treeID")

cat("====================================\n")
cat("Ringkasan Segmentasi\n")
cat("Dalponte :", nrow(hasil_dalponte), "pohon\n")
cat("Li       :", nrow(hasil_li), "pohon\n")
cat("====================================\n")

head(hasil_dalponte)
head(hasil_li)

summary(hasil_dalponte$area_m2)
summary(hasil_li$area_m2)

############################################################
# TUGAS 4 - EKSTRAKSI METRICS
############################################################
# DALPONTE
metrics_dalponte_table <- crown_metrics(
  las_dalponte,
  func = ~list(
    z_max = max(Z, na.rm=TRUE),
    z_std = sd(Z, na.rm=TRUE),
    z_med = median(Z, na.rm=TRUE),
    i_mean = mean(Intensity, na.rm=TRUE),
    i_max = max(Intensity, na.rm=TRUE)
  )
)

metrics_dalponte_table <- st_drop_geometry(metrics_dalponte_table)

# LI
metrics_li_table <- crown_metrics(
  las_li,
  func = ~list(
    z_max = max(Z, na.rm=TRUE),
    z_std = sd(Z, na.rm=TRUE),
    z_med = median(Z, na.rm=TRUE),
    i_mean = mean(Intensity, na.rm=TRUE),
    i_max = max(Intensity, na.rm=TRUE)
  )
)

metrics_li_table <- st_drop_geometry(metrics_li_table)

# Lihat hasil
head(metrics_dalponte_table)
head(metrics_li_table)

write.csv(metrics_dalponte_table, "metrics_dalponte.csv", row.names=FALSE)
write.csv(metrics_li_table, "metrics_li.csv", row.names=FALSE)

############################################################
# TUGAS 5 - HITUNG AGB (Black 2004)
############################################################
# Rumus umum: AGB = a × (Height^b)
# Menggunakan z_max sebagai tinggi pohon

a <- 0.0673
b <- 2.5

metrics_dalponte <- metrics_dalponte %>%
  mutate(AGB = a * (z_max^b))

metrics_li <- metrics_li %>%
  mutate(AGB = a * (z_max^b))

summary(metrics_dalponte$AGB)
summary(metrics_li$AGB)

# Total biomassa (kg)
total_agb_dalponte <- sum(metrics_dalponte$AGB, na.rm = TRUE)
total_agb_li <- sum(metrics_li$AGB, na.rm = TRUE)

# Rata-rata biomassa per pohon (kg)
mean_agb_dalponte <- mean(metrics_dalponte$AGB, na.rm = TRUE)
mean_agb_li <- mean(metrics_li$AGB, na.rm = TRUE)

# Tampilkan hasil
total_agb_dalponte
total_agb_li

mean_agb_dalponte
mean_agb_li

############################################################
# TUGAS 6 - VISUALISASI AGB
############################################################
# Visualisasi AGB Dalponte 2016
library(patchwork)
make_crowns <- function(las_obj){
  df <- las_obj@data %>%
    filter(!is.na(treeID) & treeID > 0)
  
  pts <- st_as_sf(df, coords = c("X","Y"), remove = FALSE)
  
  crowns <- pts %>%
    group_by(treeID) %>%
    summarise(geometry = st_convex_hull(st_union(geometry)),
              .groups = "drop")
    return(crowns)
}

crowns_dalponte <- make_crowns(las_dalponte)
crowns_li       <- make_crowns(las_li)

crowns_dalponte <- crowns_dalponte %>%
  left_join(metrics_dalponte %>% select(treeID, AGB), by="treeID") %>%
  mutate(AGB_ton = AGB/1000)

crowns_li <- crowns_li %>%
  left_join(metrics_li %>% select(treeID, AGB), by="treeID") %>%
  mutate(AGB_ton = AGB/1000)

global_min <- min(c(crowns_dalponte$AGB_ton,
                    crowns_li$AGB_ton), na.rm=TRUE)

global_max <- max(c(crowns_dalponte$AGB_ton,
                    crowns_li$AGB_ton), na.rm=TRUE)

p1 <- ggplot(crowns_dalponte) +
  geom_sf(aes(fill = AGB_ton), color = NA) +
  scale_fill_viridis_c(
    name = "AGB (Ton)",
    option = "C",
    limits = c(global_min, global_max)
  ) +
  labs(
    title = "Peta Distribusi Biomassa Diatas Tanah (AGB)",
    subtitle = "Estimasi berbasis LiDAR Individu Pohon (ITD)-Dalponte (2016)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face="bold", size=11, hjust=0.5),
    plot.subtitle = element_text(size=9, hjust=0.5),
    legend.position = "right"
  ) +
  annotate("text",
           x = Inf, y = -Inf,
           label = "Lokasi: Perm Krai, Rusia",
           hjust = 1.1, vjust = -1,
           size = 3)
p2 <- ggplot(crowns_li) +
  geom_sf(aes(fill = AGB_ton), color = NA) +
  scale_fill_viridis_c(
    name = "AGB (Ton)",
    option = "C",
    limits = c(global_min, global_max)
  ) +
  labs(
    title = "Peta Distribusi Biomassa Diatas Tanah (AGB)",
    subtitle = "Estimasi berbasis LiDAR Individu Pohon (ITD)-Li (2012)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face="bold", size=11, hjust=0.5),
    plot.subtitle = element_text(size=9, hjust=0.5),
    legend.position = "right"
  ) +
  annotate("text",
           x = Inf, y = -Inf,
           label = "Lokasi: Perm Krai, Rusia",
           hjust = 1.1, vjust = -1,
           size = 3)
final_plot <- (p1 | p2) +
  plot_annotation(
    caption = "Gambar 6. Visualisasi hasil AGB menggunakan allometrik Black (2004)"
  ) &
  theme(
    plot.caption = element_text(face="italic",
                                size=12,
                                hjust=0.5)
  )

print(final_plot)

ggsave("AGB_Visualisasi_Dalponte_vs_Li.png",
       final_plot,
       width = 12,
       height = 6,
       dpi = 300)


combined_plot <- (p1 | p2) +
  plot_annotation(
    caption = "Gambar 6. Visualisasi hasil AGB menggunakan allometrik Black (2004)"
  )

print(combined_plot)

