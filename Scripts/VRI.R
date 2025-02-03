library(sf)
st_layers("/Users/Bronwyn/Documents/local-git/MPB/Input/VRI_JNP.gdb")

dist <- read_sf("/Users/Bronwyn/Documents/local-git/MPB/Input/VRI_JNP.gdb", layer ="JNP_VRI_DBO_Disturbance")
dist2 <- read_sf("/Users/Bronwyn/Documents/local-git/MPB/Input/VRI_JNP.gdb", layer ="JNP_VRI_DBO_Disturbance2")

