from numpy import pi

# ----------- GEO -------------------------------------------------------------
# mean earth radius in meter as defined by the International Union of
# Geodesy and Geophysics. Used for the spherical kd-tree and other things.
EARTH_RADIUS_M = 6371009.0
EARTH_RADIUS_KM = 6371009.0/1000.0
EARTH_CIRCUM_M = 2*pi*EARTH_RADIUS_KM
EARTH_CIRCUM_KM = 2*pi*EARTH_RADIUS_KM
DEG2M = EARTH_CIRCUM_M/360.0
DEG2KM = EARTH_CIRCUM_KM/360.0
M2DEG = 1.0/DEG2M
KM2DEG = 1.0/DEG2KM
