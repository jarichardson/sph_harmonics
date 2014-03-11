#example mapping script
#minmax gggrx_dg.llg
WEST=-60
EAST=-25
NORTH=55
SOUTH=25
INSPACING=0.25
GRDSPACING=0.1

LON0=$[$WEST+($EAST-$WEST)/2]
LAT0=$[$SOUTH+($NORTH-$SOUTH)/2]
LAT1=$SOUTH
LAT2=$NORTH

LATGRATICULE=5
LONGRATICULE=5
CONTOUR=50

IN="marius/gggrx_dg.llg"
GRD="shall_v_deep/boug_unweighted.grd"
HILL="shall_v_deep/boug_unweighted.hill.grd"

#shallow.251-320.grd
#full.2-320.grd
#deep.2-250.grd

XY="domes.xy"
CPT="shall_v_deep/grav_anom_moon.cpt"
EPS="felsic.boug_unweighted.eps"

#gmtinfo $IN
#xyz2grd $IN -R$WEST/$EAST/$SOUTH/$NORTH -I$INSPACING -G$GRD
#surface $IN -R$WEST/$EAST/$SOUTH/$NORTH  -I$GRDSPACING  -G$GRD -T0.25 -C0.1 -V
gmt makecpt -Chaxby -D -T-400/400/100 >$CPT
gmt grdgradient $GRD -A0/270 -Ne0.2 -G$HILL

gmt grdimage $GRD -Y3i --PROJ_ELLIPSOID=Moon -R$WEST/$EAST/$SOUTH/$NORTH -C$CPT -I$HILL -JB$LON0/$LAT0/$LAT1/$LAT2/6i -Ba$LONGRATICULEg$LONGRATICULE/a$LATGRATICULEg$LATGRATICULE -K -P > $EPS
gmt grdcontour $GRD --PROJ_ELLIPSOID=Moon -J -R -C$CONTOUR -Wcthick -K -O >> $EPS
gmt psxy $XY -R -J -O -K --PROJ_ELLIPSOID=Moon -Sx0.3c -Wthick -V >> $EPS

gmt psscale  -D3i/-0.5i/9c/0.65ch --FONT_ANNOT_PRIMARY=10p --FONT_LABEL=12p --MAP_LABEL_OFFSET=0.3c -O -C$CPT -I -E -B100:"Bouguer Gravity Anomaly (2560 kg/m3)":/:mGal: >> $EPS

gmt ps2raster $EPS -A -Tg -P 
rm $EPS

