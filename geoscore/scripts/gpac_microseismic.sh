#!/bin/bash

if [ $# -ne 1 ]; then
   echo "usage: ARGS (1) GEOS standard output file"
   exit 1
fi 
if [ -e $1 ]; then
    echo "Found $1 GEOS standard output file"
else
    echo "Cannot find $1 GEOS standard output file"
    exit 1
fi

#ASSIGN DEFAULT ENVIRONMENTAL VARIABLES
if [ -z "$GEOS_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GEOS_SCRIPT_PATH=~/apps/geos_git/gpaca/scripts
    if [ -d "$GEOS_SCRIPT_PATH" ];then
	echo "$0 using GEOS_SCRIPTH_PATH $GEOS_SCRIPT_PATH"
    else
	export GEOS_SCRIPT_PATH=/usr/gapps/GEOS/scripts
	if [ -d "$GEOS_SCRIPT_PATH" ];then
	        echo "$0 cannot find GEOS_SCRIPT_PATH $GEOS_SCRIPT_PATH"
		    exit 1
		    fi
    fi
fi

GSO=$1
if [ -d dfiles ]; then
    echo "Found dfiles subdirectory"
else
    echo "Creating dfiles subdirectory"
    mkdir dfiles
fi

grep ^Tape $GSO | awk '{if($4 > 0){print $2" "$3}}' > dfiles/beachpts_betadist_points.dat
grep ^AseismicTape $GSO | awk '{if($4 > 0){print $2" "$3}}' > dfiles/beachpts_betaaseis_points.dat
grep irac $GSO > sources
perl -e 'my $fname = shift; open(IN,"<$fname"); while(<IN>){if(/^Aseis/){print "0\n";}elsif(/^Tape/){print "1\n";}} close IN;' -f $GSO > aseis

#get points and threshold points
perl $GEOS_SCRIPT_PATH/wpp2points.pl sources > pts

#gutenberg-richter
cat pts | awk '{print ($5/(3600*24*365))" "(log(1000000*$4)/log(10) - 9.1)/1.5}' > mags
perl $GEOS_SCRIPT_PATH/gutenberg_richter.pl mags 9.83e-5 -4 -3.5 -3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 > gr

#get XMDV file
paste -d ' ' pts aseis > pts_aseis
cat pts_aseis | awk '{if($4 > 0){print $1" "$2" "$3" "(log(1000000*$4)/log(10) - 9.1)/1.5" "$5" "$6}}' > pts1k
cat pts1k | awk '{if($6 > 0){print $1" "$2" "$3" "$4" "$5" "$6}}' > pts1k_seis
perl $GEOS_SCRIPT_PATH/points2xmdv.pl pts1k x y z M0 time seismic > pts1k.okc

#get plots
if [ -e "lune.csh" ]; then
    echo "Using pre-existing lune.csh file"
else
    echo "Creating lune.csh"
    LFN=lune.csh
    echo "/opt/local/lib/gmt4/bin/gmtset PAPER_MEDIA letter MEASURE_UNIT inch BASEMAP_TYPE plain PLOT_DEGREE_FORMAT D TICK_LENGTH 0.3c LABEL_FONT_SIZE 10 ANOT_FONT_SIZE 10  HEADER_FONT 1 ANOT_FONT 1 LABEL_FONT 1 HEADER_FONT_SIZE 18 FRAME_PEN 2p TICK_PEN 2p" > $LFN
    echo "/opt/local/lib/gmt4/bin/psbasemap -JH0/3i -R-30/30/-90/90 -Ba10f5g10:\" \":/a10f5g10:\" \":wesn -G200 -K -V -P -X2 -Y1 > lune_hammer_iplot1.ps" >> $LFN
    echo "/opt/local/lib/gmt4/bin/psxy ./dfiles//beach_patch_01.lonlat -G120 -J -R -K -O -V >>lune_hammer_iplot1.ps" >> $LFN
    echo "/opt/local/lib/gmt4/bin/psxy ./dfiles//beach_patch_02.lonlat -G255 -J -R -K -O -V >>lune_hammer_iplot1.ps" >> $LFN
    echo "/opt/local/lib/gmt4/bin/psbasemap -JH0/3i -R-30/30/-90/90 -Ba10f5g10:\" \":/a10f5g10:\" \":wesn -K -O -V >> lune_hammer_iplot1.ps" >> $LFN
    echo "/opt/local/lib/gmt4/bin/psxy ./dfiles//beach_arc_01.lonlat -W2p,0,-- -J -R -K -O -V >>lune_hammer_iplot1.ps" >> $LFN
    echo "/opt/local/lib/gmt4/bin/psxy ./dfiles//beach_arc_02.lonlat -W2p,0,-- -J -R -K -O -V >>lune_hammer_iplot1.ps" >> $LFN
    echo "/opt/local/lib/gmt4/bin/psxy ./dfiles//beach_arc_03.lonlat -W2p,0,-- -J -R -K -O -V >>lune_hammer_iplot1.ps" >> $LFN
    echo "/opt/local/lib/gmt4/bin/psxy ./dfiles//beach_arc_04.lonlat -W2p,0,-- -J -R -K -O -V >>lune_hammer_iplot1.ps" >> $LFN
    echo "/opt/local/lib/gmt4/bin/psxy ./dfiles//beach_arc_05.lonlat -W2p,0,-- -J -R -K -O -V >>lune_hammer_iplot1.ps" >> $LFN
    echo "/opt/local/lib/gmt4/bin/psxy ./dfiles//beach_arc_06.lonlat -W2p,0,-- -J -R -K -O -V >>lune_hammer_iplot1.ps" >> $LFN
    echo "/opt/local/lib/gmt4/bin/psxy ./dfiles//beachpts_betadist_points.dat -N -Sc4p -W0.5p,0/0/0 -G0/255/255 -J -R -K -O -V >>lune_hammer_iplot1.ps" >> $LFN
    echo "/opt/local/lib/gmt4/bin/psxy -N -Sc4p -W0.5p,0/0/0 -G0/255/255 -R0/1/0/1 -JX1i -Xa3.5 -Ya6 -K -O -V >>lune_hammer_iplot1.ps<<EOF" >> $LFN
    echo "0 1.2" >> $LFN
    echo "EOF" >> $LFN
    echo "/opt/local/lib/gmt4/bin/pstext -N -R0/1/0/1 -JX1i -Xa3.5 -Ya6 -K -O -V >>lune_hammer_iplot1.ps<<EOF" >> $LFN
    echo " 0.2 1.2 12 0 1 LM betadist (n=175)" >> $LFN
    echo "EOF" >> $LFN
    echo "/opt/local/lib/gmt4/bin/pstext -R0/1/0/1 -JX1i -Xa-1 -Ya8.7 -O -V >>lune_hammer_iplot1.ps<<EOF" >> $LFN
    echo " -1 -1 11 0 1 LM TEST" >> $LFN
    echo "EOF" >> $LFN
fi
#perl $GEOS_SCRIPT_PATH/gammadelta2lune.pl betadist betaaseis
perl $GEOS_SCRIPT_PATH/gammadelta2lune.pl betadist
