#!/bin/bash

#ARGS (1) sphere file (herbold format) (2) unit sphere mesh (3) output ABAQUS file
if [ $# -ne 3 ]; then
   echo "usage: ARGS (1) sphere file (herbold format) (2) unit sphere mesh (3) output ABAQUS file"
   exit 1
fi 

#ASSIGN DEFAULT ENVIRONMENTAL VARIABLES
if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/gpac_trunk/scripts
    if [ -d "$GPAC_SCRIPT_PATH" ];then
	echo "$0 using GPAC_SCRIPTH_PATH $GPAC_SCRIPT_PATH"
    else
	export GPAC_SCRIPT_PATH=/usr/gapps/GEOS/scripts
	if [ -d "$GPAC_SCRIPT_PATH" ];then
	    echo "$0 using GPAC_SCRIPTH_PATH $GPAC_SCRIPT_PATH"
	else
	    echo "$0 cannot find GPAC_SCRIPT_PATH $GPAC_SCRIPT_PATH"
	    exit 1
	fi
    fi
fi

#ASSIGN DEFAULT ENVIRONMENTAL VARIABLES
if [ -z "$GPAC_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_EXEC=~/apps/gpac_trunk/src/GPAC.x
    if [ -e "$GPAC_EXEC" ];then
	echo "$0 using GPAC_EXEC $GPAC_EXEC"
    else
	echo "$0 cannot find GPAC_EXEC $GPAC_EXEC"
	exit 1
    fi
fi

if [ -z "$VISIT_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export VISIT_EXEC=/Applications/VisIt.app/Contents/MacOS/VisIt
    if [ -e "$VISIT_EXEC" ];then
	echo "$0 using VISIT_EXEC $VISIT_EXEC"
    else
	export VISIT_EXEC=/usr/gapps/visit/bin/visit
	if [ -e "$VISIT_EXEC" ];then
	    echo "$0 using VISIT_EXEC $VISIT_EXEC"
	else
	    echo "$0 cannot find VISIT_EXEC $VISIT_EXEC"
	    exit 1
	fi
    fi
fi

#REQUIRED SCRIPTS
#(1) gpac_shear_permeability_input.pl
#(2) gpac_surface.sh (see for more required scripts)
#(3) 

# -v 1.5.5 -nowin -cli -l srun -p pbatch -np 80 -s script.py -verbose -geometry 512x512+0+0 -switch ib -x prism34

#set problem dimensions
X=5
Y=5
Z=1

#write initial GPAC and APGEN input scripts
perl $GPAC_SCRIPT_PATH/gpac_shear_permeability_input.pl 0 > meshx.xml
perl $GPAC_SCRIPT_PATH/gpac_shear_permeability_input.pl 1 > meshy.xml
echo "0.03 0.02" > statparams

#run Monte Carlo
for i in {0..0}
do
    #1) GENERATE SURFACE PROFILE
    bash $GPAC_SCRIPT_PATH/gpac_surface.sh statparams 5 5 1.2 $X $Y $Z `perl -e 'my $x = shift; $x *= 0.1;print $x;' -f $Z` 1
    mv mesh.inp mesh.geom

    #2) RUN GPAC (X-DIRECTION)
    cp meshx.xml mesh.xml

    $GPAC_EXEC -i mesh.xml -m mesh.geom
     
    #3) ANALYZE RESULTS (X-DIRECTION)
    perl $GPAC_SCRIPT_PATH/gpac_query_shear_and_tau.pl mesh.geom `pwd` tmp.session 0 2 0 > tmp.py
    #top displacement (averaged)
    $VISIT_EXEC -nowin -cli -s tmp.py > tmp0

    perl $GPAC_SCRIPT_PATH/gpac_query_shear_and_tau.pl mesh.geom `pwd` tmp.session 0 2 1 > tmp.py
    #bottom displacement (averaged)
    $VISIT_EXEC -nowin -cli -s tmp.py > tmp1

    paste -d ' ' tmp0 tmp1 > shear_tau
    mv shear_tau tmp0

    perl $GPAC_SCRIPT_PATH/gpac_query_shear_and_tau.pl mesh.geom `pwd` tmp.session 0 2 2 > tmp.py
    #common plane force (sum of tau*area)
    $VISIT_EXEC -nowin -cli -s tmp.py > tmp1

    paste -d ' ' tmp0 tmp1 > shear_tau
    mv shear_tau tmp0

    perl $GPAC_SCRIPT_PATH/gpac_query_permeability.pl mesh.geom `pwd` tmp.session 0 2 > tmp.py
    #get the flow in the x-direction
    $VISIT_EXEC -nowin -cli -s tmp.py > tmp1

    paste -d ' ' tmp0 tmp1 > shear_tau
    mv shear_tau tmp0

    #4) RUN GPAC (Y-DIRECTION)
    cp meshy.xml mesh.xml
    $GPAC_EXEC -i mesh.xml -m mesh.geom
    
    #5) ANALYZE RESULTS (Y-DIRECTION)
    perl $GPAC_SCRIPT_PATH/gpac_query_permeability.pl mesh.geom `pwd` tmp.session 1 2 > tmp.py
    $VISIT_EXEC -nowin -cli -s tmp.py > tmp1

    paste -d ' ' tmp0 tmp1 > shear_tau
    mv shear_tau dataset$i
done

exit
