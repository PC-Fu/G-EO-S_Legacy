#!/bin/bash

#ARGS (1) Minimum burn-in time in years (2) range of burn-in time perturbation in years (3) post burn-in simulation duration in years (4) number of realizations (5) directory of pore pressures (6) number of pore pressure files
if [ $# -ne 6 ]; then
   echo "Needed: ppvis.session in CWD and pp_reXX in directory of pore pressures - this directory should be where ppvis.session points for saving visualization files!"
   echo "usage $0 : (1) Minimum burn-in time in years (2) range of burn-in time perturbation in years (3) post burn-in simulation duration in years (4) number of realizations (5) directory of pore pressures (6) number of pore pressure files"
   exit 1
fi 

if [ -d $5 ]; then
    echo "Using pp_reXX in $5"
else
    echo "$5 -> directory does not exist"
    exit 1
fi

if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/geos_git/gpaccore/scripts
    if [ -d $GPAC_SCRIPT_PATH ]; then
        echo "$0 setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
    else
	export GPAC_SCRIPT_PATH=~/apps/geos_git/gpaca/scripts
	if [ -d $GPAC_SCRIPT_PATH ]; then
            echo "$0 setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
	else
            echo "$0 EXITING: cannot find GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
            exit 1
	fi
    fi
fi

if [ -z "$GPAC_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_EXEC=~/apps/geos_git/gpaccore/GPAC.opt
    if [ -e $GPAC_EXEC ]; then
        echo "$0 setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
    else
	export GPAC_EXEC=~/apps/geos_git/gpaca/GPAC.opt
	if [ -e $GPAC_EXEC ]; then
            echo "$0 setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
	else
            echo "$0 EXITING: cannot find GPAC_EXEC: $GPAC_EXEC"
            exit 1
	fi
    fi
fi

if [ -z "$VISIT" ]; then
    #NOTE: this only provides a value for the current script
    export VISIT=/Applications/VisIt.app/Contents/MacOS/VisIt
    if [ -e $VISIT ]; then
        echo "$0 setting default for VISIT: $VISIT"
    else
        echo "$0 EXITING: cannot find VISIT: $VISIT"
        exit 1
    fi
fi


#PRINT THE VISIT PYTHON SCRIPT
echo "----> WRITING PYTHON FILE <----"
echo "import sys" > tmp.py
echo "RestoreSession(\"ppvis.session\",0)" >> tmp.py
echo "nplots = TimeSliderGetNStates()" >> tmp.py
echo "start=0" >> tmp.py
echo "end = nplots" >> tmp.py
echo "for i in range(start, end):" >> tmp.py
echo "    SetTimeSliderState(i)" >> tmp.py
echo "    Query(\"Max\")" >> tmp.py
echo "    print \"%g\" % GetQueryOutputValue()" >> tmp.py
echo "    SaveWindow()" >> tmp.py
echo "sys.exit(\"\")" >> tmp.py

#GO THROUGH EACH PP FILE
for i in `seq 1 $6`
do
    if [ -e maxpps$i ]; then
	echo "----> REMOVING maxpps$i <----"
	rm maxpps$i
    fi
    echo "----> REALIZATION $i <----"

    #remove old files
    rm pp ms_* restart_*

    #create new pore pressure symbolic link
    ln -s $5/pp_re$i pp
    #$GPAC_EXEC -i run.xml
    $GPAC_SCRIPT_PATH/gpac_run_multiple_faults.sh $1 $2 $3 $4
    tar -czf re$i.tar.gz ms_*
    if [ -e mags ]; then
	mv mags mags$i
    fi

    #visualize and query
    $VISIT -nowin -cli -s tmp.py >> maxpps$i
    pushd $5
    export j=1
    for file in `ls -1 visit????.png`
    do
	mv $file re$i.$j.png
	j=$((j+1))
    done
    #files are moved .... don't need this
    #rm visit????.png
    popd
done
