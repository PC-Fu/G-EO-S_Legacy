#!/bin/bash

#ARGS (1) file 1 (2) file 2
if [ $# -ne 1 ]; then
   echo "usage: ARGS (1) XML file"
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
if [ -z "$XMLVALIDATOR_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export XMLVALIDATOR_EXEC=~/apps/gpac_trunk/external_libs/bin/xmllint
    if [ -e "$XMLVALIDATOR_EXEC" ];then
	echo "$0 using XMLVALIDATOR_EXEC $XMLVALIDATOR_EXEC"
    else
	export XMLVALIDATOR_EXEC=/usr/gapps/GEOS/external_libs_icc_dbg/bin/xmllint
	if [ -e "$XMLVALIDATOR_EXEC" ];then
	    echo "$0 using XMLVALIDATOR_EXEC $XMLVALIDATOR_EXEC"
	else
	    echo "$0 cannot find XMLVALIDATOR_EXEC $XMLVALIDATOR_EXEC"
	    exit 1
	fi
    fi
fi

if [ -z "$XML_SCHEMA_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export XML_SCHEMA_PATH=~/apps/geos_git/gpaccore/src/schema/gpac.xsd
    if [ -e "$XML_SCHEMA_PATH" ]; then
	echo "$0 using XML_SCHEMA_PATH $XML_SCHEMA_PATH"
    else
	export XML_SCHEMA_PATH=/usr/gapps/GEOS/schema/gpac.xsd
	if [ -f "$XML_SCHEMA_PATH" ];then
	    echo "$0 using XML_SCHEMA_PATH $XML_SCHEMA_PATH"
	else
	    echo "$0 cannot find XML_SCHEMA_PATH $XML_SCHEMA_PATH"
	    exit 1
	fi
    fi
fi
if [ ! -f $1 ]; then
    echo $PWD
    echo "$1 (1) does not exist - cannot validate!"
    exit 1
fi

cp $1 xxxxxx_$1

if [ -f "$GPAC_SCRIPT_PATH/xml_insert_include.pl" ]; then
echo "...including xml"
perl $GPAC_SCRIPT_PATH/xml_insert_include.pl xxxxxx_$1 > $1
else
echo "cannot find $GPAC_SCRIPT_PATH/xml_insert_include.pl"
fi

$XMLVALIDATOR_EXEC --noout --schema $XML_SCHEMA_PATH $1
mv xxxxxx_$1 $1

