#!/bin/bash
if [ -z "$GPAC_SCRIPT_PATH" ]; then
#NOTE: this only provides a value for the current script
export GPAC_SCRIPT_PATH=~/apps/gpac_trunk/scripts
echo "setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
fi

if [ "1" -eq "1" ]; then

##############################################################
cd Material
bash make_materials.sh
perl $GPAC_SCRIPT_PATH/gpac_constitutive_model_query_class.pl Materials make_materials.sh
perl $GPAC_SCRIPT_PATH/gpac_extract_xml.pl MaterialType make_materials.sh > ../tmp.xsd
cd ../

##############################################################
cd Interface
bash make_interfaces.sh
perl $GPAC_SCRIPT_PATH/gpac_constitutive_model_query_class.pl Interfaces make_interfaces.sh
perl $GPAC_SCRIPT_PATH/gpac_extract_xml.pl ContactType make_interfaces.sh >> ../tmp.xsd
cd ../

##############################################################
perl -e 'open(IN,"<../schema/gpac.xsd");my $xon = 0; while(<IN>){my $line = $_; if(/BEGIN=====/){$xon = 1; print $line;}elsif(/END=====/){$xon = 0; open(INSERT, "<tmp.xsd"); while(<INSERT>) {print $_;} close INSERT; print $line;}elsif($xon == 0){print $line;}} close IN;' > tmp2.xsd
mv tmp2.xsd ../schema/gpac.xsd

fi

##############################################################
cd CohesiveZone
bash make_cohesive_zones.sh
perl $GPAC_SCRIPT_PATH/gpac_constitutive_model_query_class.pl CohesiveZones make_cohesive_zones.sh
cd ../
