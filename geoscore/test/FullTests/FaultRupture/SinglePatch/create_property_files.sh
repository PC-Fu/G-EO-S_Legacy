#!/usr/bash

perl create_property_files.pl 0 -1080 -10 80 -1000 10 frictionCoefficient 0.7 >> SinglePatch.xml
perl create_property_files.pl 0 -1080 -10 80 -1000 10 currentFrictionCoefficient 0.7 >> SinglePatch.xml
perl create_property_files.pl 0 -1080 -10 80 -1000 10 A 0.0065 >> SinglePatch.xml
perl create_property_files.pl 0 -1080 -10 80 -1000 10 B 0.015 >> SinglePatch.xml
perl create_property_files.pl 0 -1080 -10 80 -1000 10 Dc 0.0000225 >> SinglePatch.xml
perl create_property_files.pl 0 -1080 -10 80 -1000 10 alpha 0.25 >> SinglePatch.xml
perl create_property_files.pl 0 -1080 -10 80 -1000 10 shearRateStar 0.000001 >> SinglePatch.xml
perl create_property_files.pl 0 -1080 -10 80 -1000 10 shearSlipRateAB 0.000001 >> SinglePatch.xml
perl create_property_files.pl 0 -1080 -10 80 -1000 10 stateStar 0.1 >> SinglePatch.xml
