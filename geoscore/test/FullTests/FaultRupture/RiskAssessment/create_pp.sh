#!/bin/bash
perl ~/apps/gpac_trunk/scripts/gpac_pore_pressure_file.pl 700 1 4000 3000 -6000 3000 4000 500 pscale tpp xpp ypp zpp > pp
perl ~/apps/gpac_trunk/scripts/gpac_pore_pressure_file.pl 700 1 4000 3000 -6000 3000 4000 500 pscale_fault tpp_fault xpp ypp zpp > pp_fault
