*HEADING
cubit(/g/g22/scottj/fault.geom.inp): 05/31/2012: 15:20:36
version: 13.0
**
********************************** P A R T S **********************************
*PART, NAME=Part-Default
**
********************************** N O D E S **********************************
*NODE, NSET=ALLNODES
     1,    80,   0,   -1080
     2,     0,   0,   -1080
     3,    80,   0,   -1000
     4,     0,   0,   -1000
     5,   240,   0,   -1080
     6,   160,   0,   -1080
     7,   240,   0,   -1000
     8,   160,   0,   -1000
**
********************************** E L E M E N T S ****************************
*ELEMENT, TYPE=S4R, ELSET=EB1
       1,       1,       2,      4,      3
       2,       6,       1,      3,      8
       3,       5,       6,      8,      7
**
********************************** P R O P E R T I E S ************************
*SHELL SECTION, ELSET=EB1, SECTION INTEGRATION=SIMPSON, MATERIAL=Default-Steel
1.000000e+00
**
*END PART
**
**
**
********************************** E N D   P A R T S **********************************
**
**
********************************** A S S E M B L Y ************************************
**
*ASSEMBLY, NAME=ASSEMBLY1
**
*INSTANCE, NAME=Part-Default_1, PART=Part-Default
*END INSTANCE
**
*END ASSEMBLY
**
**
**
*MATERIAL, NAME = Default-Steel
*ELASTIC, TYPE=ISOTROPIC
2.068000e+05, 2.900000e-01
*DENSITY
7.000000e-06
*CONDUCTIVITY,TYPE=ISO
4.500000e-02
*SPECIFIC HEAT
5.000000e+02
**
**
************************************** H I S T O R Y *************************************
**
*PREPRINT
**
**************************************** S T E P 1 ***************************************
*STEP,INC=100,NAME=Default Set
**
*STATIC
1, 1, 1e-05, 1
**
**
**
**
*END STEP
