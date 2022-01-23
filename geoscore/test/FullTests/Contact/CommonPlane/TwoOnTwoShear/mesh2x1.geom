*HEADING
cubit(tmp.inp): 10/04/2011: 07:57:17
version: 13.0
**
********************************** P A R T S **********************************
*PART, NAME=Part-Default
**
********************************** N O D E S **********************************
*NODE, NSET=ALLNODES
       1, -1.000000e+00, -1.000000e+00, 1.5
       2, -1.000000e+00, -1.000000e+00, 0.5
       3, -1.000000e+00, 0.000000e+00, 0.5
       4, -1.000000e+00, 0.000000e+00, 1.5
       5, 0.000000e+00, -1.000000e+00, 1.5
       6, 0.000000e+00, -1.000000e+00, 0.5
       7, 0.000000e+00, 0.000000e+00, 0.5
       8, 0.000000e+00, 0.000000e+00, 1.5
       9, -1.000000e+00, 1.000000e+00, 0.5
      10, -1.000000e+00, 1.000000e+00, 1.5
      11, 0.000000e+00, 1.000000e+00, 0.5
      12, 0.000000e+00, 1.000000e+00, 1.5
      13, -1.000000e+00, -2.000000e+00, 0.5
      14, -1.000000e+00, -2.000000e+00, -0.5
      15, -1.000000e+00, -1.000000e+00, -0.5
      16, -1.000000e+00, -1.000000e+00, 0.5
      17, 0.000000e+00, -2.000000e+00, 0.5
      18, 0.000000e+00, -2.000000e+00, -0.5
      19, 0.000000e+00, -1.000000e+00, -0.5
      20, 0.000000e+00, -1.000000e+00, 0.5
      21, -1.000000e+00, 0.000000e+00, -0.5
      22, -1.000000e+00, 0.000000e+00, 0.5
      23, 0.000000e+00, 0.000000e+00, -0.5
      24, 0.000000e+00, 0.000000e+00, 0.5
**
********************************** E L E M E N T S ****************************
*ELEMENT, TYPE=C3D8R, ELSET=EB1
       1,       1,       2,       3,       4,       5,       6,       7,       8
       2,       4,       3,       9,      10,       8,       7,      11,      12
*ELEMENT, TYPE=C3D8R, ELSET=EB2
       3,      13,      14,      15,      16,      17,      18,      19,      20
       4,      16,      15,      21,      22,      20,      19,      23,      24
**
********************************** N O D E S E T S **********************************
*NSET, NSET=NS1
      1,        4,       5,       8,      10,      12,
*NSET, NSET=BLOCK1
 1  2  3  4  5  6  7  8  9 10 11 12
*NSET, NSET=BLOCK2
 13 14 15 16 17 18 19 20 21 22 23 24
*NSET, NSET=NS2
      14,      15,      18,      19,      21,      23,
**
********************************** P R O P E R T I E S ************************
*SOLID SECTION, ELSET=EB1, MATERIAL=Default-Steel
*SOLID SECTION, ELSET=EB2, MATERIAL=Default-Steel
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
