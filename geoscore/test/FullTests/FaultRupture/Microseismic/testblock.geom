*HEADING
                                                                                
*NODE,NSET=PART1
1,-1.,-1.,-1.
2,-1.,-1., 1.
3,-1., 1.,-1.
4,-1., 1., 1.
5, 1.,-1.,-1.
6, 1.,-1., 1.
7, 1., 1.,-1.
8, 1., 1., 1.
**** WARNING - MATERIAL    1 UNDEFINED FOR THE FOLLOWING ELEMENTS ****
*ORIENTATION,NAME=SOR1,DEFINITION=COORDINATES
0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00
0,0.000000E+00
*SOLID SECTION,ELSET=PM1,MATERIAL=unknown,ORIENTATION=SOR1
*ELEMENT,TYPE=C3D8,ELSET=PM1
1,1,5,7,3,2,6,8,4
*NSET,NSET=xneg
1,2,3,4
*NSET,NSET=xpos
5,6,7,8
