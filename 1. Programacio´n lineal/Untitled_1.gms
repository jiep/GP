VARIABLES
X1, X2, X3, X4, X5, X6;

POSITIVES VARIABLES
X1, X2, X3, X4, X5, X6, Z;

EQUATIONS
OBJ, R1, R2, R3, R4

OBJ.. Z =E= 150*0.95*(x1+x2+x3+x4+x5+x6) - (115*x1+120*x2+115*x3+120*x4+114*x5+115*x6) - 5*(x1+x2+x3+x4+x5+x6);
R1..  x1 + x2 + x3 =L= 225;
R2..  x4 + x5 + x6 =L= 350;
R3..  2.8*x1+0.1*x2+1.5*x3-4*x4-0.8*x5-1.1x6 =L= 0;
R4..  -5.3*x1-3.1*x2-4.5*x3+1*x4-2.2*x5+1.9*x6 =L= 0;