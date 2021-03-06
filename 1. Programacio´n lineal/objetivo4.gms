VARIABLES
          X1, X2, X3, X4, X5, X6
          X_1, X_2, X_3, X_4, X_5, X_6
          BENEFICIOS
          INGRESOS
          COSTES_PRODUCCION
          COSTES_REFINADO
          PRODUCTO
          PRODUCTO_SIN_REFINAR
          COSTES_ALMACENES

POSITIVE VARIABLES
X1, X2, X3, X4, X5, X6, X_1, X_2, X_3, X_4, X_5, X_6;

BINARY VARIABLES
DELTA1, DELTA2, DELTA3, DELTA4, DELTA5, DELTA6
DELTA_I1, DELTA_I2, DELTA_I3, DELTA_I4, DELTA_I5, DELTA_I6
DELTA_II1, DELTA_II2, DELTA_II3, DELTA_II4, DELTA_II5, DELTA_II6;


EQUATIONS
OBJ, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R_PRODUCTO, R_PRODUCTO_SIN_REFINAR, R_INGRESOS, R_COSTES_PRODUCCION, R_COSTES_REFINADO
R_MAX_TRES_ACEITES
R_MIN_15_TON_11
R_MIN_15_TON_12
R_MIN_15_TON_21
R_MIN_15_TON_22
R_MIN_15_TON_31
R_MIN_15_TON_32
R_MIN_15_TON_41
R_MIN_15_TON_42
R_MIN_15_TON_51
R_MIN_15_TON_52
R_MIN_15_TON_61
R_MIN_15_TON_62
R_SELECION_1
R_SELECION_2
R_SELECION_3
R_SELECION_4
R_SELECION_5
R_SELECION_6
R_IF_THEN
R_COSTES_ALMACENES;

R1..  X1 =E= 0.95*X_1;
R2..  X2 =E= 0.95*X_2;
R3..  X3 =E= 0.95*X_3;
R4..  X4 =E= 0.95*X_4;
R5..  X5 =E= 0.95*X_5;
R6..  X6 =E= 0.95*X_6;


R_PRODUCTO..               PRODUCTO =E= X1 + X2 + X3 + X4 + X5 + X6;
R_PRODUCTO_SIN_REFINAR..   PRODUCTO_SIN_REFINAR =E= X_1 + X_2 + X_3 + X_4 + X_5 + X_6;
R_INGRESOS..               INGRESOS =E= 150*PRODUCTO;
R_COSTES_PRODUCCION..      COSTES_PRODUCCION =E= 115*X_1 + 120*X_2 + 115*X_3 + 120*X_4 + 114*X_5 + 115*X_6;
R_COSTES_REFINADO..        COSTES_REFINADO =E= 5*PRODUCTO_SIN_REFINAR;
R_COSTES_ALMACENES..       COSTES_ALMACENES =E= (650*DELTA_I1 + 450*DELTA_II1) + (650*DELTA_I2 + 450*DELTA_II2) + (650*DELTA_I3 + 450*DELTA_II3) + (650*DELTA_I4 + 450*DELTA_II4) + (650*DELTA_I5 + 450*DELTA_II5) + (650*DELTA_I6 + 450*DELTA_II6);

R7..  8.8*X1+6.1*X2+7.5*X3+2*X4+5.2*X5+4.9*X6 -3*PRODUCTO =G= 0;
R8..  8.8*X1+6.1*X2+7.5*X3+2*X4+5.2*X5+4.9*X6 -6*PRODUCTO =L= 0;
R9..  X_1 + X_2 + X_3 =L= 225;
R10.. X_4 + X_5 + X_6 =L= 350;

R_MAX_TRES_ACEITES..  DELTA1 + DELTA2 + DELTA3 + DELTA4 + DELTA5 + DELTA6 =L= 3;

R_MIN_15_TON_11..     15*DELTA1 =L= X_1;
R_MIN_15_TON_12..     X_1 =L= 160*DELTA_I1 + 80*DELTA_II1;
R_MIN_15_TON_21..     15*DELTA2 =L= X_2;
R_MIN_15_TON_22..     X_2=L= 160*DELTA_I2 + 80*DELTA_II2;
R_MIN_15_TON_31..     15*DELTA3 =L= X_3;
R_MIN_15_TON_32..     X_3 =L= 160*DELTA_I3 + 80*DELTA_II3;

R_MIN_15_TON_41..     15*DELTA4 =L= X_4;
R_MIN_15_TON_42..     X_4 =L= 160*DELTA_I4 + 80*DELTA_II4;
R_MIN_15_TON_51..     15*DELTA5 =L= X_5;
R_MIN_15_TON_52..     X_5 =L= 160*DELTA_I5 + 80*DELTA_II5;
R_MIN_15_TON_61..     15*DELTA6 =L= X_6;
R_MIN_15_TON_62..     X_6 =L= 160*DELTA_I6 + 80*DELTA_II6;

R_IF_THEN..           DELTA1 + DELTA2 + DELTA3 =L= 3*(1-DELTA6);

R_SELECION_1..        DELTA_I1 + DELTA_II1 =E= DELTA1;
R_SELECION_2..        DELTA_I2 + DELTA_II2 =E= DELTA2;
R_SELECION_3..        DELTA_I3 + DELTA_II3 =E= DELTA3;
R_SELECION_4..        DELTA_I4 + DELTA_II4 =E= DELTA4;
R_SELECION_5..        DELTA_I5 + DELTA_II5 =E= DELTA5;
R_SELECION_6..        DELTA_I6 + DELTA_II6 =E= DELTA6;


OBJ.. BENEFICIOS =E= INGRESOS - COSTES_PRODUCCION - COSTES_REFINADO - COSTES_ALMACENES


MODEL LINEAL /ALL/;
SOLVE LINEAL USING MIP MAXIMIZING BENEFICIOS;
