{
  /*
   * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
   */
  const double xcoord = xcoordGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
  const double ycoord = ycoordGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
  const double zcoord = zcoordGF[CCTK_GFINDEX3D(cctkGH, i0,i1,i2)];
  /*
   * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
   */
  const double FDPart3_0 = ((a)*(a));
  const double FDPart3_1 = ((zcoord)*(zcoord));
  const double FDPart3_2 = ((xcoord)*(xcoord));
  const double FDPart3_3 = ((ycoord)*(ycoord));
  const double FDPart3_4 = FDPart3_2 + FDPart3_3;
  const double FDPart3_5 = FDPart3_1 + FDPart3_4;
  const double FDPart3_6 = (1.0/(FDPart3_5));
  const double FDPart3_7 = FDPart3_1*FDPart3_6;
  const double FDPart3_8 = FDPart3_0*FDPart3_7 + FDPart3_5;
  const double FDPart3_9 = (1.0/(FDPart3_8));
  const double FDPart3_10 = sqrt(FDPart3_5);
  const double FDPart3_12 = 2*FDPart3_10*M;
  const double FDPart3_14 = FDPart3_12*FDPart3_9 + 1;
  const double FDPart3_15 = (1.0/(FDPart3_14));
  const double FDPart3_16 = M*xcoord;
  const double FDPart3_20 = 2*FDPart3_15*FDPart3_9;
  const double FDPart3_22 = 1 - FDPart3_7;
  const double FDPart3_24 = FDPart3_0*((FDPart3_14)*(FDPart3_14))*((FDPart3_22)*(FDPart3_22));
  const double FDPart3_25 = ((FDPart3_8)*(FDPart3_8));
  const double FDPart3_26 = FDPart3_0*FDPart3_12*FDPart3_22*FDPart3_9 + FDPart3_0 + FDPart3_5;
  const double FDPart3_27 = FDPart3_14*FDPart3_22*FDPart3_26*FDPart3_8 - FDPart3_24*FDPart3_8;
  const double FDPart3_28 = (1.0/((FDPart3_27)*(FDPart3_27)*(FDPart3_27)));
  const double FDPart3_29 = FDPart3_14*FDPart3_22*FDPart3_26 - FDPart3_24;
  const double FDPart3_30 = (1.0/(FDPart3_14*FDPart3_22*FDPart3_25*FDPart3_26*FDPart3_28*FDPart3_29 - FDPart3_24*FDPart3_25*FDPart3_28*FDPart3_29));
  const double FDPart3_31 = (1.0/((FDPart3_27)*(FDPart3_27)));
  const double FDPart3_32 = FDPart3_29*FDPart3_30*FDPart3_31*FDPart3_8;
  const double FDPart3_33 = FDPart3_14*FDPart3_32;
  const double FDPart3_34 = xcoord*ycoord;
  const double FDPart3_35 = (1.0/(FDPart3_4));
  const double FDPart3_36 = FDPart3_22*FDPart3_35;
  const double FDPart3_37 = (1.0/(FDPart3_10));
  const double FDPart3_39 = -FDPart3_14*FDPart3_37*a;
  const double FDPart3_41 = 2*FDPart3_32*FDPart3_34*FDPart3_36*FDPart3_39;
  const double FDPart3_43 = (1.0/((FDPart3_4)*(FDPart3_4)));
  const double FDPart3_44 = FDPart3_22*FDPart3_26*FDPart3_43;
  const double FDPart3_46 = (1.0/(FDPart3_22));
  const double FDPart3_47 = FDPart3_30*(FDPart3_14*FDPart3_22*FDPart3_25*FDPart3_26*FDPart3_31 - FDPart3_24*FDPart3_25*FDPart3_31);
  const double FDPart3_49 = FDPart3_46*FDPart3_47/((FDPart3_5)*(FDPart3_5)*(FDPart3_5));
  const double FDPart3_50 = FDPart3_33*FDPart3_6;
  const double FDPart3_54 = FDPart3_37*zcoord;
  const double FDPart3_56 = -FDPart3_14*FDPart3_32*FDPart3_36*FDPart3_54*a;
  const double FDPart3_57 = pow(FDPart3_5, -3.0/2.0);
  const double FDPart3_58 = FDPart3_1*FDPart3_57;
  const double FDPart3_59 = FDPart3_37 - FDPart3_58;
  const double FDPart3_60 = FDPart3_46*FDPart3_47*FDPart3_57*FDPart3_59;
  const double FDPart3_62 = FDPart3_46*((FDPart3_59)*(FDPart3_59));
  const double FDPart3_63 = (1.0/((FDPart3_5)*(FDPart3_5)));
  const double FDPart3_64 = FDPart3_1*FDPart3_2*FDPart3_63;
  const double FDPart3_65 = sqrt(FDPart3_14);
  const double FDPart3_67 = acos(FDPart3_54);
  const double FDPart3_68 = cos(2*FDPart3_67);
  const double FDPart3_70 = 2*FDPart3_2;
  const double FDPart3_71 = 2*FDPart3_3;
  const double FDPart3_72 = 2*FDPart3_1 + FDPart3_70 + FDPart3_71;
  const double FDPart3_73 = FDPart3_0*FDPart3_68 + FDPart3_0 + FDPart3_72;
  const double FDPart3_74 = (1.0/(4*FDPart3_10*M + FDPart3_73));
  const double FDPart3_75 = FDPart3_65*FDPart3_74;
  const double FDPart3_76 = FDPart3_75*M;
  const double FDPart3_78 = 4*FDPart3_46*FDPart3_76;
  const double FDPart3_79 = (1.0/(FDPart3_73));
  const double FDPart3_81 = 16*FDPart3_0*FDPart3_79;
  const double FDPart3_82 = FDPart3_76*FDPart3_81;
  const double FDPart3_83 = (1.0/((FDPart3_73)*(FDPart3_73)));
  const double FDPart3_84 = FDPart3_0*FDPart3_68 + FDPart3_0 - FDPart3_72;
  const double FDPart3_86 = FDPart3_36*FDPart3_65*FDPart3_83*FDPart3_84;
  const double FDPart3_87 = FDPart3_37*FDPart3_86*a;
  const double FDPart3_88 = 4*FDPart3_16*FDPart3_87*ycoord;
  const double FDPart3_89 = ((a)*(a)*(a));
  const double FDPart3_90 = FDPart3_16*FDPart3_75;
  const double FDPart3_91 = FDPart3_90*ycoord;
  const double FDPart3_92 = 16*FDPart3_36*FDPart3_58*FDPart3_79*FDPart3_89*FDPart3_91;
  const double FDPart3_94 = FDPart3_83*FDPart3_84*(FDPart3_12 + FDPart3_73);
  const double FDPart3_95 = 4*FDPart3_76*FDPart3_94;
  const double FDPart3_96 = ((a)*(a)*(a)*(a));
  const double FDPart3_98 = FDPart3_22*FDPart3_43*FDPart3_83*(4*FDPart3_0*FDPart3_10*FDPart3_68*(FDPart3_0 + FDPart3_10*(2*FDPart3_10 + M)) + 4*FDPart3_0*FDPart3_5*(2*FDPart3_10 - M) + 8*pow(FDPart3_5, 5.0/2.0) + FDPart3_96*(FDPart3_10 - M)*cos(4*FDPart3_67) + FDPart3_96*(3*FDPart3_10 + M));
  const double FDPart3_99 = FDPart3_10*FDPart3_76*FDPart3_98;
  const double FDPart3_102 = FDPart3_1*FDPart3_63*FDPart3_91;
  const double FDPart3_104 = 8*FDPart3_79*FDPart3_89;
  const double FDPart3_105 = FDPart3_104*FDPart3_58*FDPart3_76;
  const double FDPart3_106 = 4*FDPart3_6*FDPart3_94;
  const double FDPart3_107 = ((zcoord)*(zcoord)*(zcoord));
  const double FDPart3_108 = FDPart3_54*FDPart3_59;
  const double FDPart3_109 = 4*FDPart3_108*FDPart3_46;
  const double FDPart3_110 = 8*FDPart3_0*FDPart3_79;
  const double FDPart3_111 = FDPart3_75*M*ycoord;
  const double FDPart3_112 = FDPart3_104*FDPart3_36*FDPart3_59*zcoord;
  const double FDPart3_113 = FDPart3_1*FDPart3_3*FDPart3_63;
  alphaGF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = sqrt(FDPart3_15);
  betaU0GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 2*FDPart3_15*FDPart3_16*FDPart3_9;
  betaU1GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_20*M*ycoord;
  betaU2GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_20*M*zcoord;
  gammaDD00GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_1*FDPart3_2*FDPart3_49 + FDPart3_2*FDPart3_33*FDPart3_6 + FDPart3_3*FDPart3_32*FDPart3_44 - FDPart3_41;
  gammaDD01GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_1*FDPart3_34*FDPart3_49 + FDPart3_2*FDPart3_32*FDPart3_36*FDPart3_39 - FDPart3_3*FDPart3_32*FDPart3_36*FDPart3_39 - FDPart3_32*FDPart3_34*FDPart3_44 + FDPart3_34*FDPart3_50;
  gammaDD02GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_14*FDPart3_29*FDPart3_30*FDPart3_31*FDPart3_6*FDPart3_8*xcoord*zcoord - FDPart3_56*ycoord - FDPart3_60*xcoord*zcoord;
  gammaDD11GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_1*FDPart3_3*FDPart3_49 + FDPart3_2*FDPart3_32*FDPart3_44 + FDPart3_3*FDPart3_50 + FDPart3_41;
  gammaDD12GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_50*ycoord*zcoord + FDPart3_56*xcoord - FDPart3_60*ycoord*zcoord;
  gammaDD22GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_33*FDPart3_7 + FDPart3_47*FDPart3_62;
  KDD00GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_2*FDPart3_6*FDPart3_95 + FDPart3_64*FDPart3_78 + FDPart3_64*FDPart3_82 + FDPart3_71*FDPart3_99 + FDPart3_88 + FDPart3_92;
  KDD01GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 4*FDPart3_102*FDPart3_46 + FDPart3_102*FDPart3_81 - FDPart3_105*FDPart3_2*FDPart3_36 + FDPart3_105*FDPart3_3*FDPart3_36 + FDPart3_106*FDPart3_91 - FDPart3_12*FDPart3_34*FDPart3_75*FDPart3_98 - FDPart3_70*FDPart3_87*M + FDPart3_71*FDPart3_87*M;
  KDD02GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = 8*FDPart3_0*FDPart3_107*FDPart3_63*FDPart3_65*FDPart3_74*FDPart3_79*M*xcoord - FDPart3_108*FDPart3_110*FDPart3_90 - FDPart3_109*FDPart3_90 - FDPart3_111*FDPart3_112 + 2*FDPart3_22*FDPart3_35*FDPart3_37*FDPart3_65*FDPart3_83*FDPart3_84*M*a*ycoord*zcoord + 4*FDPart3_6*FDPart3_65*FDPart3_74*FDPart3_83*FDPart3_84*M*xcoord*zcoord*(FDPart3_12 + FDPart3_73);
  KDD11GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_113*FDPart3_78 + FDPart3_113*FDPart3_82 + FDPart3_3*FDPart3_6*FDPart3_95 + FDPart3_70*FDPart3_99 - FDPart3_88 - FDPart3_92;
  KDD12GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = FDPart3_106*FDPart3_111*zcoord + FDPart3_107*FDPart3_110*FDPart3_111*FDPart3_63 - FDPart3_108*FDPart3_110*FDPart3_111 - FDPart3_109*FDPart3_111 + FDPart3_112*FDPart3_90 - 2*FDPart3_16*FDPart3_54*FDPart3_86*a;
  KDD22GF[CCTK_GFINDEX3D(cctkGH, i0, i1, i2)] = -FDPart3_1*FDPart3_37*FDPart3_59*FDPart3_82 + 4*FDPart3_5*FDPart3_62*FDPart3_76 + FDPart3_7*FDPart3_95;
}