/*
 * Copyright (c) 1997-1999, 2003 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Mon Mar 24 02:08:02 EST 2003 */

#include "fftw-int.h"
#include "fftw.h"

/* Generated by: /homee/stevenj/cvs/fftw/gensrc/genfft -magic-alignment-check -magic-twiddle-load-all -magic-variables 4 -magic-loopi -hc2hc-forward 8 */

/*
 * This function contains 108 FP additions, 44 FP multiplications,
 * (or, 90 additions, 26 multiplications, 18 fused multiply/add),
 * 29 stack variables, and 64 memory accesses
 */
static const fftw_real K382683432 =
FFTW_KONST(+0.382683432365089771728459984030398866761344562);
static const fftw_real K923879532 =
FFTW_KONST(+0.923879532511286756128183189396788286822416626);
static const fftw_real K707106781 =
FFTW_KONST(+0.707106781186547524400844362104849039284835938);

/*
 * Generator Id's : 
 * $Id: fhf_8.c,v 1.4 2012-03-13 22:12:28 wwriggers Exp $
 * $Id: fhf_8.c,v 1.4 2012-03-13 22:12:28 wwriggers Exp $
 * $Id: fhf_8.c,v 1.4 2012-03-13 22:12:28 wwriggers Exp $
 */

void fftw_hc2hc_forward_8(fftw_real *A, const fftw_complex *W,
			  int iostride, int m, int dist)
{
     int i;
     fftw_real *X;
     fftw_real *Y;
     X = A;
     Y = A + (8 * iostride);
     {
	  fftw_real tmp105;
	  fftw_real tmp109;
	  fftw_real tmp115;
	  fftw_real tmp121;
	  fftw_real tmp108;
	  fftw_real tmp118;
	  fftw_real tmp112;
	  fftw_real tmp120;
	  ASSERT_ALIGNED_DOUBLE;
	  {
	       fftw_real tmp103;
	       fftw_real tmp104;
	       fftw_real tmp113;
	       fftw_real tmp114;
	       ASSERT_ALIGNED_DOUBLE;
	       tmp103 = X[0];
	       tmp104 = X[4 * iostride];
	       tmp105 = tmp103 + tmp104;
	       tmp109 = tmp103 - tmp104;
	       tmp113 = X[7 * iostride];
	       tmp114 = X[3 * iostride];
	       tmp115 = tmp113 - tmp114;
	       tmp121 = tmp113 + tmp114;
	  }
	  {
	       fftw_real tmp106;
	       fftw_real tmp107;
	       fftw_real tmp110;
	       fftw_real tmp111;
	       ASSERT_ALIGNED_DOUBLE;
	       tmp106 = X[2 * iostride];
	       tmp107 = X[6 * iostride];
	       tmp108 = tmp106 + tmp107;
	       tmp118 = tmp106 - tmp107;
	       tmp110 = X[iostride];
	       tmp111 = X[5 * iostride];
	       tmp112 = tmp110 - tmp111;
	       tmp120 = tmp110 + tmp111;
	  }
	  {
	       fftw_real tmp119;
	       fftw_real tmp122;
	       fftw_real tmp116;
	       fftw_real tmp117;
	       ASSERT_ALIGNED_DOUBLE;
	       X[2 * iostride] = tmp105 - tmp108;
	       tmp119 = tmp105 + tmp108;
	       tmp122 = tmp120 + tmp121;
	       X[4 * iostride] = tmp119 - tmp122;
	       X[0] = tmp119 + tmp122;
	       Y[-2 * iostride] = tmp121 - tmp120;
	       tmp116 = K707106781 * (tmp112 + tmp115);
	       X[3 * iostride] = tmp109 - tmp116;
	       X[iostride] = tmp109 + tmp116;
	       tmp117 = K707106781 * (tmp115 - tmp112);
	       Y[-iostride] = tmp117 - tmp118;
	       Y[-3 * iostride] = tmp118 + tmp117;
	  }
     }
     X = X + dist;
     Y = Y - dist;
     for (i = 2; i < m; i = i + 2, X = X + dist, Y = Y - dist, W = W + 7) {
	  fftw_real tmp29;
	  fftw_real tmp65;
	  fftw_real tmp92;
	  fftw_real tmp97;
	  fftw_real tmp63;
	  fftw_real tmp75;
	  fftw_real tmp78;
	  fftw_real tmp87;
	  fftw_real tmp40;
	  fftw_real tmp98;
	  fftw_real tmp68;
	  fftw_real tmp89;
	  fftw_real tmp52;
	  fftw_real tmp70;
	  fftw_real tmp73;
	  fftw_real tmp86;
	  ASSERT_ALIGNED_DOUBLE;
	  {
	       fftw_real tmp23;
	       fftw_real tmp91;
	       fftw_real tmp28;
	       fftw_real tmp90;
	       ASSERT_ALIGNED_DOUBLE;
	       tmp23 = X[0];
	       tmp91 = Y[-7 * iostride];
	       {
		    fftw_real tmp25;
		    fftw_real tmp27;
		    fftw_real tmp24;
		    fftw_real tmp26;
		    ASSERT_ALIGNED_DOUBLE;
		    tmp25 = X[4 * iostride];
		    tmp27 = Y[-3 * iostride];
		    tmp24 = c_re(W[3]);
		    tmp26 = c_im(W[3]);
		    tmp28 = (tmp24 * tmp25) - (tmp26 * tmp27);
		    tmp90 = (tmp26 * tmp25) + (tmp24 * tmp27);
	       }
	       tmp29 = tmp23 + tmp28;
	       tmp65 = tmp23 - tmp28;
	       tmp92 = tmp90 + tmp91;
	       tmp97 = tmp91 - tmp90;
	  }
	  {
	       fftw_real tmp57;
	       fftw_real tmp76;
	       fftw_real tmp62;
	       fftw_real tmp77;
	       ASSERT_ALIGNED_DOUBLE;
	       {
		    fftw_real tmp54;
		    fftw_real tmp56;
		    fftw_real tmp53;
		    fftw_real tmp55;
		    ASSERT_ALIGNED_DOUBLE;
		    tmp54 = X[7 * iostride];
		    tmp56 = Y[0];
		    tmp53 = c_re(W[6]);
		    tmp55 = c_im(W[6]);
		    tmp57 = (tmp53 * tmp54) - (tmp55 * tmp56);
		    tmp76 = (tmp55 * tmp54) + (tmp53 * tmp56);
	       }
	       {
		    fftw_real tmp59;
		    fftw_real tmp61;
		    fftw_real tmp58;
		    fftw_real tmp60;
		    ASSERT_ALIGNED_DOUBLE;
		    tmp59 = X[3 * iostride];
		    tmp61 = Y[-4 * iostride];
		    tmp58 = c_re(W[2]);
		    tmp60 = c_im(W[2]);
		    tmp62 = (tmp58 * tmp59) - (tmp60 * tmp61);
		    tmp77 = (tmp60 * tmp59) + (tmp58 * tmp61);
	       }
	       tmp63 = tmp57 + tmp62;
	       tmp75 = tmp57 - tmp62;
	       tmp78 = tmp76 - tmp77;
	       tmp87 = tmp76 + tmp77;
	  }
	  {
	       fftw_real tmp34;
	       fftw_real tmp66;
	       fftw_real tmp39;
	       fftw_real tmp67;
	       ASSERT_ALIGNED_DOUBLE;
	       {
		    fftw_real tmp31;
		    fftw_real tmp33;
		    fftw_real tmp30;
		    fftw_real tmp32;
		    ASSERT_ALIGNED_DOUBLE;
		    tmp31 = X[2 * iostride];
		    tmp33 = Y[-5 * iostride];
		    tmp30 = c_re(W[1]);
		    tmp32 = c_im(W[1]);
		    tmp34 = (tmp30 * tmp31) - (tmp32 * tmp33);
		    tmp66 = (tmp32 * tmp31) + (tmp30 * tmp33);
	       }
	       {
		    fftw_real tmp36;
		    fftw_real tmp38;
		    fftw_real tmp35;
		    fftw_real tmp37;
		    ASSERT_ALIGNED_DOUBLE;
		    tmp36 = X[6 * iostride];
		    tmp38 = Y[-iostride];
		    tmp35 = c_re(W[5]);
		    tmp37 = c_im(W[5]);
		    tmp39 = (tmp35 * tmp36) - (tmp37 * tmp38);
		    tmp67 = (tmp37 * tmp36) + (tmp35 * tmp38);
	       }
	       tmp40 = tmp34 + tmp39;
	       tmp98 = tmp34 - tmp39;
	       tmp68 = tmp66 - tmp67;
	       tmp89 = tmp66 + tmp67;
	  }
	  {
	       fftw_real tmp46;
	       fftw_real tmp71;
	       fftw_real tmp51;
	       fftw_real tmp72;
	       ASSERT_ALIGNED_DOUBLE;
	       {
		    fftw_real tmp43;
		    fftw_real tmp45;
		    fftw_real tmp42;
		    fftw_real tmp44;
		    ASSERT_ALIGNED_DOUBLE;
		    tmp43 = X[iostride];
		    tmp45 = Y[-6 * iostride];
		    tmp42 = c_re(W[0]);
		    tmp44 = c_im(W[0]);
		    tmp46 = (tmp42 * tmp43) - (tmp44 * tmp45);
		    tmp71 = (tmp44 * tmp43) + (tmp42 * tmp45);
	       }
	       {
		    fftw_real tmp48;
		    fftw_real tmp50;
		    fftw_real tmp47;
		    fftw_real tmp49;
		    ASSERT_ALIGNED_DOUBLE;
		    tmp48 = X[5 * iostride];
		    tmp50 = Y[-2 * iostride];
		    tmp47 = c_re(W[4]);
		    tmp49 = c_im(W[4]);
		    tmp51 = (tmp47 * tmp48) - (tmp49 * tmp50);
		    tmp72 = (tmp49 * tmp48) + (tmp47 * tmp50);
	       }
	       tmp52 = tmp46 + tmp51;
	       tmp70 = tmp46 - tmp51;
	       tmp73 = tmp71 - tmp72;
	       tmp86 = tmp71 + tmp72;
	  }
	  {
	       fftw_real tmp41;
	       fftw_real tmp64;
	       fftw_real tmp85;
	       fftw_real tmp88;
	       ASSERT_ALIGNED_DOUBLE;
	       tmp41 = tmp29 + tmp40;
	       tmp64 = tmp52 + tmp63;
	       Y[-4 * iostride] = tmp41 - tmp64;
	       X[0] = tmp41 + tmp64;
	       {
		    fftw_real tmp95;
		    fftw_real tmp96;
		    fftw_real tmp93;
		    fftw_real tmp94;
		    ASSERT_ALIGNED_DOUBLE;
		    tmp95 = tmp92 - tmp89;
		    tmp96 = tmp63 - tmp52;
		    X[6 * iostride] = -(tmp95 - tmp96);
		    Y[-2 * iostride] = tmp96 + tmp95;
		    tmp93 = tmp89 + tmp92;
		    tmp94 = tmp86 + tmp87;
		    X[4 * iostride] = -(tmp93 - tmp94);
		    Y[0] = tmp94 + tmp93;
	       }
	       tmp85 = tmp29 - tmp40;
	       tmp88 = tmp86 - tmp87;
	       Y[-6 * iostride] = tmp85 - tmp88;
	       X[2 * iostride] = tmp85 + tmp88;
	       {
		    fftw_real tmp81;
		    fftw_real tmp99;
		    fftw_real tmp84;
		    fftw_real tmp100;
		    fftw_real tmp82;
		    fftw_real tmp83;
		    ASSERT_ALIGNED_DOUBLE;
		    tmp81 = tmp65 - tmp68;
		    tmp99 = tmp97 - tmp98;
		    tmp82 = tmp73 - tmp70;
		    tmp83 = tmp75 + tmp78;
		    tmp84 = K707106781 * (tmp82 - tmp83);
		    tmp100 = K707106781 * (tmp82 + tmp83);
		    Y[-7 * iostride] = tmp81 - tmp84;
		    X[3 * iostride] = tmp81 + tmp84;
		    X[5 * iostride] = -(tmp99 - tmp100);
		    Y[-iostride] = tmp100 + tmp99;
	       }
	       {
		    fftw_real tmp69;
		    fftw_real tmp101;
		    fftw_real tmp80;
		    fftw_real tmp102;
		    fftw_real tmp74;
		    fftw_real tmp79;
		    ASSERT_ALIGNED_DOUBLE;
		    tmp69 = tmp65 + tmp68;
		    tmp101 = tmp98 + tmp97;
		    tmp74 = tmp70 + tmp73;
		    tmp79 = tmp75 - tmp78;
		    tmp80 = K707106781 * (tmp74 + tmp79);
		    tmp102 = K707106781 * (tmp79 - tmp74);
		    Y[-5 * iostride] = tmp69 - tmp80;
		    X[iostride] = tmp69 + tmp80;
		    X[7 * iostride] = -(tmp101 - tmp102);
		    Y[-3 * iostride] = tmp102 + tmp101;
	       }
	  }
     }
     if (i == m) {
	  fftw_real tmp1;
	  fftw_real tmp19;
	  fftw_real tmp4;
	  fftw_real tmp18;
	  fftw_real tmp8;
	  fftw_real tmp14;
	  fftw_real tmp11;
	  fftw_real tmp15;
	  fftw_real tmp2;
	  fftw_real tmp3;
	  ASSERT_ALIGNED_DOUBLE;
	  tmp1 = X[0];
	  tmp19 = X[4 * iostride];
	  tmp2 = X[2 * iostride];
	  tmp3 = X[6 * iostride];
	  tmp4 = K707106781 * (tmp2 - tmp3);
	  tmp18 = K707106781 * (tmp2 + tmp3);
	  {
	       fftw_real tmp6;
	       fftw_real tmp7;
	       fftw_real tmp9;
	       fftw_real tmp10;
	       ASSERT_ALIGNED_DOUBLE;
	       tmp6 = X[iostride];
	       tmp7 = X[5 * iostride];
	       tmp8 = (K923879532 * tmp6) - (K382683432 * tmp7);
	       tmp14 = (K382683432 * tmp6) + (K923879532 * tmp7);
	       tmp9 = X[3 * iostride];
	       tmp10 = X[7 * iostride];
	       tmp11 = (K382683432 * tmp9) - (K923879532 * tmp10);
	       tmp15 = (K923879532 * tmp9) + (K382683432 * tmp10);
	  }
	  {
	       fftw_real tmp5;
	       fftw_real tmp12;
	       fftw_real tmp21;
	       fftw_real tmp22;
	       ASSERT_ALIGNED_DOUBLE;
	       tmp5 = tmp1 + tmp4;
	       tmp12 = tmp8 + tmp11;
	       X[3 * iostride] = tmp5 - tmp12;
	       X[0] = tmp5 + tmp12;
	       tmp21 = tmp11 - tmp8;
	       tmp22 = tmp19 - tmp18;
	       Y[-2 * iostride] = tmp21 - tmp22;
	       Y[-iostride] = tmp21 + tmp22;
	  }
	  {
	       fftw_real tmp17;
	       fftw_real tmp20;
	       fftw_real tmp13;
	       fftw_real tmp16;
	       ASSERT_ALIGNED_DOUBLE;
	       tmp17 = tmp14 + tmp15;
	       tmp20 = tmp18 + tmp19;
	       Y[0] = -(tmp17 + tmp20);
	       Y[-3 * iostride] = tmp20 - tmp17;
	       tmp13 = tmp1 - tmp4;
	       tmp16 = tmp14 - tmp15;
	       X[2 * iostride] = tmp13 - tmp16;
	       X[iostride] = tmp13 + tmp16;
	  }
     }
}

static const int twiddle_order[] = { 1, 2, 3, 4, 5, 6, 7 };
fftw_codelet_desc fftw_hc2hc_forward_8_desc = {
     "fftw_hc2hc_forward_8",
     (void (*)()) fftw_hc2hc_forward_8,
     8,
     FFTW_FORWARD,
     FFTW_HC2HC,
     179,
     7,
     twiddle_order,
};
