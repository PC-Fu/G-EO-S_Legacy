//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2014, Lawrence Livermore National Security, LLC.
//
//  Produced at the Lawrence Livermore National Laboratory
//
//  Written by:
//
//  Randolph Settgast		Stuart Walsh
//  Scott Johnson		Pengcheng Fu
//  Joshua White
//
//  LLNL-CODE-656690
//  GMOD-B, Version 1.0
//
//  All rights reserved.
//
//  This file is part of GMOD-B. For details, please contact Scott Johnson (scott.johnson@alum.mit.edu) or Randolph Settgast (rrsettgast@gmail.com).
//  Please also read "Additional BSD Notice" below.
//
//  Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright notice, this list of conditions and the disclaimer below.
//  * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the disclaimer (as noted below) in the 
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of the LLNS/LLNL nor the names of its contributors may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  Additional BSD Notice
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" {
  #include <stdio.h>

  typedef double __C_real; // could be float.
  #define C const
  #define bigReal 1.E38  // FLT_MAX, DBL_MAX if double above

  #include <limits.h>

  typedef struct{__C_real x; __C_real y;} __C_point;
  typedef struct{__C_point min; __C_point max;} __C_box;
  typedef long long __C_hp;
  typedef struct{long x; long y;} __C_ipoint;
  typedef struct{long mn; long mx;} __C_rng;
  typedef struct{__C_ipoint ip; __C_rng rx; __C_rng ry; short in;} vertex;

  static void bd(__C_real * X, __C_real y){*X = *X<y ? *X:y;}

  static void bu(__C_real * X, __C_real y){*X = *X>y ? *X:y;}

  static void range(__C_box& B, __C_point * x, int c){
    while(c--){
      bd(&B.min.x, x[c].x); bu(&B.max.x, x[c].x);
      bd(&B.min.y, x[c].y); bu(&B.max.y, x[c].y);
    }
  }

  static void fit(__C_box& B, __C_point * x, int cx, vertex * ix, __C_real sclx, __C_real scly, __C_real mid, int fudge)
  {
    {
      int c=cx;
      while(c--){
        ix[c].ip.x = (((long)((x[c].x - B.min.x)*sclx - mid)&(~7))|fudge)|(c&1);//--
        ix[c].ip.y = ((long)((x[c].y - B.min.y)*scly - mid)&(~7))|fudge;//--
      }
    }
    ix[0].ip.y += cx&1;
    ix[cx] = ix[0];
    {
      int c=cx; while(c--) {

        __C_rng temp1a = {ix[c].ip.x,ix[c+1].ip.x};
        __C_rng temp1b = {ix[c+1].ip.x,ix[c].ip.x};
        ix[c].rx = ix[c].ip.x < ix[c+1].ip.x ? temp1a : temp1b ;

        __C_rng temp2a = {ix[c].ip.y,ix[c+1].ip.y};
        __C_rng temp2b = {ix[c+1].ip.y,ix[c].ip.y};
        ix[c].ry = ix[c].ip.y < ix[c+1].ip.y ? temp2a : temp2b ;
        ix[c].in=0;
      }
    }
  }

  static __C_hp area(__C_ipoint a, __C_ipoint p, __C_ipoint q)
  {
    return (__C_hp)p.x*q.y - (__C_hp)p.y*q.x + (__C_hp)a.x*(p.y - q.y) + (__C_hp)a.y*(q.x - p.x);
  }

  static int ovl(__C_rng p, __C_rng q){return (p.mn < q.mx) && (q.mn < p.mx);}

  static void cntrib(__C_ipoint f, __C_ipoint t, short w, __C_hp& s)
  {
    //printf("s(48)=%Ld -> ", s);
    s+=(__C_hp)w*(t.x-f.x)*(t.y+f.y)/2;
    //printf("s(55)=%Ld t.x=%ld f.x=%ld t.y=%ld f.y=%ld\n", s, t.x, t.y, f.x, f.y);
  }

  static void inness(vertex * P, int cP, vertex * Q, int cQ, __C_hp& ss)
  {
    __C_hp s=0;
    int c=cQ;
    __C_ipoint p = P[0].ip;
    while(c--)
      if(Q[c].rx.mn < p.x && p.x < Q[c].rx.mx)
      {
        int sgn = 0 < area(p, Q[c].ip, Q[c+1].ip);
        s += (sgn != (Q[c].ip.x < Q[c+1].ip.x)) ? 0 : (sgn?-1:1);
      }
    {
      int j;
      for(j=0; j<cP; ++j)
      {
        if(s)
          cntrib(P[j].ip, P[j+1].ip, s, ss);
        s += P[j].in; //printf("s(74)=%Ld\n", s);
      }
    }
  }

  static void cross(vertex * a, vertex * b, vertex * c, vertex * d,
                   double a1, double a2, double a3, double a4, __C_hp& s)
  {
    __C_real r1=a1/((__C_real)a1+a2);
    __C_real r2 = a3/((__C_real)a3+a4);

    __C_ipoint temp1 = {a->ip.x + static_cast<long>(r1*(b->ip.x-a->ip.x)),
                        a->ip.y + static_cast<long>(r1*(b->ip.y-a->ip.y))};
    cntrib( temp1, b->ip, 1, s);

    __C_ipoint temp2 = {c->ip.x + static_cast<long>(r2*(d->ip.x - c->ip.x)),
                        c->ip.y + static_cast<long>(r2*(d->ip.y - c->ip.y))};
    cntrib(d->ip, temp2, 1, s);
    ++a->in;
    --c->in;
  }

  __C_real inter(__C_point * a, int na, __C_point * b, int nb){

//    vertex ipa[na + 1];
//    vertex ipb[nb + 1];
    // variable length arrays are NOT part of ISO c++11
    vertex ipa[20];
    vertex ipb[20];

    __C_box B =
    {
    { bigReal, bigReal },
    { -bigReal, -bigReal } };
    double ascale;
    if (na < 3 || nb < 3)
      return 0;
    range(B, a, na);
    range(B, b, nb);
    {
      const __C_real gamut = 500000000., mid = gamut / 2.;
      __C_real rngx = B.max.x - B.min.x, sclx = gamut / rngx, rngy = B.max.y - B.min.y, scly = gamut
          / rngy;
      fit(B, a, na, ipa, sclx, scly, mid, 0);
      fit(B, b, nb, ipb, sclx, scly, mid, 2);
      ascale = sclx * scly;
    }
    {
      __C_hp s = 0;
      int j, k;
      for (j = 0; j < na; ++j)
        for (k = 0; k < nb; ++k)
          if (ovl(ipa[j].rx, ipb[k].rx) && ovl(ipa[j].ry, ipb[k].ry))
          {
            __C_hp a1 = -area(ipa[j].ip, ipb[k].ip, ipb[k + 1].ip), a2 = area(ipa[j + 1].ip, ipb[k].ip,
                                                                          ipb[k + 1].ip);
            {
              int o = a1 < 0;
              if (o == (a2 < 0))
              {
                __C_hp a3 = area(ipb[k].ip, ipa[j].ip, ipa[j + 1].ip), a4 = -area(ipb[k + 1].ip,
                                                                              ipa[j].ip,
                                                                              ipa[j + 1].ip);
                if ((a3 < 0) == (a4 < 0))
                {
                  if (o)
                    cross(&ipa[j], &ipa[j + 1], &ipb[k], &ipb[k + 1], a1, a2, a3, a4, s);
                  else
                    cross(&ipb[k], &ipb[k + 1], &ipa[j], &ipa[j + 1], a3, a4, a1, a2, s);
                }
              }
            }
          }
      {
        inness(ipa, na, ipb, nb, s);
        inness(ipb, nb, ipa, na, s);
      }
      //printf("s(141)=%Ld\n", s);
      return s / ascale;
    }
  }
}
