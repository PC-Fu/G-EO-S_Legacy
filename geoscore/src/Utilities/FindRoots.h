/*
 * FindRoots.h
 *
 *  Created on: Jun 4, 2013
 *      Author: johnson346
 */

#ifndef FINDROOTS_H_
#define FINDROOTS_H_

class FindRoots
{
public:
  static realT FindRoot(realT(*f)(const realT, const realT *),
                        const realT* params,
                        const realT lower,
                        const realT upper,
                        const realT tol)
  {
    const localIndex MAXITERATIONS = 100;

    realT a = lower;
    realT b = upper;
    realT fa = f(a, params);
    realT fb = f(b, params);
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    {
      return std::numeric_limits<realT>::max();
    }

    realT fc = fb;
    realT c = 0.0;
    realT d = 0.0;
    realT e = 0.0;
    for (localIndex iter = 0; iter < MAXITERATIONS; iter++)
    {
      if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
      {
        c = a;
        fc = fa;
        e = d = b - a;
      }
      if (fabs(fc) < fabs(fb))
      {
        a = b;
        b = c;
        c = a;
        fa = fb;
        fb = fc;
        fc = fa;
      }

      //Set convergence tolerance
      const realT tol1 = std::numeric_limits<realT>::min() * 2.0 * fabs(b) + 0.5 * tol;

      const realT xm = 0.5 * (c - b);

      if (fabs(xm) <= tol1 || isEqual(fb, 0.0))
        return b;

      if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
      {
        const realT s = fb / fa; /* Attempt inverse quadratic interpolation. */
        realT p, q, r;
        if (isEqual(a, c))
        {
          p = 2.0 * xm * s;
          q = 1.0 - s;
        }
        else
        {
          q = fa / fc;
          r = fb / fc;
          p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
          q = (q - 1.0) * (r - 1.0) * (s - 1.0);
        }
        q = p > 0.0 ? -q : q;
        p = fabs(p);
        const realT min1 = 3.0 * xm * q - fabs(tol1 * q);
        const realT min2 = fabs(e * q);
        if (2.0 * p < (min1 < min2 ? min1 : min2))
        {
          e = d; /* Accept interpolation. */
          d = p / q;
        }
        else
        {
          d = xm; /* Interpolation failed, use bisection. */
          e = d;
        }
      }
      else
      { /* Bounds decreasing too slowly, use bisection. */
        d = xm;
        e = d;
      }
      a = b; /* Move last best guess to a. */
      fa = fb;
      if (fabs(d) > tol1) /* Evaluate new trial root. */
        b += d;
      else
        b += xm > 0 ? fabs(tol1) : -fabs(tol1);
      fb = f(b, params);
    } //for i to MAXITERATIONS
    return std::numeric_limits<realT>::max();
  }
};

#endif /* FINDROOTS_H_ */
