using DoubleDouble;
using System;
using System.Collections.ObjectModel;

namespace DoubleDoubleIntegrate {
    public static class GaussLegendreIntegral {
        public static ddouble Integrate(Func<ddouble, ddouble> f, ddouble a, ddouble b, int n) {
            ddouble r = b - a;

            if (!ddouble.IsFinite(r)) {
                throw new ArgumentException("Invalid integation interval.", $"{nameof(a)},{nameof(b)}");
            }

            if (n < GaussLegendrePoints.MinPoints || n > GaussLegendrePoints.MaxPoints) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }

            ReadOnlyCollection<(ddouble x, ddouble w)> ps = GaussLegendrePoints.Table[n];

            ddouble s = ddouble.Zero;

            foreach ((ddouble x, ddouble w) in ps) {
                ddouble x_shifted = x * r + a;

                s += w * f(x_shifted);
            }

            s *= r;

            return s;
        }
    }
}
