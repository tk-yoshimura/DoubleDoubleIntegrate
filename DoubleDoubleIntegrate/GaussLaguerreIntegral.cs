using DoubleDouble;
using System;
using System.Collections.ObjectModel;

namespace DoubleDoubleIntegrate {
    public static class GaussLaguerreIntegral {
        public static ddouble Integrate(Func<ddouble, ddouble> f, int n, bool f_expscaled = false) {
            if (n < GaussLegendrePoints.MinPoints || n > GaussLegendrePoints.MaxPoints) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }

            ReadOnlyCollection<(ddouble x, ddouble w, ddouble wexp)> ps = GaussLaguerrePoints.Table[n];

            ddouble s = ddouble.Zero;

            if (f_expscaled) {
                foreach ((ddouble x, ddouble w, _) in ps) {
                    s += w * f(x);
                }
            }
            else {
                foreach ((ddouble x, _, ddouble w) in ps) {
                    s += w * f(x);
                }
            }

            return s;
        }
    }
}
