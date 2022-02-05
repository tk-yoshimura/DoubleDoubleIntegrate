using DoubleDouble;
using System;
using System.Collections.ObjectModel;
using System.Linq;

namespace DoubleDoubleIntegrate {
    public static class GaussHermiteIntegral {
        public static ddouble Integrate(Func<ddouble, ddouble> f, int n, bool f_expscaled = false) {
            if (n < GaussLegendrePoints.MinPoints || n > GaussLegendrePoints.MaxPoints) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }

            ReadOnlyCollection<(ddouble x, ddouble w, ddouble wexp)> ps = GaussHermitePoints.Table[n];

            ddouble s = ddouble.Zero;

            if ((n & 1) == 0) {
                if (f_expscaled) {
                    foreach ((ddouble x, ddouble w, _) in ps) {
                        s += w * (f(x) + f(-x));
                    }
                }
                else {
                    foreach ((ddouble x, _, ddouble w) in ps) {
                        s += w * (f(x) + f(-x));
                    }
                }
            }
            else {
                if (f_expscaled) {
                    s += ps[0].w * f(0);

                    foreach ((ddouble x, ddouble w, _) in ps.Skip(1)) {
                        s += w * (f(x) + f(-x));
                    }
                }
                else {
                    s += ps[0].wexp * f(0);

                    foreach ((ddouble x, _, ddouble w) in ps.Skip(1)) {
                        s += w * (f(x) + f(-x));
                    }
                }
            }

            return s;
        }
    }
}
