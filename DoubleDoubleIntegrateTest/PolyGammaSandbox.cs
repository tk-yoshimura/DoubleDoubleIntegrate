using DoubleDouble;
using DoubleDoubleIntegrate;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace DoubleDoubleIntegrateTest {
    [TestClass]
    public class PolyGammaSandbox {
        [TestMethod]
        public void PolyGammaIntegrate() {
            int n = 32;
            ddouble z = 64;

            for (int pts = 16; pts <= 64; pts++) {
                Console.WriteLine(pts);

                ddouble y = PolyGammaNearZero(n, z, pts, Math.Max(1, 8d / n));
                ddouble y_limit = PolyGammaLimit(n, z);

                Console.WriteLine(y);
                Console.WriteLine(y_limit);
            }
        }

        public ddouble PolyGammaNearZero(int n, ddouble x, int pts, ddouble scale) {
            if (n < 0) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }
            if (n == 0) {
                return ddouble.Digamma(x);
            }

            if (x < 0) { 
                throw new ArgumentOutOfRangeException(nameof(x));
            }

            if (x < 1) {
                ddouble v = PolyGammaNearZero(n, x + 1, pts, scale);
                ddouble y = v + ddouble.Gamma(n + 1) / ddouble.Pow(x, n + 1);

                return y;
            }

            ddouble ir, it;

            ddouble polygamma_ir(ddouble t){
                ddouble y = ddouble.Pow(t, n) * ddouble.Exp(-x * t) / (1 - ddouble.Exp(-t));

                return y;
            }; 

            ddouble polygamma_it(ddouble u){
                ddouble v = (u + scale * n) / x;
                ddouble y = ddouble.Pow(v, n) / (1 - ddouble.Exp(-v));

                return y;
            };

            ir = GaussLegendreIntegral.Integrate(polygamma_ir, ddouble.Zero, scale * n / x, pts);
            it = ddouble.Exp(-scale * n) / x * GaussLaguerreIntegral.Integrate(polygamma_it, pts, f_expscaled: true);

            ddouble i = ir + it;

            Console.WriteLine(ir);
            Console.WriteLine(it);

            Console.WriteLine(i);

            return i;
        }

        public ddouble PolyGammaLimit(int n, ddouble x) {

            ddouble inv_x2 = 1 / (x * x), c = ddouble.Pow(x, -n);
            ddouble v = c * ddouble.Gamma(n) * (2 * x + n) / (2 * x);
            ddouble u = c * ddouble.Gamma(n + 2) / 2 * inv_x2;
            ddouble dv = ddouble.BernoulliSequence[1] * u;
 
            v += dv;

            for (int k = 2; k < ddouble.BernoulliSequence.Count; k++) {
                u *= inv_x2 * checked((n + 2 * k - 2) * (n + 2 * k - 1)) / checked((2 * k) * (2 * k - 1));
                dv = ddouble.BernoulliSequence[k] * u;
                ddouble next_v = v + dv;

                if (v == next_v) {
                    break;
                }
                if (ddouble.IsNaN(next_v)) {
                    return ddouble.Zero;
                }

                v = next_v;
            }

            return v;
        }
    }
}
