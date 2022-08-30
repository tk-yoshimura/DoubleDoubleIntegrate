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
            ddouble z = 32;

            for (int pts = 16; pts <= 64; pts++) {
                Console.WriteLine(pts);

                ddouble y = PolyGamma(n, z, pts, Math.Max(1, 8d / n));

                Console.WriteLine(y);
            }
        }

        public ddouble PolyGamma(int n, ddouble x, int pts, ddouble scale) {
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
                ddouble v = PolyGamma(n, x + 1, pts, scale);
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
    }
}
