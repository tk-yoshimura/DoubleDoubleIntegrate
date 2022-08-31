using DoubleDouble;
using DoubleDoubleIntegrate;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace DoubleDoubleIntegrateTest {
    [TestClass]
    public class PolyLogSandbox {
        [TestMethod]
        public void PolyLogIntegrate() {
            int n = 3;
            ddouble z = 2;

            for (int pts = 16; pts <= 64; pts++) {
                Console.WriteLine(pts);

                ddouble y = PolyLogIntegral(n, z, pts);

                Console.WriteLine(y);
            }
        }

        public ddouble PolyLogIntegral(int n, ddouble x, int pts) {
            if (n <= 1) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }
            if (x < 1) { 
                throw new ArgumentOutOfRangeException(nameof(x));
            }

            ddouble ir, it, h = PolylogIntegrandPeak(n, (double)x) * 6;

            ddouble polygamma_ir(ddouble t){
                ddouble y = ddouble.Pow(t, n - 1) / (ddouble.Exp(t) / x + 1);

                return y;
            }; 

            ddouble polygamma_it(ddouble u){
                ddouble v = u + h;
                ddouble y = ddouble.Pow(v, n - 1) / (1 / x + ddouble.Exp(-v));

                return y;
            };

            ir = GaussLegendreIntegral.Integrate(polygamma_ir, ddouble.Zero, h, pts);
            it = ddouble.Exp(-h) * GaussLaguerreIntegral.Integrate(polygamma_it, pts, f_expscaled: true);

            ddouble i = ir + it;

            Console.WriteLine(ir);
            Console.WriteLine(it);

            Console.WriteLine(i);

            return i;
        }

        static double PolylogIntegrandPeak(int n, double x) {
            if (n <= 0) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }
            if (n == 1) {
                return 0;
            }

            double t = 0.375 + 0.625 * Math.Log2(x + 4);

            for (int i = 0; i < 8; i++) {
                double xexp = Math.Exp(-t) * x;
                double d = (n - 1) * (xexp + 1) - t;
                double dv = (n - 1) * xexp + 1;
                double dt = d / dv;

                t += dt;

                if (Math.Abs(dt / t) < 1e-15) {
                    break;
                }
            }

            return t;
        }
    }
}
