using DoubleDouble;
using DoubleDoubleIntegrate;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace DoubleDoubleIntegrateTest {
    [TestClass]
    public class PolyLogSandbox {
        [TestMethod]
        public void PolyLogIntegrate() {
            int n = 2;
            ddouble z = 1024;

            for (int pts = 16; pts <= 64; pts++) {
                Console.WriteLine(pts);

                ddouble y = PolyLogIntegralR2(n, z, pts);

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

            ddouble h0 = PolylogIntegrandPeak(n, (double)x), h1 = h0 * 3;

            ddouble polygamma_ir(ddouble t){
                ddouble y = ddouble.Pow(t, n - 1) / (ddouble.Exp(t) / x + 1);

                return y;
            }; 

            ddouble polygamma_it(ddouble u){
                ddouble v = u + h1;
                ddouble y = ddouble.Pow(v, n - 1) / (1 / x + ddouble.Exp(-v));

                return y;
            };

            ddouble ir0 = GaussLegendreIntegral.Integrate(polygamma_ir, ddouble.Zero, h0, pts);
            ddouble ir1 = GaussLegendreIntegral.Integrate(polygamma_ir, h0, h1, pts);
            ddouble it = ddouble.Exp(-h1) * GaussLaguerreIntegral.Integrate(polygamma_it, pts, f_expscaled: true);

            ddouble i = ir0 + ir1 + it;

            Console.WriteLine(ir0);
            Console.WriteLine(ir1);
            Console.WriteLine(it);

            Console.WriteLine(i);

            return i;
        }

        public ddouble PolyLogIntegralR2(int n, ddouble x, int pts) {
            if (n <= 1) {
                throw new ArgumentOutOfRangeException(nameof(n));
            }
            if (x < 0.5) { 
                throw new ArgumentOutOfRangeException(nameof(x));
            }

            ddouble h = PolylogIntegrandPeak(n, (double)x);

            ddouble polygamma_ir(ddouble t){
                ddouble y = ddouble.Pow(t, n - 1) / (ddouble.Exp(t) / x + 1);

                return y;
            }; 


            int k;
            ddouble ir = 0, ir0, ir1;
            ddouble thr = 1.25 * ddouble.Exp(-h);

            ir0 = GaussLegendreIntegral.Integrate(polygamma_ir, ddouble.Zero, h, pts);
            ir += ir0;

            Console.WriteLine($"0h: {ir0}");

            for (k = 1; ; k++) { 
                ir1 = GaussLegendreIntegral.Integrate(polygamma_ir, k * h, (k + 1) * h, pts);
                ir += ir1;

                Console.WriteLine($"{k}h: {ir1}");

                if (!(ir1 / ir0 >= thr) || k >= 32) {
                    break;
                }

                ir0 = ir1;
            }

            ddouble sh = (k + 1) * h;

            ddouble polygamma_it(ddouble u){
                ddouble v = u + sh;
                ddouble y = ddouble.Pow(v, n - 1) / (1 / x + ddouble.Exp(-v));

                return y;
            };
            
            ddouble it = ddouble.Exp(-sh) * GaussLaguerreIntegral.Integrate(polygamma_it, pts, f_expscaled: true);

            ddouble i = ir + it;

            Console.WriteLine($"{k + 1}h+: {it}");

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
