using DoubleDouble;
using DoubleDoubleIntegrate;

namespace DoubleDoubleIntegrateSandbox {
    internal class Program {
        static void Main() {
            //PolyGammaNearZero(1, 1d / 3, 20, Math.Max(1, 8d / 1));

            //{
            //    ddouble z = 32;
            //
            //    for (int n = 2; n <= 16; n++) {
            //        Console.WriteLine($"n={n}");
            //
            //        ddouble y_limit = PolyGammaLimit(n, z);
            //        Console.WriteLine($"{y_limit:e16}");
            //
            //        for (int pts = 16; pts <= 24; pts++) {
            //            ddouble y = PolyGammaNearZero(n, z, pts, Math.Max(1, 8d / n));
            //
            //            Console.WriteLine($"pts={pts:D2}, {y:e16}");
            //        }
            //    }
            //}

            for (ddouble z = 32; z <= 32; z *= 2) {
                Console.WriteLine($"z={z}");

                for (int n = 1; n <= 16; n++) {
                    Console.WriteLine($"n={n}");

                    for (int pts = 16; pts <= 24; pts++) {
                        ddouble y = PolyGammaNearZero(n, z, pts, Math.Max(1, 8d / n));

                        Console.WriteLine($"pts={pts:D2}, {y:e20}");
                    }
                }
            }

            {
                System.Collections.ObjectModel.ReadOnlyCollection<(ddouble x, ddouble w)> points = GaussLegendrePoints.Table[36];

                foreach ((ddouble x, ddouble w) in points) {
                    Console.WriteLine($"({x:e20}, {w:e20}), ");
                }
            }

            Console.WriteLine("");

            {
                System.Collections.ObjectModel.ReadOnlyCollection<(ddouble x, ddouble w, ddouble wexp)> points = GaussLaguerrePoints.Table[36];

                foreach ((ddouble x, ddouble w, _) in points) {
                    Console.WriteLine($"({x:e20}, {w:e20}), ");
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }

        public static ddouble PolyGammaNearZero(int n, ddouble x, int pts, ddouble scale) {
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

            ddouble polygamma_ir(ddouble t) {
                ddouble y = ddouble.Pow(t, n) * ddouble.Exp(-x * t) / (1 - ddouble.Exp(-t));
                return y;
            };
            ddouble polygamma_it(ddouble u) {
                ddouble v = (u + scale * n) / x;
                ddouble y = ddouble.Pow(v, n) / (1 - ddouble.Exp(-v));
                return y;
            };

            ir = GaussLegendreIntegral.Integrate(polygamma_ir, ddouble.Zero, scale * n / x, pts);
            it = ddouble.Exp(-scale * n) / x * GaussLaguerreIntegral.Integrate(polygamma_it, pts, f_expscaled: true);

            ddouble i = ir + it;

            //Console.WriteLine(ir);
            //Console.WriteLine(it);
            //
            //Console.WriteLine(i);

            return i;
        }

        public static ddouble PolyGammaLimit(int n, ddouble x, int max_terms = 10, double etoi = 2e-16) {
            ddouble inv_x2 = 1 / (x * x), c = ddouble.Pow(x, -n);
            ddouble v = c * ddouble.Gamma(n) * (2 * x + n) / (2 * x);
            ddouble u = c * ddouble.Gamma(n + 2) / 2 * inv_x2;
            ddouble dv = ddouble.BernoulliSequence[1] * u;
            v += dv;

            for (int k = 2; k <= max_terms; k++) {
                u *= inv_x2 * checked((n + 2 * k - 2) * (n + 2 * k - 1)) / checked((2 * k) * (2 * k - 1));
                dv = ddouble.BernoulliSequence[k] * u;
                ddouble next_v = v + dv;

                Console.WriteLine(dv);

                if (ddouble.Abs(dv) < ddouble.Abs(v) * etoi) {
                    return v;
                }
                if (ddouble.IsNaN(next_v)) {
                    break;
                }

                v = next_v;
            }

            return ddouble.NaN;
        }
    }
}