using DoubleDouble;
using System;
using System.Collections.ObjectModel;

namespace DoubleDoubleIntegrate {
    public enum GaussKronrodOrder : int {
        G3K7 = 3,
        G4K9 = 4,
        G7K15 = 7,
        G8K17 = 8,
        G15K31 = 15,
        G16K33 = 16,
        G31K63 = 31,
        G32K65 = 32
    }

    public static class GaussKronrodIntegral {

        public static (ddouble value, ddouble error) Integrate(Func<ddouble, ddouble> f, ddouble a, ddouble b, GaussKronrodOrder order = GaussKronrodOrder.G7K15) {
            ReadOnlyCollection<(ddouble x, ddouble wg, ddouble wk)> ps = GaussKronrodPoints.Table[order];

            ddouble sg = ddouble.Zero, sk = ddouble.Zero;
            ddouble r = b - a;

            if (!ddouble.IsFinite(r)) {
                throw new ArgumentException("Invalid integation interval.", $"{nameof(a)},{nameof(b)}");
            }

            for (int i = 0; i < ps.Count; i++) {
                ddouble x = ps[i].x;
                ddouble x_shifted = x * r + a;

                ddouble y = f(x_shifted);

                sk += ps[i].wk * y;

                if ((i & 1) == 1) {
                    sg += ps[i].wg * y;
                }
            }

            sk *= r;
            sg *= r;

            ddouble error = ddouble.Abs(sk - sg);

            return (sk, error);
        }

        private static (ddouble value, ddouble error) AdaptiveIntegrateFiniteInterval(Func<ddouble, ddouble> f, ddouble a, ddouble b, ddouble eps, GaussKronrodOrder order = GaussKronrodOrder.G7K15, int depth = 8) {
            (ddouble value, ddouble error) = Integrate(f, a, b, order);

            if (!(error > eps) || depth <= 0) {
                return (value, error);
            }

            ddouble c = (a + b) / 2, eps_half = eps / 2;

            (ddouble value1, ddouble error1) = AdaptiveIntegrateFiniteInterval(f, a, c, eps_half, order, depth - 1);
            (ddouble value2, ddouble error2) = AdaptiveIntegrateFiniteInterval(f, c, b, eps_half, order, depth - 1);

            return (value1 + value2, error1 + error2);
        }

        private static (ddouble value, ddouble error) AdaptiveIntegrateInfiniteInterval(Func<ddouble, ddouble> f, ddouble a, ddouble b, ddouble eps, GaussKronrodOrder order = GaussKronrodOrder.G7K15, int depth = 8) {
            if (ddouble.IsNaN(a) || ddouble.IsNaN(b)) {
                throw new ArgumentException("Invalid integation interval.", $"{nameof(a)},{nameof(b)}");
            }

            if (a > b) {
                (ddouble value, ddouble error) = AdaptiveIntegrateInfiniteInterval(f, b, a, eps, order, depth);

                return (-value, error);
            }

            if (ddouble.IsInfinity(a) && ddouble.IsInfinity(b)) {
                if (a.Sign == b.Sign) {
                    throw new ArgumentException("Invalid integation interval.", $"{nameof(a)},{nameof(b)}");
                }

                ddouble g(ddouble t) {
                    if (ddouble.IsZero(t)) {
                        return ddouble.Zero;
                    }

                    ddouble u = (1 - t) / t;

                    return (f(u) + f(-u)) / (t * t);
                }

                return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth);
            }

            if (ddouble.IsFinite(a) && ddouble.IsInfinity(b)) {
                if (ddouble.IsZero(a)) {
                    ddouble g(ddouble t) {
                        if (ddouble.IsZero(t)) {
                            return ddouble.Zero;
                        }

                        ddouble u = (1 - t) / t;

                        return f(u) / (t * t);
                    }

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth);
                }
                else {
                    ddouble g(ddouble t) {
                        if (ddouble.IsZero(t)) {
                            return ddouble.Zero;
                        }

                        ddouble u = (1 - t) / t;

                        return f(a + u) / (t * t);
                    }

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth);
                }
            }

            if (ddouble.IsInfinity(a) && ddouble.IsFinite(b)) {
                if (ddouble.IsZero(b)) {
                    ddouble g(ddouble t) {
                        if (ddouble.IsZero(t)) {
                            return ddouble.Zero;
                        }

                        ddouble u = (t - 1) / t;

                        return f(u) / (t * t);
                    }

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth);
                }
                else {
                    ddouble g(ddouble t) {
                        if (ddouble.IsZero(t)) {
                            return ddouble.Zero;
                        }

                        ddouble u = (t - 1) / t;

                        return f(b + u) / (t * t);
                    }

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth);
                }
            }

            throw new ArgumentException("Invalid integation interval.", $"{nameof(a)},{nameof(b)}");
        }

        public static (ddouble value, ddouble error) AdaptiveIntegrate(Func<ddouble, ddouble> f, ddouble a, ddouble b, ddouble eps, GaussKronrodOrder order = GaussKronrodOrder.G7K15, int depth = 8) {
            if (ddouble.IsFinite(a) && ddouble.IsFinite(b)) {
                return AdaptiveIntegrateFiniteInterval(f, a, b, eps, order, depth);
            }
            else {
                return AdaptiveIntegrateInfiniteInterval(f, a, b, eps, order, depth);
            }
        }
    }
}
