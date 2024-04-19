using DoubleDouble;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;

namespace DoubleDoubleIntegrate {
    public enum GaussKronrodOrder : int {
        G3K7 = 3,
        G4K9 = 4,
        G7K15 = 7,
        G8K17 = 8,
        G15K31 = 15,
        G16K33 = 16,
        G31K63 = 31,
        G32K65 = 32,
        G63K127 = 63,
        G64K129 = 64,
        G127K255 = 127,
        G128K257 = 128,
    }

    public static class GaussKronrodIntegral {

        public static (ddouble value, ddouble error) Integrate(Func<ddouble, ddouble> f, ddouble a, ddouble b, GaussKronrodOrder order = GaussKronrodOrder.G31K63) {
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

        private static (ddouble value, ddouble error, long eval_points) UnlimitedIntegrateFiniteInterval(Func<ddouble, ddouble> f, ddouble a, ddouble b, ddouble eps, GaussKronrodOrder order) {
            Stack<(ddouble a, ddouble b, ddouble eps)> stack = new();
            stack.Push((a, b, eps));

            long eval_points_sum = 0;
            ddouble value_sum = 0d, error_sum = 0d;

            while (stack.Count > 0) {
                (a, b, eps) = stack.Pop();

                (ddouble value, ddouble error) = Integrate(f, a, b, order);

                long eval_points = 1 + 2 * (int)order;
                eval_points_sum += eval_points;

                if (!(error > eps)) {
                    value_sum += value;
                    error_sum += error;
                    continue;
                }

                ddouble c = (a + b) / 2, eps_half = eps / 2;
                stack.Push((a, c, eps_half));
                stack.Push((c, b, eps_half));
            }

            return (value_sum, error_sum, eval_points_sum);
        }

        private static (ddouble value, ddouble error, long eval_points) LimitedDepthIntegrateFiniteInterval(Func<ddouble, ddouble> f, ddouble a, ddouble b, ddouble eps, GaussKronrodOrder order, int maxdepth) {
            Debug.Assert(maxdepth >= 0);

            Stack<(ddouble a, ddouble b, ddouble eps, int depth)> stack = new();
            stack.Push((a, b, eps, maxdepth));

            long eval_points_sum = 0;
            ddouble value_sum = 0d, error_sum = 0d;

            while (stack.Count > 0) {
                (a, b, eps, int depth) = stack.Pop();

                (ddouble value, ddouble error) = Integrate(f, a, b, order);

                long eval_points = 1 + 2 * (int)order;
                eval_points_sum += eval_points;

                if (!(error > eps) || depth <= 0) {
                    value_sum += value;
                    error_sum += error;
                    continue;
                }

                ddouble c = (a + b) / 2, eps_half = eps / 2;
                depth -= 1;
                stack.Push((a, c, eps_half, depth));
                stack.Push((c, b, eps_half, depth));
            }

            return (value_sum, error_sum, eval_points_sum);
        }

        private static (ddouble value, ddouble error, long eval_points) LimitedEvalIntegrateFiniteInterval(Func<ddouble, ddouble> f, ddouble a, ddouble b, ddouble eps, GaussKronrodOrder order, long discontinue_eval_points) {
            Debug.Assert(discontinue_eval_points >= 0);

            PriorityQueue<(ddouble a, ddouble b, ddouble eps), long> queue = new();
            queue.Enqueue((a, b, eps), 0);

            long eval_points_sum = 0;
            ddouble value_sum = 0d, error_sum = 0d;

            while (queue.Count > 0) {
                (a, b, eps) = queue.Dequeue();

                (ddouble value, ddouble error) = Integrate(f, a, b, order);

                long eval_points = 1 + 2 * (int)order;
                eval_points_sum += eval_points;

                if (!(error > eps) || eval_points_sum > discontinue_eval_points) {
                    value_sum += value;
                    error_sum += error;
                    continue;
                }

                ddouble c = (a + b) / 2, eps_half = eps / 2;
                long priority = double.ILogB((double)error);
                queue.Enqueue((a, c, eps_half), -priority);
                queue.Enqueue((c, b, eps_half), -priority);
            }

            return (value_sum, error_sum, eval_points_sum);
        }

        private static (ddouble value, ddouble error, long eval_points) LimitedDepthAndEvalIntegrateFiniteInterval(Func<ddouble, ddouble> f, ddouble a, ddouble b, ddouble eps, GaussKronrodOrder order, int maxdepth, long discontinue_eval_points) {
            Debug.Assert(maxdepth >= 0);
            Debug.Assert(discontinue_eval_points >= 0);

            PriorityQueue<(ddouble a, ddouble b, ddouble eps, int depth), long> queue = new();
            queue.Enqueue((a, b, eps, maxdepth), 0);

            long eval_points_sum = 0;
            ddouble value_sum = 0d, error_sum = 0d;

            while (queue.Count > 0) {
                (a, b, eps, int depth) = queue.Dequeue();

                (ddouble value, ddouble error) = Integrate(f, a, b, order);

                long eval_points = 1 + 2 * (int)order;
                eval_points_sum += eval_points;

                if (!(error > eps) || depth <= 0 || eval_points_sum > discontinue_eval_points) {
                    value_sum += value;
                    error_sum += error;
                    continue;
                }

                ddouble c = (a + b) / 2, eps_half = eps / 2;
                long priority = double.ILogB((double)error);
                depth -= 1;
                queue.Enqueue((a, c, eps_half, depth), -priority);
                queue.Enqueue((c, b, eps_half, depth), -priority);
            }

            return (value_sum, error_sum, eval_points_sum);
        }

        private static (ddouble value, ddouble error, long eval_points) AdaptiveIntegrateFiniteInterval(Func<ddouble, ddouble> f, ddouble a, ddouble b, ddouble eps, GaussKronrodOrder order, int maxdepth, long discontinue_eval_points) {
            if (maxdepth >= 0 && discontinue_eval_points >= 0) {
                return LimitedDepthAndEvalIntegrateFiniteInterval(f, a, b, eps, order, maxdepth, discontinue_eval_points);
            }
            if (maxdepth >= 0) {
                return LimitedDepthIntegrateFiniteInterval(f, a, b, eps, order, maxdepth);
            }
            if (discontinue_eval_points >= 0) {
                return LimitedEvalIntegrateFiniteInterval(f, a, b, eps, order, discontinue_eval_points);
            }

            return UnlimitedIntegrateFiniteInterval(f, a, b, eps, order);
        }

        private static (ddouble value, ddouble error, long eval_points) AdaptiveIntegrateInfiniteInterval(Func<ddouble, ddouble> f, ddouble a, ddouble b, ddouble eps, GaussKronrodOrder order, int depth, long discontinue_eval_points) {
            if (ddouble.IsNaN(a) || ddouble.IsNaN(b)) {
                throw new ArgumentException("Invalid integation interval.", $"{nameof(a)},{nameof(b)}");
            }

            if (a > b) {
                (ddouble value, ddouble error, long eval_points) = AdaptiveIntegrateInfiniteInterval(f, b, a, eps, order, depth, discontinue_eval_points);

                return (-value, error, eval_points);
            }

            if (ddouble.IsInfinity(a) && ddouble.IsInfinity(b)) {
                if (ddouble.Sign(a) == ddouble.Sign(b)) {
                    throw new ArgumentException("Invalid integation interval.", $"{nameof(a)},{nameof(b)}");
                }

                ddouble g(ddouble t) {
                    if (ddouble.IsZero(t)) {
                        return ddouble.Zero;
                    }

                    ddouble u = (1 - t) / t;

                    return (f(u) + f(-u)) / (t * t);
                }

                return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth, discontinue_eval_points);
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

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth, discontinue_eval_points);
                }
                else {
                    ddouble g(ddouble t) {
                        if (ddouble.IsZero(t)) {
                            return ddouble.Zero;
                        }

                        ddouble u = (1 - t) / t;

                        return f(a + u) / (t * t);
                    }

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth, discontinue_eval_points);
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

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth, discontinue_eval_points);
                }
                else {
                    ddouble g(ddouble t) {
                        if (ddouble.IsZero(t)) {
                            return ddouble.Zero;
                        }

                        ddouble u = (t - 1) / t;

                        return f(b + u) / (t * t);
                    }

                    return AdaptiveIntegrateFiniteInterval(g, 0, 1, eps, order, depth, discontinue_eval_points);
                }
            }

            throw new ArgumentException("Invalid integation interval.", $"{nameof(a)},{nameof(b)}");
        }

        public static (ddouble value, ddouble error, long eval_points) AdaptiveIntegrate(Func<ddouble, ddouble> f, ddouble a, ddouble b, ddouble eps, GaussKronrodOrder order = GaussKronrodOrder.G31K63, int maxdepth = -1, long discontinue_eval_points = -1) {
            if (maxdepth < -1) {
                throw new ArgumentOutOfRangeException(nameof(maxdepth), "Invalid param. maxdepth=-1: infinite, maxdepth>=0: finite");
            }
            if (discontinue_eval_points < -1) {
                throw new ArgumentOutOfRangeException(nameof(discontinue_eval_points), "Invalid param. discontinue_eval_points=-1: infinite, discontinue_eval_points>=0: finite");
            }
            if (!(eps >= 0d)) {
                throw new ArgumentOutOfRangeException(nameof(eps), "Invalid param. eps must be nonnegative value");
            }

            if (ddouble.IsFinite(a) && ddouble.IsFinite(b)) {
                return AdaptiveIntegrateFiniteInterval(f, a, b, eps, order, maxdepth, discontinue_eval_points);
            }
            else {
                return AdaptiveIntegrateInfiniteInterval(f, a, b, eps, order, maxdepth, discontinue_eval_points);
            }
        }
    }
}
