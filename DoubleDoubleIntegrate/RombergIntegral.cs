using DoubleDouble;
using System;

namespace DoubleDoubleIntegrate {

    public static class RombergIntegral {

        public static ddouble Integrate(Func<ddouble, ddouble> f, ddouble a, ddouble b, int precision_level = 10) {
            ddouble h = b - a;

            if (!ddouble.IsFinite(h)) {
                throw new ArgumentOutOfRangeException($"{nameof(a)},{nameof(b)}");
            }

            if (precision_level < 1 || precision_level > 24) {
                throw new ArgumentOutOfRangeException(nameof(precision_level));
            }

            int max_div = 1 << precision_level;
            ddouble min_h = h / max_div;
            ddouble[] v = new ddouble[max_div + 1];
            RichardsonExtrapolation conv = new();

            for (int i = 0; i <= max_div; i++) {
                v[i] = f(a + i * min_h);
            }

            ddouble t = h * (v[0] + v[max_div]) / 2, new_t;
            conv.Inject(t);

            for (int s = max_div; s > 1; s /= 2) {
                new_t = 0;
                for (int i = s / 2; i < max_div; i += s) {
                    new_t += v[i];
                }

                h /= 2;
                t = t / 2 + h * new_t;

                conv.Inject(t);
            }

            return conv.ConvergenceValue;
        }
    }
}
