using DoubleDouble;
using System;
using System.Collections.Generic;

namespace DoubleDoubleIntegrate {

    internal class RichardsonExtrapolation {
        static readonly List<ddouble> rs = new() { ddouble.NaN };
        readonly List<ddouble[]> values = new();

        public void Inject(ddouble new_value) {
            if (SeriesCount <= 0) {
                values.Add(new ddouble[] { new_value });
                return;
            }

            ddouble[] t = values[SeriesCount - 1], t_next = new ddouble[SeriesCount + 1];

            t_next[0] = new_value;

            for (int i = 1; i <= SeriesCount; i++) {
                t_next[i] = t_next[i - 1] + (t_next[i - 1] - t[i - 1]) * R(i);
            }

            values.Add(t_next);
        }

        private static ddouble R(int i) {
            for (int k = rs.Count; k <= i; k++) {
                ddouble r = 1d / (ddouble.Ldexp(1d, k * 2) - 1);

                rs.Add(r);
            }

            return rs[i];
        }

        public IEnumerable<ddouble> Series {
            get {
                for (int i = 0; i < values.Count; i++) {
                    yield return values[i][i];
                }
            }
        }

        public ddouble ConvergenceValue {
            get {
                if (SeriesCount <= 0) {
                    throw new InvalidOperationException();
                }

                return values[SeriesCount - 1][SeriesCount - 1];
            }
        }

        public int SeriesCount => values.Count;

        public ddouble Epsilon {
            get {
                if (SeriesCount <= 1) {
                    throw new InvalidOperationException();
                }

                return ddouble.Abs(values[SeriesCount - 1][SeriesCount - 1] - values[SeriesCount - 2][SeriesCount - 2]);
            }
        }
    }
}
