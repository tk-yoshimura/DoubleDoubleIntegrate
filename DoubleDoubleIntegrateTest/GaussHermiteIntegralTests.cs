using DoubleDouble;
using DoubleDoubleIntegrate;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace DoubleDoubleIntegrateTest {
    [TestClass]
    public class GaussHermiteIntegralTests {
        [TestMethod]
        public void IntegrateExpTest() {
            for (int n = 4; n <= 64; n++) {
                ddouble y = GaussHermiteIntegral.Integrate((x) => ddouble.Exp(-x * x), n, f_expscaled: false);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(ddouble.Sqrt(ddouble.Pi) - GaussHermiteIntegral.Integrate((x) => ddouble.Exp(-x * x), 32, f_expscaled: false)), 1e-29);
            Assert.AreEqual(0d, (double)(ddouble.Sqrt(ddouble.Pi) - GaussHermiteIntegral.Integrate((x) => ddouble.Exp(-x * x), 33, f_expscaled: false)), 1e-29);

            for (int n = 4; n <= 64; n++) {
                ddouble y = GaussHermiteIntegral.Integrate((x) => 1, n, f_expscaled: true);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(ddouble.Sqrt(ddouble.Pi) - GaussHermiteIntegral.Integrate((x) => 1, 32, f_expscaled: true)), 1e-29);
            Assert.AreEqual(0d, (double)(ddouble.Sqrt(ddouble.Pi) - GaussHermiteIntegral.Integrate((x) => 1, 33, f_expscaled: true)), 1e-29);
        }

        [TestMethod]
        public void IntegrateXExpTest() {
            for (int n = 4; n <= 64; n++) {
                ddouble y = GaussHermiteIntegral.Integrate((x) => ddouble.Exp(-x * x) * x * x, n, f_expscaled: false);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(ddouble.Sqrt(ddouble.Pi) / 2 - GaussHermiteIntegral.Integrate((x) => ddouble.Exp(-x * x) * x * x, 32, f_expscaled: false)), 1e-29);
            Assert.AreEqual(0d, (double)(ddouble.Sqrt(ddouble.Pi) / 2 - GaussHermiteIntegral.Integrate((x) => ddouble.Exp(-x * x) * x * x, 33, f_expscaled: false)), 1e-29);

            for (int n = 4; n <= 64; n++) {
                ddouble y = GaussHermiteIntegral.Integrate((x) => ddouble.Exp(-x * x / 2) * x * x, n, f_expscaled: false);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(ddouble.Sqrt(ddouble.Pi * 2) - GaussHermiteIntegral.Integrate((x) => ddouble.Exp(-x * x / 2) * x * x, 63, f_expscaled: false)), 1e-25);
            Assert.AreEqual(0d, (double)(ddouble.Sqrt(ddouble.Pi * 2) - GaussHermiteIntegral.Integrate((x) => ddouble.Exp(-x * x / 2) * x * x, 64, f_expscaled: false)), 1e-25);

            for (int n = 4; n <= 64; n++) {
                ddouble y = GaussHermiteIntegral.Integrate((x) => x * x, n, f_expscaled: true);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(ddouble.Sqrt(ddouble.Pi) / 2 - GaussHermiteIntegral.Integrate((x) => x * x, 32, f_expscaled: true)), 1e-29);
            Assert.AreEqual(0d, (double)(ddouble.Sqrt(ddouble.Pi) / 2 - GaussHermiteIntegral.Integrate((x) => x * x, 33, f_expscaled: true)), 1e-29);
        }
    }
}
