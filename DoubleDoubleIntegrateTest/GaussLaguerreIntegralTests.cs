using DoubleDouble;
using DoubleDoubleIntegrate;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace DoubleDoubleIntegrateTest {
    [TestClass]
    public class GaussLaguerreIntegralTests {
        [TestMethod]
        public void IntegrateExpTest() {
            for (int n = 4; n <= 64; n++) {
                ddouble y = GaussLaguerreIntegral.Integrate((x) => ddouble.Exp(-x), n, f_expscaled: false);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(1 - GaussLaguerreIntegral.Integrate((x) => ddouble.Exp(-x), 32, f_expscaled: false)), 1e-29);

            for (int n = 4; n <= 64; n++) {
                ddouble y = GaussLaguerreIntegral.Integrate((x) => 1, n, f_expscaled: true);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(1 - GaussLaguerreIntegral.Integrate((x) => 1, 32, f_expscaled: true)), 1e-29);
        }

        [TestMethod]
        public void IntegrateXExpTest() {
            for (int n = 4; n <= 64; n++) {
                ddouble y = GaussLaguerreIntegral.Integrate((x) => ddouble.Exp(-x) * x * x, n, f_expscaled: false);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(2 - GaussLaguerreIntegral.Integrate((x) => ddouble.Exp(-x) * x * x, 32, f_expscaled: false)), 1e-29);

            for (int n = 4; n <= 64; n++) {
                ddouble y = GaussLaguerreIntegral.Integrate((x) => ddouble.Exp(-x / 2) * x * x, n, f_expscaled: false);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(16 - GaussLaguerreIntegral.Integrate((x) => ddouble.Exp(-x / 2) * x * x, 40, f_expscaled: false)), 1e-29);

            for (int n = 4; n <= 64; n++) {
                ddouble y = GaussLaguerreIntegral.Integrate((x) => x * x, n, f_expscaled: true);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(2 - GaussLaguerreIntegral.Integrate((x) => x * x, 32, f_expscaled: true)), 1e-29);
        }
    }
}
