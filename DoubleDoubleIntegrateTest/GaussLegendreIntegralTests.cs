using DoubleDouble;
using DoubleDoubleIntegrate;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace DoubleDoubleIntegrateTest {
    [TestClass]
    public class GaussLegendreIntegralTests {
        [TestMethod]
        public void IntegrateSinTest() {
            for (int n = 4; n <= 64; n++) {
                ddouble y = GaussLegendreIntegral.Integrate(ddouble.Sin, ddouble.Zero, ddouble.PI, n);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(2 - GaussLegendreIntegral.Integrate(ddouble.Sin, ddouble.Zero, ddouble.PI, 32)), 1e-30);
        }

        [TestMethod]
        public void IntegrateExpTest() {
            for (int n = 4; n <= 64; n++) {
                ddouble y = GaussLegendreIntegral.Integrate(ddouble.Exp, 1, 4, n);

                Console.WriteLine($"{n}\t {y}");
            }

            Assert.AreEqual(0d, (double)(ddouble.E * (ddouble.Cube(ddouble.E) - 1) - GaussLegendreIntegral.Integrate(ddouble.Exp, 1, 4, 32)), 1e-29);
        }
    }
}
