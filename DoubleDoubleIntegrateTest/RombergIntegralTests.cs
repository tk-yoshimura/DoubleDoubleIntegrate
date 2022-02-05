using DoubleDouble;
using DoubleDoubleIntegrate;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace DoubleDoubleIntegrateTest {
    [TestClass()]
    public class RombergIntegralTests {
        [TestMethod()]
        public void IntegrateTest() {
            static ddouble f(ddouble x) => ddouble.Sqrt(1 - x * x);

            for (int level = 1; level <= 16; level++) {
                ddouble v = RombergIntegral.Integrate(f, 0, ddouble.Sqrt(2) / 2, level);

                Console.WriteLine($"{level}\t{v}");
            }

            {
                ddouble v = RombergIntegral.Integrate(f, 0, ddouble.Sqrt(2) / 2, 20);
                Assert.AreEqual(0, (double)((ddouble.PI + 2) / 8 - v), 1e-20);
            }
        }
    }
}