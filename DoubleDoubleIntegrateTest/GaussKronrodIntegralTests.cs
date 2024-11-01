using DoubleDouble;
using DoubleDoubleIntegrate;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace DoubleDoubleIntegrateTest {
    [TestClass]
    public class GaussKronrodIntegralTests {
        [TestMethod]
        public void IntegrateSinTest() {
            foreach (GaussKronrodOrder order in Enum.GetValues<GaussKronrodOrder>()) {
                (ddouble y, ddouble err) = GaussKronrodIntegral.Integrate(ddouble.Sin, ddouble.Zero, ddouble.Pi, order);

                Console.WriteLine($"{order}\t {y}\t {err}");
            }

            Assert.AreEqual(0d, (double)(2 - GaussKronrodIntegral.Integrate(ddouble.Sin, ddouble.Zero, ddouble.Pi, GaussKronrodOrder.G32K65).value), 1e-30);
        }

        [TestMethod]
        public void IntegrateExpTest() {
            foreach (GaussKronrodOrder order in Enum.GetValues<GaussKronrodOrder>()) {
                (ddouble y, ddouble err) = GaussKronrodIntegral.Integrate(ddouble.Exp, 1, 4, order);

                Console.WriteLine($"{order}\t {y}\t {err}");
            }

            Assert.AreEqual(0d, (double)(ddouble.E * (ddouble.Cube(ddouble.E) - 1) - GaussKronrodIntegral.Integrate(ddouble.Exp, 1, 4, GaussKronrodOrder.G32K65).value), 1e-29);
        }

        [TestMethod]
        public void AdaptiveIntegrateExpTest() {
            foreach (GaussKronrodOrder order in Enum.GetValues<GaussKronrodOrder>()) {
                (ddouble y, ddouble err, long eval_points) = GaussKronrodIntegral.AdaptiveIntegrate(ddouble.Exp, 1, 4, 1e-25, order, maxdepth: 10);

                Console.WriteLine($"{order}\t {y}\t {err}\t {eval_points}");
            }

            Assert.AreEqual(0d, (double)(ddouble.E * (ddouble.Cube(ddouble.E) - 1) - GaussKronrodIntegral.AdaptiveIntegrate(ddouble.Exp, 1, 4, 1e-25, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-25);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest1() {
            foreach (GaussKronrodOrder order in Enum.GetValues<GaussKronrodOrder>()) {
                (ddouble y, ddouble err, long eval_points) = GaussKronrodIntegral.AdaptiveIntegrate(x => ddouble.Exp(-x), 0, ddouble.PositiveInfinity, 1e-25, order, maxdepth: 10);

                Console.WriteLine($"{order}\t {y}\t {err}\t {eval_points}");
            }

            Assert.AreEqual(0d, (double)(1 - GaussKronrodIntegral.AdaptiveIntegrate(x => ddouble.Exp(-x), 0, ddouble.PositiveInfinity, 1e-25, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-25);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest2() {
            foreach (GaussKronrodOrder order in Enum.GetValues<GaussKronrodOrder>()) {
                (ddouble y, ddouble err, long eval_points) = GaussKronrodIntegral.AdaptiveIntegrate(x => ddouble.Exp(-x), 1, ddouble.PositiveInfinity, 1e-25, order, maxdepth: 10);

                Console.WriteLine($"{order}\t {y}\t {err}\t {eval_points}");
            }

            ddouble expected = "1.353352832366126918939994949724844034076315459095758814681588724e-1";

            Assert.AreEqual(0d, (double)(expected - GaussKronrodIntegral.AdaptiveIntegrate(x => ddouble.Exp(-x), 2, ddouble.PositiveInfinity, 1e-25, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-25);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest3() {
            Assert.AreEqual(0d, (double)(-1 - GaussKronrodIntegral.AdaptiveIntegrate(x => ddouble.Exp(-x), ddouble.PositiveInfinity, 0, 1e-25, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-25);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest4() {
            ddouble expected = "-1.353352832366126918939994949724844034076315459095758814681588724e-1";

            Assert.AreEqual(0d, (double)(expected - GaussKronrodIntegral.AdaptiveIntegrate(x => ddouble.Exp(-x), ddouble.PositiveInfinity, 2, 1e-25, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-25);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest5() {
            Assert.AreEqual(0d, (double)(1 - GaussKronrodIntegral.AdaptiveIntegrate(ddouble.Exp, ddouble.NegativeInfinity, 0, 1e-25, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-25);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest6() {
            ddouble expected = ddouble.E * ddouble.E;

            Assert.AreEqual(0d, (double)(expected - GaussKronrodIntegral.AdaptiveIntegrate(ddouble.Exp, ddouble.NegativeInfinity, 2, 1e-25, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-25);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest7() {
            Assert.AreEqual(0d, (double)(-1 - GaussKronrodIntegral.AdaptiveIntegrate(ddouble.Exp, 0, ddouble.NegativeInfinity, 1e-25, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-25);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest8() {
            ddouble expected = -ddouble.E * ddouble.E;

            Assert.AreEqual(0d, (double)(expected - GaussKronrodIntegral.AdaptiveIntegrate(ddouble.Exp, 2, ddouble.NegativeInfinity, 1e-25, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-25);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest9() {
            ddouble expected = ddouble.Sqrt(ddouble.Pi);

            Assert.AreEqual(0d, (double)(expected - GaussKronrodIntegral.AdaptiveIntegrate(x => ddouble.Exp(-x * x), ddouble.NegativeInfinity, ddouble.PositiveInfinity, 1e-25, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-25);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest10() {
            ddouble expected = -ddouble.Sqrt(ddouble.Pi);

            Assert.AreEqual(0d, (double)(expected - GaussKronrodIntegral.AdaptiveIntegrate(x => ddouble.Exp(-x * x), ddouble.PositiveInfinity, ddouble.NegativeInfinity, 1e-25, GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-25);
        }

        [TestMethod]
        public void AdaptiveIntegrateInfiniteExpTest11() {
            ddouble expected = -ddouble.Sqrt(ddouble.Pi);

            Assert.AreEqual(0d, (double)(expected - GaussKronrodIntegral.AdaptiveIntegrate(x => ddouble.Exp(-x * x), ddouble.PositiveInfinity, ddouble.NegativeInfinity, eps: 0, order: GaussKronrodOrder.G32K65, maxdepth: 10).value), 1e-30);
        }

        [TestMethod]
        public void LimitedEvalIntegrateExpTest() {
            (ddouble y, ddouble err, long eval_points) = GaussKronrodIntegral.AdaptiveIntegrate(x => ddouble.Exp(ddouble.Sin(16 * x)), 0, 4, 1e-25, discontinue_eval_points: 1024);

            Console.WriteLine(y);
            Console.WriteLine(err);
            Console.WriteLine(eval_points);

            Assert.IsTrue(ddouble.Abs(y - "5.09923585008566829954993402287") < 1e-29);
            Assert.IsTrue(err < 1e-25);
        }

        [TestMethod]
        public void LimitedDepthIntegrateExpTest() {
            (ddouble y, ddouble err, long eval_points) = GaussKronrodIntegral.AdaptiveIntegrate(x => ddouble.Exp(ddouble.Sin(16 * x)), 0, 4, 1e-25, maxdepth: 10);

            Console.WriteLine(y);
            Console.WriteLine(err);
            Console.WriteLine(eval_points);

            Assert.IsTrue(ddouble.Abs(y - "5.09923585008566829954993402287") < 1e-29);
            Assert.IsTrue(err < 1e-25);
        }

        [TestMethod]
        public void LimitedEvalAndDepthIntegrateExpTest() {
            (ddouble y, ddouble err, long eval_points) = GaussKronrodIntegral.AdaptiveIntegrate(x => ddouble.Exp(ddouble.Sin(16 * x)), 0, 4, 1e-25, maxdepth: 10, discontinue_eval_points: 1024);

            Console.WriteLine(y);
            Console.WriteLine(err);
            Console.WriteLine(eval_points);

            Assert.IsTrue(ddouble.Abs(y - "5.09923585008566829954993402287") < 1e-29);
            Assert.IsTrue(err < 1e-25);
        }

        [TestMethod]
        public void UnlimitedIntegrateExpTest() {
            (ddouble y, ddouble err, long eval_points) = GaussKronrodIntegral.AdaptiveIntegrate(x => ddouble.Exp(ddouble.Sin(16 * x)), 0, 4, 1e-25);

            Console.WriteLine(y);
            Console.WriteLine(err);
            Console.WriteLine(eval_points);

            Assert.IsTrue(ddouble.Abs(y - "5.09923585008566829954993402287") < 1e-29);
            Assert.IsTrue(err < 1e-25);
        }

        [TestMethod]
        public void LimitedEvalIntegrateExpIncrementTest() {
            for (int discontinue_eval_points = 16; discontinue_eval_points <= 1024; discontinue_eval_points *= 2) {
                (ddouble y, ddouble err, long eval_points) = GaussKronrodIntegral.AdaptiveIntegrate(x => ddouble.Exp(ddouble.Sin(16 * x)), 0, 4, 1e-25, discontinue_eval_points: discontinue_eval_points);

                Console.WriteLine(y);
                Console.WriteLine(err);
                Console.WriteLine(eval_points);
            }
        }
    }
}
