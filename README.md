# DoubleDoubleIntegrate
 Double-Double Numerical Integration Implements

## Requirement
.NET 6.0

## Install

[Download DLL](https://github.com/tk-yoshimura/DoubleDoubleIntegrate/releases)  
[Download Nuget](https://www.nuget.org/packages/tyoshimura.doubledouble.integrate/)  

- Import DoubleDouble(https://github.com/tk-yoshimura/DoubleDouble)

## Usage
```csharp
static ddouble f(ddouble x) => ddouble.Sqrt(1 - x * x);

ddouble v = RombergIntegral.Integrate(f, 0, ddouble.Sqrt(2) / 2, level: 20);
Assert.AreEqual(0, (double)((ddouble.PI + 2) / 8 - v), 1e-20);
```

## Licence
[MIT](https://github.com/tk-yoshimura/DoubleDoubleIntegrate/blob/main/LICENSE)

## Author

[T.Yoshimura](https://github.com/tk-yoshimura)
