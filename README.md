# DoubleDoubleIntegrate
 Double-Double Numerical Integration Implements 

## Requirement
.NET 7.0

## Install

[Download DLL](https://github.com/tk-yoshimura/DoubleDoubleIntegrate/releases)  
[Download Nuget](https://www.nuget.org/packages/tyoshimura.doubledouble.integrate/)  

- Import DoubleDouble(https://github.com/tk-yoshimura/DoubleDouble)

## Usage
```csharp
// Gauss-Legendre Integrate 32 Points: sin(t) t=0 to pi
GaussLegendreIntegral.Integrate(
    ddouble.Sin, 
    ddouble.Zero, ddouble.PI, 
    n: 32
);

// Gauss-Kronrod Adaptive Integrate 7-15: exp(t) t=1 to 4
GaussKronrodIntegral.AdaptiveIntegrate(
    ddouble.Exp, 
    1, 4, 
    eps: 1e-25, 
    order: GaussKronrodOrder.G7K15, 
    depth: 10
);

// Gauss-Kronrod Adaptive Integrate 32-65: exp(-t^2) t=-inf to +inf
GaussKronrodIntegral.AdaptiveIntegrate(
    x => ddouble.Exp(-x * x), 
    ddouble.NegativeInfinity, ddouble.PositiveInfinity, 
    eps: 1e-25, 
    order: GaussKronrodOrder.G32K65, 
    depth: 10
);

// Romberg Integrate: sqrt(1 - t^2) t=0 to sqrt(2)/2
RombergIntegral.Integrate(
    x => ddouble.Sqrt(1 - x * x), 
    0, ddouble.Sqrt(2) / 2, 
    level: 20
);
```

## Licence
[MIT](https://github.com/tk-yoshimura/DoubleDoubleIntegrate/blob/main/LICENSE)

## Author

[T.Yoshimura](https://github.com/tk-yoshimura)
