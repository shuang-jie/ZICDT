# ZICDT: Directed Tree Structure Learning for Zero-Inflated Compositional Nodes

Each node is a whole composition on the probability simplex. The method
scores candidate edges with a KL divergence (finite on exact zeros),
models the child conditional mean as a mixture of a baseline and a
parent-driven, column-stochastic transition map, fits parameters by EM,
and finds the globally optimal directed tree/forest via Chu-Liu/Edmonds.
The edge penalty is chosen by cross-validation.

## See also

Useful links:

- <https://github.com/sz333024/ZICDT>

- Report bugs at <https://github.com/sz333024/ZICDT/issues>

## Author

**Maintainer**: Shuangjie Zhang <sz333024@gmail.com>
