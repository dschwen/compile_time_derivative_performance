#include "CompileTimeDerivatives.h"
#include <chrono>
#include <execution>
#include <algorithm>
#include <numeric>

using namespace libMesh;

void for_loop()
{
  std::cout << "for_loop()\n";
  using namespace CompileTimeDerivatives;

  enum
  {
    dX
  };

  Real v = 0.0;
  const auto x = makeRef<dX>(v);
  const auto result = x * (1.0 - x) - (x * log(x) + (1.0 - x) * log(1.0 - x));

  Real r0 = 0, r1 = 0, r2 = 0;
  Real s0 = 0, s1 = 0, s2 = 0;

  std::chrono::time_point<std::chrono::system_clock> start, end;

  // evaluate expression template
  start = std::chrono::system_clock::now();
  for (v = 0.01; v <= 0.99; v += 1e-9)
  {
    r0 += result();
    r1 += result.D<dX>()();
    r2 += result.D<dX>().D<dX>()();
  }
  end = std::chrono::system_clock::now();
  std::cout << "Expression template : " << std::chrono::duration<double>(end - start).count() << "s\n";

  // evaluate native expression
  start = std::chrono::system_clock::now();
  for (v = 0.01; v <= 0.99; v += 1e-9)
  {
    s0 += v * (1.0 - v) - (v * std::log(v) + (1.0 - v) * std::log(1.0 - v));
    s1 += -2.0 * v - std::log(v) + std::log(1.0 - v) - (v - 1.0) / (1.0 - v);
    s2 += -2.0 + 1.0 / (v - 1.0) - 1.0 / v;
  }
  end = std::chrono::system_clock::now();
  std::cout << "Hard coded expression: " << std::chrono::duration<double>(end - start).count() << "s\n";

  int n = 0;
  for (v = 0.01; v <= 0.99; v += 1e-9)
    n++;

  std::cout << "Iterations: " << n << '\n';
  std::cout << "Diffs: " << (r0 - s0) << ", " << (r1 - s1) << ", " << (r2 - s2) << '\n';
  std::cout << "Vals: " << r0 << ", " << r1 << ", " << r2 << '\n';
}

void transform_reduce()
{
  std::cout << "transform_reduce()\n";

  enum
  {
    dX
  };

  std::vector<Real> in(98000);
  for (unsigned int i = 1000; i < 99000; ++i)
    in[i - 1000] = i / 100000.0;

  std::chrono::time_point<std::chrono::system_clock> start, end;

  // evaluate expression template
  start = std::chrono::system_clock::now();

  int n = std::transform_reduce(std::execution::par,
                                in.begin(), in.end(),
                                0, std::plus<int>(),
                                [](const Real v0)
                                {
                                  Real v;
                                  int s = 0;
                                  for (v = v0; v <= v0 + 1.0 / 100000.0; v += 1e-9)
                                    s++;
                                  return s;
                                });

  Real r0 = std::transform_reduce(std::execution::par,
                                  in.begin(), in.end(),
                                  0.0, std::plus<Real>(),
                                  [](const Real v0)
                                  {
                                    using namespace CompileTimeDerivatives;
                                    Real v;
                                    const auto x = makeRef<dX>(v);
                                    const auto result = x * (1.0 - x) - (x * log(x) + (1.0 - x) * log(1.0 - x));
                                    Real s = 0.0;
                                    for (v = v0; v <= v0 + 1.0 / 100000.0; v += 1e-9)
                                      s += result();
                                    return s;
                                  });
  Real r1 = std::transform_reduce(std::execution::par,
                                  in.begin(), in.end(),
                                  0.0, std::plus<Real>(),
                                  [](const Real v0)
                                  {
                                    using namespace CompileTimeDerivatives;
                                    Real v;
                                    const auto x = makeRef<dX>(v);
                                    const auto result = x * (1.0 - x) - (x * log(x) + (1.0 - x) * log(1.0 - x));
                                    Real s = 0.0;
                                    for (v = v0; v <= v0 + 1.0 / 100000.0; v += 1e-9)
                                      s += result.D<dX>()();
                                    return s;
                                  });
  Real r2 = std::transform_reduce(std::execution::par,
                                  in.begin(), in.end(),
                                  0.0, std::plus<Real>(),
                                  [](const Real v0)
                                  {
                                    using namespace CompileTimeDerivatives;
                                    Real v;
                                    const auto x = makeRef<dX>(v);
                                    const auto result = x * (1.0 - x) - (x * log(x) + (1.0 - x) * log(1.0 - x));
                                    Real s = 0.0;
                                    for (v = v0; v <= v0 + 1.0 / 100000.0; v += 1e-9)
                                      s += result.D<dX>().D<dX>()();
                                    return s;
                                  });

  end = std::chrono::system_clock::now();
  std::cout << "Expression template : " << std::chrono::duration<double>(end - start).count() << "s\n";

  // evaluate native expression
  start = std::chrono::system_clock::now();

  Real s0 = std::transform_reduce(std::execution::par,
                                  in.begin(), in.end(),
                                  0.0, std::plus<Real>(),
                                  [](const Real v0)
                                  {
                                    Real v;
                                    Real s = 0.0;
                                    for (v = v0; v <= v0 + 1.0 / 100000.0; v += 1e-9)
                                      s += v * (1.0 - v) - (v * std::log(v) + (1.0 - v) * std::log(1.0 - v));
                                    return s;
                                  });
  Real s1 = std::transform_reduce(std::execution::par,
                                  in.begin(), in.end(),
                                  0.0, std::plus<Real>(),
                                  [](const Real v0)
                                  {
                                    Real v;
                                    Real s = 0.0;
                                    for (v = v0; v <= v0 + 1.0 / 100000.0; v += 1e-9)
                                      s += -2.0 * v - std::log(v) + std::log(1.0 - v) - (v - 1.0) / (1.0 - v);
                                    return s;
                                  });
  Real s2 = std::transform_reduce(std::execution::par,
                                  in.begin(), in.end(),
                                  0.0, std::plus<Real>(),
                                  [](const Real v0)
                                  {
                                    Real v;
                                    Real s = 0.0;
                                    for (v = v0; v <= v0 + 1.0 / 100000.0; v += 1e-9)
                                      s += -2.0 + 1.0 / (v - 1.0) - 1.0 / v;
                                    return s;
                                  });

  end = std::chrono::system_clock::now();
  std::cout << "Hard coded expression: " << std::chrono::duration<double>(end - start).count() << "s\n";

  std::cout << "Iterations: " << n << '\n';
  std::cout << "Diffs: " << (r0 - s0) << ", " << (r1 - s1) << ", " << (r2 - s2) << '\n';
  std::cout << "Vals: " << r0 << ", " << r1 << ", " << r2 << '\n';
}

int main()
{
  for_loop();
  transform_reduce();

  return 0;
}
