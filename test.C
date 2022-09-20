#include "CompileTimeDerivatives.h"
#include <chrono>
#include <execution>
#include <algorithm>
#include <numeric>

using namespace libMesh;

int main()
{
  enum
  {
    dX
  };

  std::vector<Real> in(980);
  for (unsigned int i = 10; i < 990; ++i)
    in[i - 10] = i / 1000.0;

  std::chrono::time_point<std::chrono::system_clock> start, end;

  // evaluate expression template
  start = std::chrono::system_clock::now();

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
                                    for (v = v0; v <= v0 + 1.0 / 1000.0; v += 1e-9)
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
                                    for (v = v0; v <= v0 + 1.0 / 1000.0; v += 1e-9)
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
                                    for (v = v0; v <= v0 + 1.0 / 1000.0; v += 1e-9)
                                      s += result.D<dX>().D<dX>()();
                                    return s;
                                  });

  end = std::chrono::system_clock::now();
  std::cout << std::chrono::duration<double>(end - start).count() << "s\n";

  // evaluate native expression
  start = std::chrono::system_clock::now();

  Real s0 = std::transform_reduce(std::execution::par,
                                  in.begin(), in.end(),
                                  0.0, std::plus<Real>(),
                                  [](const Real v0)
                                  {
                                    Real v;
                                    Real s = 0.0;
                                    for (v = v0; v <= v0 + 1.0 / 1000.0; v += 1e-9)
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
                                    for (v = v0; v <= v0 + 1.0 / 1000.0; v += 1e-9)
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
                                    for (v = v0; v <= v0 + 1.0 / 1000.0; v += 1e-9)
                                      s += -2.0 + 1.0 / (v - 1.0) - 1.0 / v;
                                    return s;
                                  });

  end = std::chrono::system_clock::now();
  std::cout << std::chrono::duration<double>(end - start).count() << "s\n";

  std::cout << "Diffs: " << (r0 - s0) << ", " << (r1 - s1) << ", " << (r2 - s2) << '\n';

  return 0;
}
