#include "CompileTimeDerivatives.h"
#include <chrono>

using namespace libMesh;
using namespace CompileTimeDerivatives;

int main()
{
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
  std::cout << std::chrono::duration<double>(end - start).count() << "s\n";

  // evaluate native expression
  start = std::chrono::system_clock::now();
  for (v = 0.01; v <= 0.99; v += 1e-9)
  {
    s0 += v * (1.0 - v) - (v * std::log(v) + (1.0 - v) * std::log(1.0 - v));
    s1 += -2.0 * v - std::log(v) + std::log(1.0 - v) - (v - 1.0) / (1.0 - v);
    s2 += -2.0 + 1.0 / (v - 1.0) - 1.0 / v;
  }
  end = std::chrono::system_clock::now();
  std::cout << std::chrono::duration<double>(end - start).count() << "s\n";

  std::cout << "Diffs: " << (r0 - s0) << ", " << (r1 - s1) << ", " << (r2 - s2) << '\n';

  return 0;
}
