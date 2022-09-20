#include "CompileTimeDerivatives.h"
#include <chrono>
#include <execution>
#include <algorithm>

using namespace libMesh;
using namespace CompileTimeDerivatives;

class FloatRange
{
public:
  class iterator
  {
  public:
    iterator(Real i, Real di) : _i(i), _di(di) {}

    Real operator*() const { return _i; }

    const iterator &operator++()
    {
      _i += _di;
      return *this;
    }

    iterator operator++(int)
    {
      iterator returnval(*this);
      ++_i;
      return returnval;
    }

    bool operator==(const iterator &j) const
    {
      return (_i >= j._i);
    }

    bool operator!=(const iterator &j) const
    {
      return !(*this == j);
    }

  private:
    Real _i, _di;
  };

  FloatRange(Real begin, Real end, Real step) : _begin(begin, step),
                                                _end(end, step)
  {
  }

  iterator begin() const { return _begin; }
  iterator end() const { return _end; }

private:
  iterator _begin, _end;
};

int main()
{
  enum
  {
    dX
  };

  const auto range = FloatRange(0.01, 0.99, 1e-9);

  Real r0 = 0, r1 = 0, r2 = 0;
  Real s0 = 0, s1 = 0, s2 = 0;

  std::chrono::time_point<std::chrono::system_clock> start, end;

  // evaluate expression template
  start = std::chrono::system_clock::now();

  for (auto v : range)
  {
    const auto x = makeRef<dX>(v);
    const auto result = x * (1.0 - x) - (x * log(x) + (1.0 - x) * log(1.0 - x));

    r0 += result();
    r1 += result.D<dX>()();
    r2 += result.D<dX>().D<dX>()();
  }
  end = std::chrono::system_clock::now();
  std::cout << std::chrono::duration<double>(end - start).count() << "s\n";

  // evaluate native expression
  start = std::chrono::system_clock::now();
  for (auto v : range)
  {
    const auto x = makeRef<dX>(v);
    const auto result = x * (1.0 - x) - (x * log(x) + (1.0 - x) * log(1.0 - x));

    s0 += v * (1.0 - v) - (v * std::log(v) + (1.0 - v) * std::log(1.0 - v));
    s1 += -2.0 * v - std::log(v) + std::log(1.0 - v) - (v - 1.0) / (1.0 - v);
    s2 += -2.0 + 1.0 / (v - 1.0) - 1.0 / v;
  }
  end = std::chrono::system_clock::now();
  std::cout << std::chrono::duration<double>(end - start).count() << "s\n";

  std::cout << "Diffs: " << (r0 - s0) << ", " << (r1 - s1) << ", " << (r2 - s2) << '\n';

  return 0;
}
