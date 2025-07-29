// sobol.cpp
#include "sobol.h"
#include <vector>
#include <cstddef>
#include <cmath>

// Direction numbers for the first 2 dimensions (can be extended for more dims)
static const unsigned long mdeg[] = {1, 2};
static const unsigned long ip[]   = {0, 1};
static const unsigned long m[]    = {1, 1, 3};

namespace sobol {

SobolGenerator::SobolGenerator(unsigned long dim)
    : dimension(dim), count(0), X(dim, 0), V(dim)
{
    for (unsigned long d = 0; d < dim; ++d) {
        V[d].resize(32);
        if (d == 0) {
            for (unsigned long i = 0; i < 32; ++i)
                V[d][i] = 1UL << (31 - i);
        } else if (d == 1) {
            V[d][0] = 1UL << 31;
            for (unsigned long i = 1; i < 32; ++i)
                V[d][i] = V[d][i-1] ^ (V[d][i-1] >> 1);
        } else {
            for (unsigned long i = 0; i < 32; ++i)
                V[d][i] = 1UL << (31 - i); // fallback, can be replaced with more direction numbers
        }
    }
}

void SobolGenerator::nextPoint(std::vector<double>& point) {
    unsigned long c = count++;
    unsigned long value;
    for (unsigned long d = 0; d < dimension; ++d) {
        if (c == 0) {
            X[d] = 0;
        } else {
            unsigned long i = 0, c1 = c;
            while (c1 & 1) {
                c1 >>= 1;
                ++i;
            }
            X[d] ^= V[d][i];
        }
        point[d] = static_cast<double>(X[d]) / std::pow(2.0, 32);
    }
}

} // namespace sobol
