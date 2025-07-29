// sobol.h
#ifndef SOBOL_H
#define SOBOL_H

#include <vector>
#include <cstddef>

namespace sobol {

class SobolGenerator {
public:
    SobolGenerator(unsigned long dim = 1);
    void nextPoint(std::vector<double>& point);

private:
    unsigned long dimension;
    unsigned long count;
    std::vector<unsigned long> X;
    std::vector< std::vector<unsigned long> > V;
};

} // namespace sobol

#endif
