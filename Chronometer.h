//
// Created by mauri on 14.05.18.
//

#ifndef SIDECHANNELATTACKKANGAROO_CHRONOMETER_H
#define SIDECHANNELATTACKKANGAROO_CHRONOMETER_H

#include <NTL/ZZ.h>

class Chronometer {
public:
    long calculateTime(const NTL::ZZ &a, const NTL::ZZ &x, const NTL::ZZ &n);
};

long Chronometer::calculateTime(const NTL::ZZ &a, const NTL::ZZ &x, const NTL::ZZ &n) {
    auto begin = std::chrono::high_resolution_clock::now();
    uint32_t iterations = 10000;
    for (uint32_t i = 0; i < iterations; ++i) {
        NTL::PowerMod(a, x, n);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
    auto calcTime = duration / iterations;
    std::cout << duration << "ns total, average : " << duration / iterations << "ns." << std::endl;
    return calcTime;
}

#endif //SIDECHANNELATTACKKANGAROO_CHRONOMETER_H
