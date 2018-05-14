//
// Created by mauri on 10.05.18.
//

#ifndef LOGSPICKING_WEINERKANGAROO_H
#define LOGSPICKING_WEINERKANGAROO_H
#define jumps_info std::vector<NTL::ZZ>

#include <NTL/ZZ.h>
#include <vector>
#include "../DiscreteLog.h"

class WeinerKangaroo : public DiscreteLog {
public:
    WeinerKangaroo() {
        x = NTL::RandomBnd(lower - upper) + upper;
        alpha = this->setAlpha();
        beta = NTL::PowerMod(alpha, x, N);
        this->fillSets();
    }

    NTL::ZZ searchCollisions();

    void setBoundries(NTL::ZZ u, NTL::ZZ d);

private:

    const long length = 45;
    const long numProcess = 8;
    const long numWild = numProcess / 2 + 1; // v
    const long numTamed = numProcess / 2 - 1; // u
    const long GAMMA = 65536;
    const long OMEGA = 0xffff;
    NTL::ZZ lower = NTL::RandomBnd(N - (n * M_PI / numProcess / 2)); // a
    NTL::ZZ upper = lower + n * M_PI / numProcess / 2; // b
    NTL::ZZ median = (upper - lower) / 2; // (a+b)/2
    unsigned long k = static_cast<unsigned long>(std::floor(NTL::log(upper - lower)));
    jumps_info distanceS; // S
    jumps_info jumpSetR; // R
    std::hash<std::string> hash_fn;

    void performStep(NTL::ZZ &x, NTL::ZZ &distance, NTL::ZZ &steps);

    std::string zToString(const NTL::ZZ &z);

    inline void fillSets() {
        NTL::ZZ temp;
        distanceS.push_back(NTL::ZZ(numTamed * numWild));
        jumpSetR.push_back(NTL::PowerMod(alpha, distanceS.at(0), N));

        for (int i = 0; i < k; i++) {
            temp = distanceS.back() * 2;
            distanceS.push_back(temp);
            temp = NTL::PowerMod(alpha, temp, N);
            jumpSetR.push_back(temp);
        }
    }


    inline void initKangaroo(int ind, NTL::ZZ &x, NTL::ZZ &distance, NTL::ZZ &steps, long &index, bool &wild);

    inline NTL::ZZ setAlpha() {
        NTL::ZZ tmp_alpha = NTL::ZZ(2);
        do {
            tmp_alpha = NTL::PowerMod(tmp_alpha, 2, N);
        } while (tmp_alpha == NTL::ZZ(1));
        return tmp_alpha;
    }
};

struct DistinguishedPoint {
    NTL::ZZ distance;
    bool isWild;
    long id;
};

#endif //LOGSPICKING_WEINERKANGAROO_H
