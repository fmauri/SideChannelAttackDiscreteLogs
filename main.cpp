#include <iostream>
#include "DiscreteLog.h"
#include "WeinerKangaroo/WeinerKangaroo.h"
#include "Chronometer.h"

int main() {
    WeinerKangaroo weinerKangaroo;
    const NTL::ZZ &alpha = weinerKangaroo.getAlpha();
    const NTL::ZZ &beta = weinerKangaroo.getBeta();
    const NTL::ZZ &x = weinerKangaroo.getX();
    const NTL::ZZ &N = weinerKangaroo.getN();
    NTL::ZZ upper, lower, tmp;
    Chronometer chronometer;
    auto mainTime = chronometer.calculateTime(alpha, x, N);

    NTL::ZZ tempX = NTL::RandomBnd(N);
    auto timing = chronometer.calculateTime(alpha, tempX, N);
    upper = 0;
    lower = 0;
    do {
        if (timing > mainTime) {
            upper = tempX;
            do {
                tempX = NTL::RandomBnd(tempX);
                timing = chronometer.calculateTime(alpha, tempX, N);
            } while (timing > mainTime);
            if (timing == mainTime || tempX == x)
                std::cout << x << " = " << tempX << std::endl;
            lower = tempX;
        } else {
            lower = tempX;
            do {
                tempX = NTL::RandomBnd(N) + tempX;
                timing = chronometer.calculateTime(alpha, tempX, N);
            } while (timing < mainTime);
            if (timing == mainTime || tempX == x)
                std::cout << x << " = " << tempX << std::endl;
            upper = tempX;
        }
    } while (upper <= lower);
    weinerKangaroo.setBoundries(upper, lower);
    weinerKangaroo.printLog();
    std::cout << "x = " << weinerKangaroo.searchCollisions() << std::endl;
    return 0;
}

