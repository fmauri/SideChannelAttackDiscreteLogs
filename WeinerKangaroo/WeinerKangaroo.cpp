//
// Created by mauri on 10.05.18.
//

#include <omp.h>
#include <map>
#include <sstream>
#include "WeinerKangaroo.h"

NTL::ZZ WeinerKangaroo::searchCollisions() {
    std::map<NTL::ZZ, DistinguishedPoint> points;
    NTL::ZZ tmp_x, distanceWild, distanceTamed, wild, tamed, st, indexTamed, indexWild, tmp, dist;
    NTL::ZZ kangarooDistance, kangarooX, kangarooSteps, result;
    DistinguishedPoint tmp_result;

    result = NTL::ZZ(0);
    distanceTamed = 0;
    distanceWild = 0;

    //////////////////////Kangaroo/////////////////////////////////
    NTL::ZZ k_x;
    NTL::ZZ k_distance;
    NTL::ZZ k_steps;
    long k_index;
    bool k_isWild;
    /////////////////////////////////////////////////////////////
    /*
     * threads from 0 to 2 are Tamed and from 3 to 7 are Wild
     */
    omp_set_num_threads(numProcess);
#pragma omp parallel for private(distanceTamed, distanceWild, indexTamed, indexWild, tmp, k_x, k_distance, k_steps, k_index, k_isWild, tmp_result) shared(result, points) schedule(dynamic)
    for (int i = 0; i < omp_get_num_threads(); ++i) {
        initKangaroo(omp_get_thread_num(), k_x, k_distance, k_steps, k_index, k_isWild);
        while (result == 0) {
            performStep(k_x, k_distance, k_steps);
            if (k_steps > (GAMMA * 20)) { break; }
            if (k_x % GAMMA != 0) continue; // is a distinguished point ?
            auto point = points.find(k_x);
            if (point != points.end()) { // do we have a collision
                if (result == 0 && (k_isWild != point->second.isWild)) {
                    if (k_isWild) {
                        distanceWild = k_distance;
                        distanceTamed = point->second.distance;
                        indexTamed = point->second.id;
                        indexWild = k_index;
                    } else {
                        distanceWild = point->second.distance;
                        distanceTamed = k_distance;;
                        indexTamed = k_index;
                        indexWild = point->second.id;
                    }
                    tmp = (median + indexTamed * numWild - indexWild * numTamed + distanceTamed - distanceWild) % n;
                    if (result != 0 && NTL::PowerMod(alpha, tmp, N) != beta) continue;
                    result = tmp;
                }
            } else {
                tmp_result = {k_distance, k_isWild, k_index};
                points.insert({k_x, tmp_result});
            }
        }
    }
    return result != 0 ? result : NTL::ZZ(-5);
}

void WeinerKangaroo::initKangaroo(int ind, NTL::ZZ &x, NTL::ZZ &distance, NTL::ZZ &steps, long &index, bool &wild) {
    if (ind < numTamed) { //tamed init
        wild = false;
        index = ind;
        x = NTL::PowerMod(alpha, (median + index * numWild), N);
    } else { // wild init
        wild = true;
        index = ind - numTamed;
        x = NTL::MulMod(beta, NTL::PowerMod(alpha, index * numTamed, N), N);
    }
    distance = NTL::ZZ(0);
    steps = NTL::ZZ(0);
}

void WeinerKangaroo::performStep(NTL::ZZ &x, NTL::ZZ &distance, NTL::ZZ &steps) {
    std::string temp = zToString(x);
    size_t i = hash_fn(temp) % k; // x % k
    x = NTL::MulMod(x, jumpSetR.at(i), N);
    distance = NTL::AddMod(distance, distanceS.at(i), n);
    steps += 1;
}

std::string WeinerKangaroo::zToString(const NTL::ZZ &z) {
    std::stringstream buffer;
    buffer << z;
    return buffer.str();
}

void WeinerKangaroo::setBoundries(NTL::ZZ u, NTL::ZZ d) {
//    if (u < this->x || d > this->x) {
//        std::cout << "Side-channel failed " << u << " == " << this->x << " " << d << " == " << this->x << "\n";
//        throw std::invalid_argument("Side-channel attack failed");
//    }
    this->upper = std::move(u);
    this->lower = std::move(d);
    this->k = static_cast<unsigned long>(std::floor(NTL::log(this->upper - this->lower)));
    this->median = (this->upper - this->lower) / 2;
}
