/*
 * Generic tools for working with Z(x) series stored as
 * vector of (x, complex Z) pairs.
 */

#ifndef ZX_TOOLS_H
#define ZX_TOOLS_H

#include <algorithm>
#include <cmath>
#include <complex>
#include <utility>
#include <vector>

using ZxSeries = std::vector<std::pair<double, std::complex<double>>>;

inline std::vector<double> extract_x(const ZxSeries& vectorOfPairs) {
    std::vector<double> extractFirst;
    int i = 0;
    for (const auto& pair : vectorOfPairs) {
        i++;
        extractFirst.push_back(pair.first);
        if (i == int(vectorOfPairs.size()) - 1) break;
    }
    return extractFirst;
}

inline std::vector<double> extract_rZ(const ZxSeries& vectorOfPairs) {
    std::vector<double> extractSecnd;
    int i = 0;
    for (const auto& pair : vectorOfPairs) {
        i++;
        extractSecnd.push_back(std::real(pair.second));
        if (i == int(vectorOfPairs.size()) - 1) break;
    }
    return extractSecnd;
}

inline std::vector<double> extract_iZ(const ZxSeries& vectorOfPairs) {
    std::vector<double> extractSecnd;
    int i = 0;
    for (const auto& pair : vectorOfPairs) {
        i++;
        extractSecnd.push_back(std::imag(pair.second));
        if (i == int(vectorOfPairs.size()) - 1) break;
    }
    return extractSecnd;
}

inline std::pair<double, double> find_extrema(const ZxSeries& vkape, double Ome_p) {
    std::vector<double> tmp;
    double extr_1, extr_2;
    int never = 1;
    for (int ii = 0; ii < int(vkape.size()); ii++) {
        if (std::real(vkape[ii].second) >= 0) {
            tmp.push_back(vkape[ii].first - std::imag(vkape[ii].second) * Ome_p);
            never = 0;
        }
    }
    std::sort(tmp.begin(), tmp.end());
    if (never == 0) {
        extr_1 = tmp[0];
        extr_2 = tmp[tmp.size() - 1];
    } else {
        extr_1 = 666;
        extr_2 = 777;
    }
    return std::make_pair(extr_1, extr_2);
}

inline bool compare_by_abs_second(const std::pair<double, double>& a,
                                  const std::pair<double, double>& b) {
    return std::fabs(a.second) > std::fabs(b.second);
}

#endif
