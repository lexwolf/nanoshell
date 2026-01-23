/*
 * This file is part of the Nano-Shell Simulation Project.
 * 
 * Copyright (C) 2025 Alessandro Veltri
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <vector>
#include <cmath>

double interpolate(double x1, double x2, double f1, double f2, double x);

std::pair<double, double> find_zeros(const std::vector<double>& x, const std::vector<double>& fx, double epsilon = 1e-6, int max_iterations = 100) {
    double zero_1 = 0.0;
    double zero_2 = 0.0;
    int iterations = 0;

    for (int i = 0; i < x.size() - 1; i++) {
        double x1 = x[i];
        double x2 = x[i + 1];
        double f1 = fx[i];
        double f2 = fx[i + 1];

        if (f1 * f2 < 0) {
            double a = x1;
            double b = x2;
            double c = a;

            while (std::abs(b - a) > epsilon && iterations < max_iterations) {
                iterations++;
                c = (a + b) / 2.0;
                double fc = interpolate(x1, x2, f1, f2, c);

                if (f1 * fc < 0) {
                    b = c;
                }
                else {
                    a = c;
                }
            }

            if (zero_1 == 0.0) {
                zero_1 = c;
            }
            else {
                zero_2 = c;
                break;
            }
        }
    }

    if (zero_1 != 0.0 && zero_2 != 0.0) {
        return std::make_pair(zero_1, zero_2);
    }
    else {
        return std::make_pair(0.0, 0.0);
    }
}

double interpolate(double x1, double x2, double f1, double f2, double x) {
    // Perform linear interpolation
    return f1 + (f2 - f1) * (x - x1) / (x2 - x1);
}

int main() {
    std::vector<double> x = { /* Values of x */ };
    std::vector<double> fx = { /* Values of f(x) */ };

    std::pair<double, double> zeros = find_zeros(x, fx);
    std::cout << "Zero 1: " << zeros.first << std::endl;
    std::cout << "Zero 2: " << zeros.second << std::endl;

    return 0;
}
