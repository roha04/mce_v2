#include "Gauss.h"
#include <cmath>
#include <iostream>

// Координати та ваги Гауса для 3 точок (з 6 занять)
const std::array<double, 3> Gauss::points = {-std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0)};
const std::array<double, 3> Gauss::weights = {5.0/9.0, 8.0/9.0, 5.0/9.0};

// Локальні вузли на гранях 20-вузлового елемента
// Грані: 0=xi=-1, 1=xi=1, 2=eta=-1, 3=eta=1, 4=zeta=-1, 5=zeta=1
const std::array<std::array<int, 8>, 6> Gauss::faceNodes = {{
    // Грань 0: xi = -1 (вузли 0,3,7,4 та серединні 11,19,15,16)
    {0, 3, 7, 4, 11, 19, 15, 16},
    // Грань 1: xi = 1 (вузли 1,2,6,5 та серединні 9,18,13,17)
    {1, 2, 6, 5, 9, 18, 13, 17},
    // Грань 2: eta = -1 (вузли 0,1,5,4 та серединні 8,17,12,16)
    {0, 1, 5, 4, 8, 17, 12, 16},
    // Грань 3: eta = 1 (вузли 2,3,7,6 та серединні 10,19,14,18)
    {2, 3, 7, 6, 10, 19, 14, 18},
    // Грань 4: zeta = -1 (вузли 0,1,2,3 та серединні 8,9,10,11)
    {0, 1, 2, 3, 8, 9, 10, 11},
    // Грань 5: zeta = 1 (вузли 4,5,6,7 та серединні 12,13,14,15)
    {4, 5, 6, 7, 12, 13, 14, 15}
}};

void Gauss::shapeFunctions(double xi, double eta, double zeta,
                           std::array<double, 20>& N) {
    // Допоміжні величини: (1 ± xi), (1 ± eta), (1 ± zeta)
    double pxi = 1.0 + xi;   double mxi = 1.0 - xi;
    double peta = 1.0 + eta; double meta = 1.0 - eta;
    double pzeta = 1.0 + zeta; double mzeta = 1.0 - zeta;

    // --- Кутові вузли (1-8) ---
    // Для кутового вузла (xi0, eta0, zeta0), де xi0, eta0, zeta0 = ±1
    // N_i = 1/8 * (1 + xi*xi0)(1 + eta*eta0)(1 + zeta*zeta0) * (xi*xi0 + eta*eta0 + zeta*zeta0 - 2)

    // Вузол 0: (-1, -1, -1)
    N[0] = 0.125 * mxi * meta * mzeta * (-xi - eta - zeta - 2.0);

    // Вузол 1: (1, -1, -1)
    N[1] = 0.125 * pxi * meta * mzeta * (xi - eta - zeta - 2.0);

    // Вузол 2: (1, 1, -1)
    N[2] = 0.125 * pxi * peta * mzeta * (xi + eta - zeta - 2.0);

    // Вузол 3: (-1, 1, -1)
    N[3] = 0.125 * mxi * peta * mzeta * (-xi + eta - zeta - 2.0);

    // Вузол 4: (-1, -1, 1)
    N[4] = 0.125 * mxi * meta * pzeta * (-xi - eta + zeta - 2.0);

    // Вузол 5: (1, -1, 1)
    N[5] = 0.125 * pxi * meta * pzeta * (xi - eta + zeta - 2.0);

    // Вузол 6: (1, 1, 1)
    N[6] = 0.125 * pxi * peta * pzeta * (xi + eta + zeta - 2.0);

    // Вузол 7: (-1, 1, 1)
    N[7] = 0.125 * mxi * peta * pzeta * (-xi + eta + zeta - 2.0);

    // --- Серединні вузли на ребрах (9-20) ---
    // Для вузла на ребрі, де змінюється тільки одна координата,
    // N_i = 1/4 * (1 - xi^2) * (1 + eta*eta0) * (1 + zeta*zeta0)
    // (аналогічно для інших комбінацій)

    // Ребра нижньої грані (zeta = -1)
    // Вузол 8: (0, -1, -1) - ребро 1-2
    N[8] = 0.25 * (1.0 - xi*xi) * meta * mzeta;

    // Вузол 9: (1, 0, -1) - ребро 2-3
    N[9] = 0.25 * pxi * (1.0 - eta*eta) * mzeta;

    // Вузол 10: (0, 1, -1) - ребро 3-4
    N[10] = 0.25 * (1.0 - xi*xi) * peta * mzeta;

    // Вузол 11: (-1, 0, -1) - ребро 4-1
    N[11] = 0.25 * mxi * (1.0 - eta*eta) * mzeta;

    // Ребра верхньої грані (zeta = 1)
    // Вузол 12: (0, -1, 1) - ребро 5-6
    N[12] = 0.25 * (1.0 - xi*xi) * meta * pzeta;

    // Вузол 13: (1, 0, 1) - ребро 6-7
    N[13] = 0.25 * pxi * (1.0 - eta*eta) * pzeta;

    // Вузол 14: (0, 1, 1) - ребро 7-8
    N[14] = 0.25 * (1.0 - xi*xi) * peta * pzeta;

    // Вузол 15: (-1, 0, 1) - ребро 8-5
    N[15] = 0.25 * mxi * (1.0 - eta*eta) * pzeta;

    // Вертикальні ребра (xi, eta фіксовані)
    // Вузол 16: (-1, -1, 0) - ребро 1-5
    N[16] = 0.25 * mxi * meta * (1.0 - zeta*zeta);

    // Вузол 17: (1, -1, 0) - ребро 2-6
    N[17] = 0.25 * pxi * meta * (1.0 - zeta*zeta);

    // Вузол 18: (1, 1, 0) - ребро 3-7
    N[18] = 0.25 * pxi * peta * (1.0 - zeta*zeta);

    // Вузол 19: (-1, 1, 0) - ребро 4-8
    N[19] = 0.25 * mxi * peta * (1.0 - zeta*zeta);
}

void Gauss::shapeFunctionDerivatives(double xi, double eta, double zeta,
                                      std::array<double, 20>& dN_dxi,
                                      std::array<double, 20>& dN_deta,
                                      std::array<double, 20>& dN_dzeta) {
    double pxi = 1.0 + xi;   double mxi = 1.0 - xi;
    double peta = 1.0 + eta; double meta = 1.0 - eta;
    double pzeta = 1.0 + zeta; double mzeta = 1.0 - zeta;

    // --- Похідні за xi ---

    // Кутові вузли
    // d/dxi [1/8 * (1 + xi*xi0)(1 + eta*eta0)(1 + zeta*zeta0) * (xi*xi0 + eta*eta0 + zeta*zeta0 - 2)]
    // = 1/8 * xi0 * (1 + eta*eta0)(1 + zeta*zeta0) * [ (xi*xi0 + eta*eta0 + zeta*zeta0 - 2) + (1 + xi*xi0) ]

    // Вузол 0: (-1, -1, -1)
    dN_dxi[0] = -0.125 * meta * mzeta * (-2.0*xi - eta - zeta - 1.0);
    // Вузол 1: (1, -1, -1)
    dN_dxi[1] =  0.125 * meta * mzeta * (2.0*xi - eta - zeta - 1.0);
    // Вузол 2: (1, 1, -1)
    dN_dxi[2] =  0.125 * peta * mzeta * (2.0*xi + eta - zeta - 1.0);
    // Вузол 3: (-1, 1, -1)
    dN_dxi[3] = -0.125 * peta * mzeta * (-2.0*xi + eta - zeta - 1.0);
    // Вузол 4: (-1, -1, 1)
    dN_dxi[4] = -0.125 * meta * pzeta * (-2.0*xi - eta + zeta - 1.0);
    // Вузол 5: (1, -1, 1)
    dN_dxi[5] =  0.125 * meta * pzeta * (2.0*xi - eta + zeta - 1.0);
    // Вузол 6: (1, 1, 1)
    dN_dxi[6] =  0.125 * peta * pzeta * (2.0*xi + eta + zeta - 1.0);
    // Вузол 7: (-1, 1, 1)
    dN_dxi[7] = -0.125 * peta * pzeta * (-2.0*xi + eta + zeta - 1.0);

    // Серединні вузли - похідні за xi
    // Вузол 8: (0, -1, -1) -> d/dxi [1/4 * (1-xi^2) * (1-eta) * (1-zeta)]
    dN_dxi[8]  = -0.5 * xi * meta * mzeta;
    // Вузол 9: (1, 0, -1) -> d/dxi [1/4 * (1+xi) * (1-eta^2) * (1-zeta)]
    dN_dxi[9]  = 0.25 * (1.0 - eta*eta) * mzeta;
    // Вузол 10: (0, 1, -1)
    dN_dxi[10] = -0.5 * xi * peta * mzeta;
    // Вузол 11: (-1, 0, -1)
    dN_dxi[11] = -0.25 * (1.0 - eta*eta) * mzeta;
    // Вузол 12: (0, -1, 1)
    dN_dxi[12] = -0.5 * xi * meta * pzeta;
    // Вузол 13: (1, 0, 1)
    dN_dxi[13] = 0.25 * (1.0 - eta*eta) * pzeta;
    // Вузол 14: (0, 1, 1)
    dN_dxi[14] = -0.5 * xi * peta * pzeta;
    // Вузол 15: (-1, 0, 1)
    dN_dxi[15] = -0.25 * (1.0 - eta*eta) * pzeta;
    // Вузол 16: (-1, -1, 0)
    dN_dxi[16] = -0.25 * meta * (1.0 - zeta*zeta);
    // Вузол 17: (1, -1, 0)
    dN_dxi[17] =  0.25 * meta * (1.0 - zeta*zeta);
    // Вузол 18: (1, 1, 0)
    dN_dxi[18] =  0.25 * peta * (1.0 - zeta*zeta);
    // Вузол 19: (-1, 1, 0)
    dN_dxi[19] = -0.25 * peta * (1.0 - zeta*zeta);

    // --- Похідні за eta ---

    // Кутові вузли
    dN_deta[0] = -0.125 * mxi * mzeta * (-xi - 2.0*eta - zeta - 1.0);
    dN_deta[1] = -0.125 * pxi * mzeta * (xi - 2.0*eta - zeta - 1.0);
    dN_deta[2] =  0.125 * pxi * mzeta * (xi + 2.0*eta - zeta - 1.0);
    dN_deta[3] =  0.125 * mxi * mzeta * (-xi + 2.0*eta - zeta - 1.0);
    dN_deta[4] = -0.125 * mxi * pzeta * (-xi - 2.0*eta + zeta - 1.0);
    dN_deta[5] = -0.125 * pxi * pzeta * (xi - 2.0*eta + zeta - 1.0);
    dN_deta[6] =  0.125 * pxi * pzeta * (xi + 2.0*eta + zeta - 1.0);
    dN_deta[7] =  0.125 * mxi * pzeta * (-xi + 2.0*eta + zeta - 1.0);

    // Серединні вузли - похідні за eta
    dN_deta[8]  = -0.25 * (1.0 - xi*xi) * mzeta;
    dN_deta[9]  = -0.5 * eta * pxi * mzeta;
    dN_deta[10] =  0.25 * (1.0 - xi*xi) * mzeta;
    dN_deta[11] = -0.5 * eta * mxi * mzeta;
    dN_deta[12] = -0.25 * (1.0 - xi*xi) * pzeta;
    dN_deta[13] = -0.5 * eta * pxi * pzeta;
    dN_deta[14] =  0.25 * (1.0 - xi*xi) * pzeta;
    dN_deta[15] = -0.5 * eta * mxi * pzeta;
    dN_deta[16] = -0.25 * mxi * (1.0 - zeta*zeta);
    dN_deta[17] = -0.25 * pxi * (1.0 - zeta*zeta);
    dN_deta[18] =  0.25 * pxi * (1.0 - zeta*zeta);
    dN_deta[19] =  0.25 * mxi * (1.0 - zeta*zeta);

    // --- Похідні за zeta ---

    // Кутові вузли
    dN_dzeta[0] = -0.125 * mxi * meta * (-xi - eta - 2.0*zeta - 1.0);
    dN_dzeta[1] = -0.125 * pxi * meta * (xi - eta - 2.0*zeta - 1.0);
    dN_dzeta[2] = -0.125 * pxi * peta * (xi + eta - 2.0*zeta - 1.0);
    dN_dzeta[3] = -0.125 * mxi * peta * (-xi + eta - 2.0*zeta - 1.0);
    dN_dzeta[4] =  0.125 * mxi * meta * (-xi - eta + 2.0*zeta - 1.0);
    dN_dzeta[5] =  0.125 * pxi * meta * (xi - eta + 2.0*zeta - 1.0);
    dN_dzeta[6] =  0.125 * pxi * peta * (xi + eta + 2.0*zeta - 1.0);
    dN_dzeta[7] =  0.125 * mxi * peta * (-xi + eta + 2.0*zeta - 1.0);

    // Серединні вузли - похідні за zeta
    dN_dzeta[8]  = -0.25 * (1.0 - xi*xi) * meta;
    dN_dzeta[9]  = -0.25 * pxi * (1.0 - eta*eta);
    dN_dzeta[10] = -0.25 * (1.0 - xi*xi) * peta;
    dN_dzeta[11] = -0.25 * mxi * (1.0 - eta*eta);
    dN_dzeta[12] =  0.25 * (1.0 - xi*xi) * meta;
    dN_dzeta[13] =  0.25 * pxi * (1.0 - eta*eta);
    dN_dzeta[14] =  0.25 * (1.0 - xi*xi) * peta;
    dN_dzeta[15] =  0.25 * mxi * (1.0 - eta*eta);
    dN_dzeta[16] = -0.5 * zeta * mxi * meta;
    dN_dzeta[17] = -0.5 * zeta * pxi * meta;
    dN_dzeta[18] = -0.5 * zeta * pxi * peta;
    dN_dzeta[19] = -0.5 * zeta * mxi * peta;
}

void Gauss::computeDFIABG(std::array<std::array<std::array<double, 20>, 3>, 27>& DFIABG) {
    int idx = 0;
    for (int k = 0; k < 3; ++k) {
        for (int j = 0; j < 3; ++j) {
            for (int i = 0; i < 3; ++i) {
                double xi = points[i];
                double eta = points[j];
                double zeta = points[k];

                std::array<double, 20> dN_dxi, dN_deta, dN_dzeta;
                shapeFunctionDerivatives(xi, eta, zeta, dN_dxi, dN_deta, dN_dzeta);

                // DFIABG[gaussNode][coord][func]
                // coord: 0=xi, 1=eta, 2=zeta
                for (int f = 0; f < 20; ++f) {
                    DFIABG[idx][0][f] = dN_dxi[f];
                    DFIABG[idx][1][f] = dN_deta[f];
                    DFIABG[idx][2][f] = dN_dzeta[f];
                }
                ++idx;
            }
        }
    }
}

// ========== 2D Функції для граней ==========

void Gauss::shapeFunctions2D(double xi, double eta,
                              std::array<double, 8>& N) {
    double pxi = 1.0 + xi;   double mxi = 1.0 - xi;
    double peta = 1.0 + eta; double meta = 1.0 - eta;

    // Кутові вузли квадрата [-1,1]x[-1,1]
    // N_i = 1/4 * (1 + xi*xi0)(1 + eta*eta0) * (xi*xi0 + eta*eta0 - 1)
    N[0] = 0.25 * mxi * meta * (-xi - eta - 1.0);  // (-1, -1)
    N[1] = 0.25 * pxi * meta * (xi - eta - 1.0);   // (1, -1)
    N[2] = 0.25 * pxi * peta * (xi + eta - 1.0);   // (1, 1)
    N[3] = 0.25 * mxi * peta * (-xi + eta - 1.0);  // (-1, 1)

    // Серединні вузли
    N[4] = 0.5 * (1.0 - xi*xi) * meta;   // (0, -1)
    N[5] = 0.5 * pxi * (1.0 - eta*eta);  // (1, 0)
    N[6] = 0.5 * (1.0 - xi*xi) * peta;   // (0, 1)
    N[7] = 0.5 * mxi * (1.0 - eta*eta);  // (-1, 0)
}

void Gauss::shapeFunctionDerivatives2D(double xi, double eta,
                                        std::array<double, 8>& dN_dxi,
                                        std::array<double, 8>& dN_deta) {
    double pxi = 1.0 + xi;   double mxi = 1.0 - xi;
    double peta = 1.0 + eta; double meta = 1.0 - eta;

    // Похідні за xi
    dN_dxi[0] = -0.25 * meta * (-2.0*xi - eta + 0.0);
    dN_dxi[1] =  0.25 * meta * (2.0*xi - eta + 0.0);
    dN_dxi[2] =  0.25 * peta * (2.0*xi + eta + 0.0);
    dN_dxi[3] = -0.25 * peta * (-2.0*xi + eta + 0.0);
    dN_dxi[4] = -xi * meta;
    dN_dxi[5] =  0.5 * (1.0 - eta*eta);
    dN_dxi[6] = -xi * peta;
    dN_dxi[7] = -0.5 * (1.0 - eta*eta);

    // Похідні за eta
    dN_deta[0] = -0.25 * mxi * (-xi - 2.0*eta + 0.0);
    dN_deta[1] = -0.25 * pxi * (xi - 2.0*eta + 0.0);
    dN_deta[2] =  0.25 * pxi * (xi + 2.0*eta + 0.0);
    dN_deta[3] =  0.25 * mxi * (-xi + 2.0*eta + 0.0);
    dN_deta[4] = -0.5 * (1.0 - xi*xi);
    dN_deta[5] = -eta * pxi;
    dN_deta[6] =  0.5 * (1.0 - xi*xi);
    dN_deta[7] = -eta * mxi;
}

void Gauss::computeDPSITE(std::array<std::array<std::array<double, 8>, 2>, 9>& DPSITE) {
    int idx = 0;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            double xi = points[i];
            double eta = points[j];

            std::array<double, 8> dN_dxi, dN_deta;
            shapeFunctionDerivatives2D(xi, eta, dN_dxi, dN_deta);

            // DPSITE[gaussNode][coord][func]
            // coord: 0=xi, 1=eta
            for (int f = 0; f < 8; ++f) {
                DPSITE[idx][0][f] = dN_dxi[f];
                DPSITE[idx][1][f] = dN_deta[f];
            }
            ++idx;
        }
    }
}