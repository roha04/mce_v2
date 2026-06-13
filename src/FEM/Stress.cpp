#define _USE_MATH_DEFINES
#include "Stress.h"
#include "Element.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <map>

// Коефіцієнти Ляме (повторюємо з Element.cpp для зручності)
static double lambda() {
    double E = Element::YOUNG_MODULUS;
    double nu = Element::POISSON_RATIO;
    return E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
}

static double mu() {
    double E = Element::YOUNG_MODULUS;
    double nu = Element::POISSON_RATIO;
    return E / (2.0 * (1.0 + nu));
}

void Stress::computeStresses(const Mesh& mesh,
                              const std::vector<double>& displacements,
                              const std::array<std::array<std::array<double, 20>, 3>, 27>& DFIABG,
                              std::vector<std::array<double, 6>>& sigma) {
    
    int nqp = mesh.nodes.size();
    if (displacements.size() < 3 * nqp) {
        std::cerr << "Помилка: недостатньо даних про переміщення для обчислення напружень\n";
        return;
    }

    // Ініціалізація: кожному вузлу відповідає сума напружень та кількість внесків
    std::vector<std::array<double, 6>> sigmaSum(nqp, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    std::vector<int> contributions(nqp, 0);

    std::cout << "Обчислення напружень у вузлах...\n";

    // Цикл по елементах
    for (int e = 0; e < mesh.nel; ++e) {
        // Координати вузлів елемента
        std::array<Eigen::Vector3d, 20> elemNodes;
        for (int i = 0; i < 20; ++i) {
            elemNodes[i] = mesh.nodes[mesh.elements[e][i]];
        }

        // Переміщення вузлів елемента (60 значень)
        Eigen::Matrix<double, 60, 1> elemDisp;
        for (int i = 0; i < 20; ++i) {
            int idx = mesh.elements[e][i];
            elemDisp(i)       = displacements[3 * idx];
            elemDisp(i + 20)  = displacements[3 * idx + 1];
            elemDisp(i + 40)  = displacements[3 * idx + 2];
        }

        // Обчислюємо напруження в кожному вузлі Гауса та усереднюємо
        // Для простоти: обчислюємо в центрі елемента (xi=0, eta=0, zeta=0)
        // та додаємо до всіх вузлів
        double xi = 0.0, eta = 0.0, zeta = 0.0;
        
        std::array<double, 20> dN_dxi, dN_deta, dN_dzeta;
        Gauss::shapeFunctionDerivatives(xi, eta, zeta, dN_dxi, dN_deta, dN_dzeta);

        // Якобіан
        Eigen::Matrix3d J;
        J.setZero();
        for (int i = 0; i < 20; ++i) {
            J(0, 0) += dN_dxi[i]   * elemNodes[i].x();
            J(1, 0) += dN_dxi[i]   * elemNodes[i].y();
            J(2, 0) += dN_dxi[i]   * elemNodes[i].z();
            J(0, 1) += dN_deta[i]  * elemNodes[i].x();
            J(1, 1) += dN_deta[i]  * elemNodes[i].y();
            J(2, 1) += dN_deta[i]  * elemNodes[i].z();
            J(0, 2) += dN_dzeta[i] * elemNodes[i].x();
            J(1, 2) += dN_dzeta[i] * elemNodes[i].y();
            J(2, 2) += dN_dzeta[i] * elemNodes[i].z();
        }

        Eigen::Matrix3d Jinv = J.inverse();

        // Похідні за глобальними координатами в центрі елемента
        std::array<double, 20> df_dx, df_dy, df_dz;
        for (int i = 0; i < 20; ++i) {
            df_dx[i] = Jinv(0, 0) * dN_dxi[i] + Jinv(0, 1) * dN_deta[i] + Jinv(0, 2) * dN_dzeta[i];
            df_dy[i] = Jinv(1, 0) * dN_dxi[i] + Jinv(1, 1) * dN_deta[i] + Jinv(1, 2) * dN_dzeta[i];
            df_dz[i] = Jinv(2, 0) * dN_dxi[i] + Jinv(2, 1) * dN_deta[i] + Jinv(2, 2) * dN_dzeta[i];
        }

        // Обчислюємо напруження в центрі елемента
        auto stressPoint = computeStressAtPoint(df_dx, df_dy, df_dz, elemDisp);

        // Додаємо до всіх вузлів елемента
        for (int i = 0; i < 20; ++i) {
            int globalNode = mesh.elements[e][i];
            for (int c = 0; c < 6; ++c) {
                sigmaSum[globalNode][c] += stressPoint[c];
            }
            contributions[globalNode]++;
        }

        if ((e + 1) % 100 == 0 || e == 0 || e == mesh.nel - 1) {
            std::cout << "  Оброблено елементів: " << (e + 1) << " / " << mesh.nel << "\r";
            std::cout.flush();
        }
    }
    std::cout << "\n";

    // Усереднення напружень по елементах
    sigma.resize(nqp);
    for (int i = 0; i < nqp; ++i) {
        if (contributions[i] > 0) {
            for (int c = 0; c < 6; ++c) {
                sigma[i][c] = sigmaSum[i][c] / contributions[i];
            }
        } else {
            sigma[i] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        }
    }

    std::cout << "Напруження обчислені.\n";
}

std::array<double, 6> Stress::computeStressAtPoint(
    const std::array<double, 20>& df_dx,
    const std::array<double, 20>& df_dy,
    const std::array<double, 20>& df_dz,
    const Eigen::Matrix<double, 60, 1>& elemDisp) {
    
    double lam = lambda();
    double m = mu();

    // Обчислюємо деформації через похідні переміщень
    // epsilon_xx = du/dx, epsilon_yy = dv/dy, epsilon_zz = dw/dz
    // gamma_xy = du/dy + dv/dx, gamma_yz = dv/dz + dw/dy, gamma_xz = du/dz + dw/dx
    double du_dx = 0.0, du_dy = 0.0, du_dz = 0.0;
    double dv_dx = 0.0, dv_dy = 0.0, dv_dz = 0.0;
    double dw_dx = 0.0, dw_dy = 0.0, dw_dz = 0.0;

    for (int i = 0; i < 20; ++i) {
        du_dx += df_dx[i] * elemDisp(i);
        du_dy += df_dy[i] * elemDisp(i);
        du_dz += df_dz[i] * elemDisp(i);
        
        dv_dx += df_dx[i] * elemDisp(i + 20);
        dv_dy += df_dy[i] * elemDisp(i + 20);
        dv_dz += df_dz[i] * elemDisp(i + 20);
        
        dw_dx += df_dx[i] * elemDisp(i + 40);
        dw_dy += df_dy[i] * elemDisp(i + 40);
        dw_dz += df_dz[i] * elemDisp(i + 40);
    }

    // Закон Гука для ізотропного тіла (формули 2)
    double theta = du_dx + dv_dy + dw_dz; // об'ємна деформація

    std::array<double, 6> stress;
    stress[0] = lam * theta + 2.0 * m * du_dx;            // sigma_xx
    stress[1] = lam * theta + 2.0 * m * dv_dy;            // sigma_yy
    stress[2] = lam * theta + 2.0 * m * dw_dz;            // sigma_zz
    stress[3] = m * (du_dy + dv_dx);                       // sigma_xy
    stress[4] = m * (dv_dz + dw_dy);                       // sigma_yz
    stress[5] = m * (du_dz + dw_dx);                       // sigma_xz

    return stress;
}

void Stress::computePrincipalStresses(const std::vector<std::array<double, 6>>& sigma,
                                       std::vector<double>& sigma1,
                                       std::vector<double>& sigma2,
                                       std::vector<double>& sigma3) {
    
    int n = sigma.size();
    sigma1.resize(n);
    sigma2.resize(n);
    sigma3.resize(n);

    for (int i = 0; i < n; ++i) {
        double sxx = sigma[i][0];
        double syy = sigma[i][1];
        double szz = sigma[i][2];
        double sxy = sigma[i][3];
        double syz = sigma[i][4];
        double sxz = sigma[i][5];

        // Інваріанти тензора напружень (формули 49)
        double J1 = sxx + syy + szz;
        double J2 = sxx*syy + sxx*szz + syy*szz - (sxy*sxy + syz*syz + sxz*sxz);
        double J3 = sxx*syy*szz + 2.0*sxy*syz*sxz - (sxx*syz*syz + syy*sxz*sxz + szz*sxy*sxy);

        auto principals = solveCubic(J1, J2, J3);
        sigma1[i] = principals[0];
        sigma2[i] = principals[1];
        sigma3[i] = principals[2];
    }
}

void Stress::computeVonMises(const std::vector<std::array<double, 6>>& sigma,
                              std::vector<double>& vonMises) {
    
    int n = sigma.size();
    vonMises.resize(n);

    for (int i = 0; i < n; ++i) {
        double sxx = sigma[i][0];
        double syy = sigma[i][1];
        double szz = sigma[i][2];
        double sxy = sigma[i][3];
        double syz = sigma[i][4];
        double sxz = sigma[i][5];

        // sigma_vM = sqrt(0.5 * ((sxx-syy)^2 + (syy-szz)^2 + (szz-sxx)^2 + 6*(sxy^2 + syz^2 + sxz^2)))
        double term1 = (sxx - syy) * (sxx - syy);
        double term2 = (syy - szz) * (syy - szz);
        double term3 = (szz - sxx) * (szz - sxx);
        double term4 = 6.0 * (sxy * sxy + syz * syz + sxz * sxz);

        vonMises[i] = std::sqrt(0.5 * (term1 + term2 + term3 + term4));
    }
}

std::array<double, 3> Stress::solveCubic(double J1, double J2, double J3) {
    // sigma^3 - J1*sigma^2 + J2*sigma - J3 = 0
    // Використовуємо тригонометричний метод для трьох дійсних коренів

    // Зведення до вигляду t^3 + p*t + q = 0
    double p = J2 - J1 * J1 / 3.0;
    double q = 2.0 * J1 * J1 * J1 / 27.0 - J1 * J2 / 3.0 + J3;

    // Дискримінант
    double D = q * q / 4.0 + p * p * p / 27.0;

    std::array<double, 3> roots = {J1 / 3.0, J1 / 3.0, J1 / 3.0};

    if (D <= 0.0) {
        // Три дійсних корені
        double r = std::sqrt(-p * p * p / 27.0);
        double phi = std::acos(-q / (2.0 * r));
        
        double shift = J1 / 3.0;
        roots[0] = 2.0 * std::cbrt(r) * std::cos(phi / 3.0) + shift;
        roots[1] = 2.0 * std::cbrt(r) * std::cos((phi + 2.0 * M_PI) / 3.0) + shift;
        roots[2] = 2.0 * std::cbrt(r) * std::cos((phi + 4.0 * M_PI) / 3.0) + shift;
    } else {
        // Один дійсний корінь
        double u = std::cbrt(-q / 2.0 + std::sqrt(D));
        double v = std::cbrt(-q / 2.0 - std::sqrt(D));
        roots[0] = u + v + J1 / 3.0;
    }

    // Сортуємо за спаданням: sigma1 >= sigma2 >= sigma3
    std::sort(roots.begin(), roots.end(), std::greater<double>());
    
    return roots;
}