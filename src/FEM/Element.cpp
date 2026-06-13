#include "Element.h"
#include <iostream>
#include <cmath>

// Коефіцієнти Ляме
// lambda = E * nu / ((1 + nu) * (1 - 2*nu))
// mu = E / (2 * (1 + nu))
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

void Element::computeStiffnessMatrix(
    const std::array<Eigen::Vector3d, 20>& elemNodes,
    const std::array<std::array<std::array<double, 20>, 3>, 27>& DFIABG,
    Eigen::Matrix<double, 60, 60>& MGE) {
    
    MGE.setZero();
    double lam = lambda();
    double m = mu();

    // Цикл по 27 вузлах Гауса
    for (int g = 0; g < 27; ++g) {
        // Отримуємо похідні базисних функцій за локальними координатами
        // DFIABG[g][coord][func]: coord: 0=xi, 1=eta, 2=zeta
        const auto& dN_dxi = DFIABG[g][0];
        const auto& dN_deta = DFIABG[g][1];
        const auto& dN_dzeta = DFIABG[g][2];

        // 1. Формуємо якобіан DXYZABG (3x3)
        Eigen::Matrix3d J;
        J.setZero();
        for (int i = 0; i < 20; ++i) {
            double x = elemNodes[i].x();
            double y = elemNodes[i].y();
            double z = elemNodes[i].z();

            // dx/dxi, dy/dxi, dz/dxi
            J(0, 0) += dN_dxi[i] * x;
            J(1, 0) += dN_dxi[i] * y;
            J(2, 0) += dN_dxi[i] * z;

            // dx/deta, dy/deta, dz/deta
            J(0, 1) += dN_deta[i] * x;
            J(1, 1) += dN_deta[i] * y;
            J(2, 1) += dN_deta[i] * z;

            // dx/dzeta, dy/dzeta, dz/dzeta
            J(0, 2) += dN_dzeta[i] * x;
            J(1, 2) += dN_dzeta[i] * y;
            J(2, 2) += dN_dzeta[i] * z;
        }

        // 2. Обчислюємо детермінант якобіана DJ
        double detJ = J.determinant();
        if (detJ <= 0.0) {
            std::cerr << "Попередження: det(J) = " << detJ << " <= 0 у вузлі Гауса " << g << "\n";
            detJ = std::abs(detJ);
        }

        // 3. Обчислюємо обернену матрицю якобіана
        Eigen::Matrix3d Jinv = J.inverse();

        // 4. Формуємо DFIXYZ[g][func][coord] - похідні за глобальними координатами
        // dN/dx, dN/dy, dN/dz для кожної з 20 функцій
        std::array<double, 20> df_dx, df_dy, df_dz;
        for (int i = 0; i < 20; ++i) {
            df_dx[i] = Jinv(0, 0) * dN_dxi[i] + Jinv(0, 1) * dN_deta[i] + Jinv(0, 2) * dN_dzeta[i];
            df_dy[i] = Jinv(1, 0) * dN_dxi[i] + Jinv(1, 1) * dN_deta[i] + Jinv(1, 2) * dN_dzeta[i];
            df_dz[i] = Jinv(2, 0) * dN_dxi[i] + Jinv(2, 1) * dN_deta[i] + Jinv(2, 2) * dN_dzeta[i];
        }

        // Вага Гауса для цієї точки
        int ii = g % 3;
        int jj = (g / 3) % 3;
        int kk = g / 9;
        double w = Gauss::weights[ii] * Gauss::weights[jj] * Gauss::weights[kk];

        // 5. Формуємо блоки матриці жорсткості
        // Матриця MGE має 3 блоки 20x20 для (x,y,z) компонент
        // Блок (0,0): компоненти u_x - u_x
        // Блок (0,1): компоненти u_x - u_y
        // І т.д.
        
        for (int i = 0; i < 20; ++i) {
            for (int j = i; j < 20; ++j) { // симетрія - обчислюємо тільки верхній трикутник
                // Обчислення коефіцієнтів a_ij за формулами (15)
                double a11 = computeA11(i, j, df_dx, df_dy, df_dz);
                double a12 = computeA12(i, j, df_dx, df_dy);
                double a13 = computeA13(i, j, df_dx, df_dz);
                double a22 = computeA22(i, j, df_dx, df_dy, df_dz);
                double a23 = computeA23(i, j, df_dy, df_dz);
                double a33 = computeA33(i, j, df_dx, df_dy, df_dz);

                double coeff = a11; // = a11 * detJ * w
                MGE(i, j)       += a11 * detJ * w;  // ux-ux
                MGE(i, j + 20)  += a12 * detJ * w;  // ux-uy
                MGE(i, j + 40)  += a13 * detJ * w;  // ux-uz
                
                MGE(i + 20, j)      += a12 * detJ * w;  // uy-ux
                MGE(i + 20, j + 20) += a22 * detJ * w;  // uy-uy
                MGE(i + 20, j + 40) += a23 * detJ * w;  // uy-uz
                
                MGE(i + 40, j)      += a13 * detJ * w;  // uz-ux
                MGE(i + 40, j + 20) += a23 * detJ * w;  // uz-uy
                MGE(i + 40, j + 40) += a33 * detJ * w;  // uz-uz
            }
        }
    }

    // Заповнюємо нижній трикутник (симетрія)
    for (int i = 0; i < 60; ++i) {
        for (int j = i + 1; j < 60; ++j) {
            MGE(j, i) = MGE(i, j);
        }
    }
}

void Element::computeLoadVector(
    const std::array<Eigen::Vector3d, 20>& elemNodes,
    const std::array<std::array<std::array<double, 8>, 2>, 9>& DPSITE,
    const std::array<int, 8>& faceLocalNodes,
    double pressure,
    Eigen::Matrix<double, 60, 1>& FE) {
    
    FE.setZero();

    // Цикл по 9 вузлах Гауса на квадраті
    for (int g = 0; g < 9; ++g) {
        // Похідні 2D базисних функцій за xi, eta
        const auto& dN_dxi = DPSITE[g][0];   // 8 значень
        const auto& dN_deta = DPSITE[g][1];  // 8 значень

        // Обчислення dx/dxi, dy/dxi, dz/dxi та dx/deta, dy/deta, dz/deta
        // через глобальні координати вузлів грані
        Eigen::Vector3d dx_dxi(0, 0, 0), dx_deta(0, 0, 0);
        for (int i = 0; i < 8; ++i) {
            int localNode = faceLocalNodes[i];
            dx_dxi += dN_dxi[i] * elemNodes[localNode];
            dx_deta += dN_deta[i] * elemNodes[localNode];
        }

        // Обчислення нормалі: cross(dx_dxi, dx_deta)
        Eigen::Vector3d normal = dx_dxi.cross(dx_deta);
        double L = normal.norm(); // якобіан переходу

        if (L < 1e-15) {
            std::cerr << "Попередження: L = " << L << " дуже малий\n";
            continue;
        }

        // Одинична нормаль
        Eigen::Vector3d n = normal / L;

        // Навантаження: P * n
        Eigen::Vector3d load = pressure * n;

        // Вага Гауса
        int ii = g % 3;
        int jj = g / 3;
        double w = Gauss::weights[ii] * Gauss::weights[jj];

        // Значення 2D базисних функцій у вузлі Гауса
        std::array<double, 8> N;
        double xi = Gauss::points[ii];
        double eta = Gauss::points[jj];
        Gauss::shapeFunctions2D(xi, eta, N);

        // Додаємо внесок у вектор навантаження
        for (int i = 0; i < 8; ++i) {
            int localNode = faceLocalNodes[i];
            
            FE(localNode)       += N[i] * load.x() * L * w;
            FE(localNode + 20)  += N[i] * load.y() * L * w;
            FE(localNode + 40)  += N[i] * load.z() * L * w;
        }
    }
}

// ========== Реалізації a_ij ==========

double Element::computeA11(int i, int j,
                           const std::array<double, 20>& df_dx,
                           const std::array<double, 20>& df_dy,
                           const std::array<double, 20>& df_dz) {
    double lam = lambda();
    double m = mu();
    // a11 = (lambda + 2*mu) * dNi/dx * dNj/dx + mu * (dNi/dy * dNj/dy + dNi/dz * dNj/dz)
    return (lam + 2.0 * m) * df_dx[i] * df_dx[j] 
           + m * (df_dy[i] * df_dy[j] + df_dz[i] * df_dz[j]);
}

double Element::computeA12(int i, int j,
                           const std::array<double, 20>& df_dx,
                           const std::array<double, 20>& df_dy) {
    double lam = lambda();
    double m = mu();
    // a12 = lambda * dNi/dx * dNj/dy + mu * dNi/dy * dNj/dx
    return lam * df_dx[i] * df_dy[j] + m * df_dy[i] * df_dx[j];
}

double Element::computeA13(int i, int j,
                           const std::array<double, 20>& df_dx,
                           const std::array<double, 20>& df_dz) {
    double lam = lambda();
    double m = mu();
    // a13 = lambda * dNi/dx * dNj/dz + mu * dNi/dz * dNj/dx
    return lam * df_dx[i] * df_dz[j] + m * df_dz[i] * df_dx[j];
}

double Element::computeA22(int i, int j,
                           const std::array<double, 20>& df_dx,
                           const std::array<double, 20>& df_dy,
                           const std::array<double, 20>& df_dz) {
    double lam = lambda();
    double m = mu();
    // a22 = mu * dNi/dx * dNj/dx + (lambda + 2*mu) * dNi/dy * dNj/dy + mu * dNi/dz * dNj/dz
    return m * df_dx[i] * df_dx[j] 
           + (lam + 2.0 * m) * df_dy[i] * df_dy[j] 
           + m * df_dz[i] * df_dz[j];
}

double Element::computeA23(int i, int j,
                           const std::array<double, 20>& df_dy,
                           const std::array<double, 20>& df_dz) {
    double lam = lambda();
    double m = mu();
    // a23 = lambda * dNi/dy * dNj/dz + mu * dNi/dz * dNj/dy
    return lam * df_dy[i] * df_dz[j] + m * df_dz[i] * df_dy[j];
}

double Element::computeA33(int i, int j,
                           const std::array<double, 20>& df_dx,
                           const std::array<double, 20>& df_dy,
                           const std::array<double, 20>& df_dz) {
    double lam = lambda();
    double m = mu();
    // a33 = mu * (dNi/dx * dNj/dx + dNi/dy * dNj/dy) + (lambda + 2*mu) * dNi/dz * dNj/dz
    return m * (df_dx[i] * df_dx[j] + df_dy[i] * df_dy[j]) 
           + (lam + 2.0 * m) * df_dz[i] * df_dz[j];
}