#pragma once
#include <vector>
#include <Eigen/Dense>
#include "Mesh.h"
#include "Gauss.h"
#include "Assembly.h"

class Stress {
public:
    // Обчислення компонент тензора напружень у вузлах
    // sigma[6][nqp] = {sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_yz, sigma_xz}
    static void computeStresses(const Mesh& mesh, 
                                 const std::vector<double>& displacements,
                                 const std::array<std::array<std::array<double, 20>, 3>, 27>& DFIABG,
                                 std::vector<std::array<double, 6>>& sigma);

    // Обчислення головних напружень
    static void computePrincipalStresses(const std::vector<std::array<double, 6>>& sigma,
                                          std::vector<double>& sigma1,
                                          std::vector<double>& sigma2,
                                          std::vector<double>& sigma3);

    // Обчислення еквівалентного напруження за Мізесом
    static void computeVonMises(const std::vector<std::array<double, 6>>& sigma,
                                 std::vector<double>& vonMises);

private:
    // Обчислення компонент тензора напружень з переміщень в одній точці
    static std::array<double, 6> computeStressAtPoint(
        const std::array<double, 20>& df_dx,
        const std::array<double, 20>& df_dy,
        const std::array<double, 20>& df_dz,
        const Eigen::Matrix<double, 60, 1>& elemDisp);

    // Розв'язання кубічного рівняння для головних напружень
    // sigma^3 - J1*sigma^2 + J2*sigma - J3 = 0
    static std::array<double, 3> solveCubic(double J1, double J2, double J3);
};