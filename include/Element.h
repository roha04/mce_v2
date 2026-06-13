#pragma once
#include <array>
#include <Eigen/Dense>
#include "Mesh.h"
#include "Gauss.h"

class Element {
public:
    // Константа для матеріалу
    static constexpr double YOUNG_MODULUS = 2.1e11;  // E - модуль Юнга
    static constexpr double POISSON_RATIO = 0.3;     // nu - коефіцієнт Пуассона

    // Обчислення матриці жорсткості елемента MGE[60][60]
    // elemNodes - координати 20 вузлів елемента
    // DFIABG - масив похідних базисних функцій у вузлах Гауса
    static void computeStiffnessMatrix(
        const std::array<Eigen::Vector3d, 20>& elemNodes,
        const std::array<std::array<std::array<double, 20>, 3>, 27>& DFIABG,
        Eigen::Matrix<double, 60, 60>& MGE);

    // Обчислення вектора навантаження елемента FE[60]
    // for a given face (0-5) with pressure P
    // elemNodes - координати 20 вузлів елемента
    // DPSITE - масив похідних 2D базисних функцій у вузлах Гауса
    // faceLocalNodes - локальні номери вузлів на цій грані (8 вузлів)
    static void computeLoadVector(
        const std::array<Eigen::Vector3d, 20>& elemNodes,
        const std::array<std::array<std::array<double, 8>, 2>, 9>& DPSITE,
        const std::array<int, 8>& faceLocalNodes,
        double pressure,
        Eigen::Matrix<double, 60, 1>& FE);

private:
    // Допоміжні функції для обчислення коефіцієнтів a_ij (див. формули 15)
    static double computeA11(int i, int j,
                             const std::array<double, 20>& df_dx,
                             const std::array<double, 20>& df_dy,
                             const std::array<double, 20>& df_dz);

    static double computeA12(int i, int j,
                             const std::array<double, 20>& df_dx,
                             const std::array<double, 20>& df_dy);

    static double computeA13(int i, int j,
                             const std::array<double, 20>& df_dx,
                             const std::array<double, 20>& df_dz);

    static double computeA22(int i, int j,
                             const std::array<double, 20>& df_dx,
                             const std::array<double, 20>& df_dy,
                             const std::array<double, 20>& df_dz);

    static double computeA23(int i, int j,
                             const std::array<double, 20>& df_dy,
                             const std::array<double, 20>& df_dz);

    static double computeA33(int i, int j,
                             const std::array<double, 20>& df_dx,
                             const std::array<double, 20>& df_dy,
                             const std::array<double, 20>& df_dz);
};