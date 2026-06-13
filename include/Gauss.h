#pragma once
#include <array>
#include <Eigen/Dense>

// Кількість вузлів Гауса
const int NGAUSS = 3;
const int NGAUSS2 = 9;  // 3x3 для поверхневих інтегралів
const int NGAUSS3 = 27; // 3x3x3 для об'ємних інтегралів

// Структура для зберігання точки Гауса
struct GaussPoint {
    double xi;   // координата в локальній системі
    double weight; // вага
};

// Клас для обчислення квадратур Гауса та базисних функцій
class Gauss {
public:
    // Координати та ваги Гауса для 3 точок
    static const std::array<double, 3> points;
    static const std::array<double, 3> weights;

    // Обчислення 20 базисних функцій серендипового елемента
    // в точці (xi, eta, zeta) в локальних координатах [-1, 1]
    static void shapeFunctions(double xi, double eta, double zeta, 
                               std::array<double, 20>& N);

    // Похідні базисних функцій за локальними координатами
    static void shapeFunctionDerivatives(double xi, double eta, double zeta,
                                          std::array<double, 20>& dN_dxi,
                                          std::array<double, 20>& dN_deta,
                                          std::array<double, 20>& dN_dzeta);

    // Обчислення масиву DFIABG[27][3][20] - похідних базисних функцій
    // у всіх вузлах Гауса
    static void computeDFIABG(std::array<std::array<std::array<double, 20>, 3>, 27>& DFIABG);

    // --- Функції для 2D (грані елемента) ---

    // 8 базисних функцій для квадрата [-1,1]x[-1,1] (грань елемента)
    static void shapeFunctions2D(double xi, double eta,
                                  std::array<double, 8>& N);

    // Похідні 2D базисних функцій
    static void shapeFunctionDerivatives2D(double xi, double eta,
                                            std::array<double, 8>& dN_dxi,
                                            std::array<double, 8>& dN_deta);

    // Обчислення масиву DPSITE[9][2][8]
    static void computeDPSITE(std::array<std::array<std::array<double, 8>, 2>, 9>& DPSITE);

    // Номери локальних вузлів на кожній з 6 граней 20-вузлового елемента
    // (для інтегрування навантаження)
    static const std::array<std::array<int, 8>, 6> faceNodes;
};