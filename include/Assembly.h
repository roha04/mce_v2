#pragma once
#include <vector>
#include <Eigen/Dense>
#include "Mesh.h"
#include "Gauss.h"

class Assembly {
public:
    // Глобальна матриця жорсткості (стрічкове зберігання)
    // MG[3*nqp][ng] - лише наддіагональна частина
    std::vector<std::vector<double>> MG;
    
    // Глобальний вектор навантаження
    std::vector<double> F;
    
    // Розв'язок (вектор переміщень)
    std::vector<double> U;

    // Параметри
    int nqp;  // кількість вузлів
    int ng;   // півширина стрічки
    int systemSize; // 3 * nqp

    Assembly();
    
    // Ініціалізація масивів
    void initialize(int numNodes, int bandwidth);
    
    // Збірка глобальної системи з елементів
    void assemble(const Mesh& mesh,
                  const std::array<std::array<std::array<double, 20>, 3>, 27>& DFIABG,
                  const std::array<std::array<std::array<double, 8>, 2>, 9>& DPSITE);
    
    // Врахування крайових умов (метод штрафу)
    void applyBoundaryConditions(const std::vector<int>& fixedNodes, double penaltyFactor = 1e30);
    
    // Розв'язування СЛАР (метод Гауса для стрічкової матриці)
    void solve();
    
    // Отримання переміщення у вузлі
    Eigen::Vector3d getDisplacement(int nodeIdx) const;
    
    // Виведення інформації
    void printInfo() const;
};