#include "Assembly.h"
#include "Element.h"
#include <iostream>
#include <algorithm>
#include <cmath>

Assembly::Assembly() : nqp(0), ng(0), systemSize(0) {}

void Assembly::initialize(int numNodes, int bandwidth) {
    nqp = numNodes;
    ng = bandwidth;
    systemSize = 3 * nqp;
    
    MG.resize(systemSize);
    for (auto& row : MG) {
        row.resize(ng, 0.0);
    }
    
    F.resize(systemSize, 0.0);
    U.resize(systemSize, 0.0);
}

void Assembly::assemble(const Mesh& mesh,
                         const std::array<std::array<std::array<double, 20>, 3>, 27>& DFIABG,
                         const std::array<std::array<std::array<double, 8>, 2>, 9>& DPSITE) {
    
    std::cout << "Збірка глобальної системи...\n";
    
    // Попереднє обчислення DFIABG (похідних у вузлах Гауса)
    // DFIABG вже має бути заповнений
    
    // Тимчасові масиви для елемента
    std::array<Eigen::Vector3d, 20> elemNodes;
    Eigen::Matrix<double, 60, 60> MGE;
    Eigen::Matrix<double, 60, 1> FE;
    
    int totalElements = mesh.nel;
    
    // Цикл по елементах
    for (int e = 0; e < totalElements; ++e) {
        // Отримуємо координати вузлів елемента
        for (int i = 0; i < 20; ++i) {
            int globalNode = mesh.elements[e][i];
            elemNodes[i] = mesh.nodes[globalNode];
        }
        
        // Обчислюємо матрицю жорсткості елемента
        Element::computeStiffnessMatrix(elemNodes, DFIABG, MGE);
        
        // Обчислюємо вектор навантаження елемента (якщо є навантажені грані)
        FE.setZero();
        for (const auto& face : mesh.loadedFaces) {
            if (face.elemId == e) {
                Eigen::Matrix<double, 60, 1> FE_face;
                Element::computeLoadVector(elemNodes, DPSITE, 
                                           Gauss::faceNodes[face.faceId],
                                           face.pressure, FE_face);
                FE += FE_face;
            }
        }
        
        // Розсилання в глобальну матрицю та вектор
        for (int i = 0; i < 20; ++i) {
            int globalNode_i = mesh.elements[e][i];
            
            for (int compI = 0; compI < 3; ++compI) {
                int row = 3 * globalNode_i + compI;
                
                // Внесок у вектор F
                F[row] += FE(i + compI * 20);
                
                for (int j = 0; j < 20; ++j) {
                    int globalNode_j = mesh.elements[e][j];
                    
                    for (int compJ = 0; compJ < 3; ++compJ) {
                        int col = 3 * globalNode_j + compJ;
                        
                        double value = MGE(i + compI * 20, j + compJ * 20);
                        
                        // Стрічкове зберігання: зберігаємо тільки j >= i
                        if (col >= row) {
                            int bandIdx = col - row;
                            if (bandIdx < ng) {
                                MG[row][bandIdx] += value;
                            } else {
                                std::cerr << "Помилка: вихід за межі стрічки! row=" 
                                          << row << ", col=" << col << ", band=" << bandIdx << "\n";
                            }
                        }
                    }
                }
            }
        }
        
        // Індикатор прогресу
        if ((e + 1) % 100 == 0 || e == 0 || e == totalElements - 1) {
            std::cout << "  Оброблено елементів: " << (e + 1) << " / " << totalElements << "\r";
            std::cout.flush();
        }
    }
    std::cout << "\nЗбірка завершена.\n";
}

void Assembly::applyBoundaryConditions(const std::vector<int>& fixedNodes, double penaltyFactor) {
    std::cout << "Врахування крайових умов...\n";
    
    for (int nodeIdx : fixedNodes) {
        for (int comp = 0; comp < 3; ++comp) {
            int row = 3 * nodeIdx + comp;
            
            // Метод штрафу: встановлюємо діагональний елемент великим
            if (0 < ng) { // діагональний елемент завжди в band 0
                MG[row][0] = penaltyFactor;
            }
            
            // Встановлюємо відповідну компоненту вектора F в 0
            // (оскільки u=0 на закріпленій грані)
            F[row] = 0.0;
        }
    }
    
    std::cout << "  Застосовано до " << fixedNodes.size() << " закріплених вузлів.\n";
}

void Assembly::solve() {
    std::cout << "Розв'язування СЛАР...\n";
    
    int n = systemSize;
    int bandwidth = ng;
    
    // Спочатку виправляємо нульові діагональні елементи (нев'язані вузли)
    for (int i = 0; i < n; ++i) {
        if (std::abs(MG[i][0]) < 1e-30) {
            MG[i][0] = 1.0;  // Встановлюємо діагональ = 1
            F[i] = 0.0;      // І переміщення = 0
        }
    }
    
    // Прямий хід методу Гауса для стрічкової матриці
    for (int i = 0; i < n; ++i) {
        // Ділення рядка i на діагональний елемент
        double diag = MG[i][0];
        if (std::abs(diag) < 1e-30) {
            std::cerr << "Помилка: нульовий діагональний елемент у рядку " << i << "\n";
            continue;
        }
        
        F[i] /= diag;
        
        // Нормалізація рядка матриці
        int maxJ = std::min(bandwidth, n - i);
        for (int j = 1; j < maxJ; ++j) {
            MG[i][j] /= diag;
        }
        
        // Виключення Гауса для рядків нижче
        for (int k = 1; k < maxJ; ++k) {
            double factor = MG[i][k];
            if (std::abs(factor) < 1e-30) continue;
            
            int row = i + k;
            F[row] -= factor * F[i];
            
            int maxL = std::min(bandwidth - k, n - row);
            for (int l = 0; l < maxL; ++l) {
                MG[row][l] -= factor * MG[i][l + k];
            }
        }
    }
    
    // Зворотній хід
    U.assign(n, 0.0);
    U[n - 1] = F[n - 1];
    
    for (int i = n - 2; i >= 0; --i) {
        double sum = 0.0;
        int maxJ = std::min(bandwidth, n - i);
        for (int j = 1; j < maxJ; ++j) {
            sum += MG[i][j] * U[i + j];
        }
        U[i] = F[i] - sum;
    }
    
    std::cout << "СЛАР розв'язана.\n";
}

Eigen::Vector3d Assembly::getDisplacement(int nodeIdx) const {
    if (nodeIdx < 0 || nodeIdx >= nqp) {
        return Eigen::Vector3d::Zero();
    }
    return Eigen::Vector3d(
        U[3 * nodeIdx],
        U[3 * nodeIdx + 1],
        U[3 * nodeIdx + 2]
    );
}

void Assembly::printInfo() const {
    std::cout << "--- Система рівнянь ---\n";
    std::cout << "Розмір системи: " << systemSize << "\n";
    std::cout << "Півширина стрічки: " << ng << "\n";
    std::cout << "Розмір MG: " << systemSize << " x " << ng << "\n";
    std::cout << "------------------------\n";
}