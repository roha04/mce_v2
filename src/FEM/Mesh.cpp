#include "Mesh.h"
#include <iostream>
#include <cmath>
#include <algorithm>

Mesh::Mesh(double length_a, double width_b, double height_c, 
           int divs_x, int divs_y, int divs_z)
    : a(length_a), b(width_b), c(height_c), nx(divs_x), ny(divs_y), nz(divs_z),
      nqp(0), nel(0), ng(0) {}

void Mesh::generateParallelepiped() {
    nodes.clear();
    elements.clear();
    fixedNodes.clear();
    loadedFaces.clear();

    // Загальна кількість вузлів:
    // По кожній осі: 2*n + 1 точок (кутові + серединні)
    int nx2 = 2 * nx + 1;
    int ny2 = 2 * ny + 1;
    int nz2 = 2 * nz + 1;
    nqp = nx2 * ny2 * nz2;
    nel = nx * ny * nz;

    // Кроки між вузлами в глобальних координатах
    double dx = a / (2.0 * nx);
    double dy = b / (2.0 * ny);
    double dz = c / (2.0 * nz);

    // 1. Генерація всіх вузлів сітки
    nodes.resize(nqp);
    for (int k = 0; k < nz2; ++k) {
        for (int j = 0; j < ny2; ++j) {
            for (int i = 0; i < nx2; ++i) {
                int idx = nodeIndex(i, j, k);
                nodes[idx] = Eigen::Vector3d(i * dx, j * dy, k * dz);
            }
        }
    }

    // 2. Генерація елементів (20-вузлові серендипові елементи)
    elements.resize(nel);
    for (int ez = 0; ez < nz; ++ez) {
        for (int ey = 0; ey < ny; ++ey) {
            for (int ex = 0; ex < nx; ++ex) {
                int elemIdx = ex + ey * nx + ez * nx * ny;
                auto& elem = elements[elemIdx];

                // Базові індекси для кутових точок
                int i0 = 2 * ex;
                int j0 = 2 * ey;
                int k0 = 2 * ez;
                int i2 = i0 + 2;
                int j2 = j0 + 2;
                int k2 = k0 + 2;

                // --- Кутові вузли (1-8) ---
                elem[0]  = nodeIndex(i0, j0, k0); // ( -1, -1, -1)
                elem[1]  = nodeIndex(i2, j0, k0); // (  1, -1, -1)
                elem[2]  = nodeIndex(i2, j2, k0); // (  1,  1, -1)
                elem[3]  = nodeIndex(i0, j2, k0); // ( -1,  1, -1)
                elem[4]  = nodeIndex(i0, j0, k2); // ( -1, -1,  1)
                elem[5]  = nodeIndex(i2, j0, k2); // (  1, -1,  1)
                elem[6]  = nodeIndex(i2, j2, k2); // (  1,  1,  1)
                elem[7]  = nodeIndex(i0, j2, k2); // ( -1,  1,  1)

                // --- Серединні вузли на ребрах нижньої грані (z = -1) (9-12) ---
                elem[8]  = nodeIndex(i0+1, j0,   k0); // (  0, -1, -1) ребро 1-2
                elem[9]  = nodeIndex(i2,   j0+1, k0); // (  1,  0, -1) ребро 2-3
                elem[10] = nodeIndex(i0+1, j2,   k0); // (  0,  1, -1) ребро 3-4
                elem[11] = nodeIndex(i0,   j0+1, k0); // ( -1,  0, -1) ребро 4-1

                // --- Серединні вузли на ребрах верхньої грані (z = 1) (13-16) ---
                elem[12] = nodeIndex(i0+1, j0,   k2); // (  0, -1,  1) ребро 5-6
                elem[13] = nodeIndex(i2,   j0+1, k2); // (  1,  0,  1) ребро 6-7
                elem[14] = nodeIndex(i0+1, j2,   k2); // (  0,  1,  1) ребро 7-8
                elem[15] = nodeIndex(i0,   j0+1, k2); // ( -1,  0,  1) ребро 8-5

                // --- Серединні вузли на вертикальних ребрах (17-20) ---
                elem[16] = nodeIndex(i0,   j0,   k0+1); // ( -1, -1,  0) ребро 1-5
                elem[17] = nodeIndex(i2,   j0,   k0+1); // (  1, -1,  0) ребро 2-6
                elem[18] = nodeIndex(i2,   j2,   k0+1); // (  1,  1,  0) ребро 3-7
                elem[19] = nodeIndex(i0,   j2,   k0+1); // ( -1,  1,  0) ребро 4-8
            }
        }
    }

    // 3. Формування масиву закріплених вузлів (нижня грань z = 0)
    //    Всі вузли з k = 0 (перший шар вузлів)
    fixedNodes.clear();
    for (int j = 0; j < ny2; ++j) {
        for (int i = 0; i < nx2; ++i) {
            fixedNodes.push_back(nodeIndex(i, j, 0));
        }
    }

    // 4. Формування масиву навантажених граней (верхня грань z = c)
    //    Для кожного елемента з ez = nz-1, грань 5 (верхня)
    loadedFaces.clear();
    for (int ey = 0; ey < ny; ++ey) {
        for (int ex = 0; ex < nx; ++ex) {
            int elemIdx = ex + ey * nx + (nz - 1) * nx * ny;
            loadedFaces.push_back({elemIdx, 5, 0.0}); // pressure буде встановлено пізніше
        }
    }

    // 5. Обчислення півширини стрічки
    ng = computeBandwidth();

    std::cout << "Сітка успішно згенерована!\n";
    std::cout << "  Вузлів (nqp): " << nqp << "\n";
    std::cout << "  Елементів (nel): " << nel << "\n";
    std::cout << "  Півширина стрічки (ng): " << ng << "\n";
}

int Mesh::computeBandwidth() const {
    int maxBand = 0;
    for (const auto& elem : elements) {
        int minNode = *std::min_element(elem.begin(), elem.end());
        int maxNode = *std::max_element(elem.begin(), elem.end());
        int band = 3 * (maxNode - minNode + 1);
        if (band > maxBand) {
            maxBand = band;
        }
    }
    return maxBand;
}

void Mesh::printInfo() const {
    std::cout << "--- Інформація про сітку ---\n";
    std::cout << "Розміри: " << a << " x " << b << " x " << c << "\n";
    std::cout << "Розбиття: " << nx << " x " << ny << " x " << nz << "\n";
    std::cout << "Кількість вузлів (nqp): " << nodes.size() << "\n";
    std::cout << "Кількість елементів (nel): " << elements.size() << "\n";
    std::cout << "Півширина стрічки (ng): " << ng << "\n";
    std::cout << "Закріплених вузлів: " << fixedNodes.size() << "\n";
    std::cout << "Навантажених граней: " << loadedFaces.size() << "\n";
    std::cout << "----------------------------\n";

    // Виведення перших кількох вузлів для перевірки
    std::cout << "\nПерші 10 вузлів:\n";
    int showCount = std::min(10, (int)nodes.size());
    for (int i = 0; i < showCount; ++i) {
        std::cout << "  Вузол " << i << ": ("
                  << nodes[i].x() << ", " << nodes[i].y() << ", " << nodes[i].z() << ")\n";
    }

    // Виведення першого елемента для перевірки
    if (!elements.empty()) {
        std::cout << "\nПерший елемент (глобальні номери вузлів):\n  ";
        for (int i = 0; i < 20; ++i) {
            std::cout << elements[0][i] << " ";
        }
        std::cout << "\n";
    }
}