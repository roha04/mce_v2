#pragma once
#include <vector>
#include <array>
#include <Eigen/Dense>

class Mesh {
public:
    // Геометричні розміри паралелепіпеда
    double a, b, c;
    
    // Кількість розбиттів по осях X, Y, Z (для сітки)
    int nx, ny, nz;

    // Загальна кількість вузлів та елементів
    int nqp, nel, ng;

    // Масив координат вузлів (аналог AKT)
    std::vector<Eigen::Vector3d> nodes;

    // Масив зв'язності елементів (аналог NT) - 20 вузлів для серендипового елемента
    std::vector<std::array<int, 20>> elements;

    // Масив закріплених вузлів (ZU) - вузли нижньої грані (z=0)
    std::vector<int> fixedNodes;

    // Масив навантажених граней (ZP) - верхня грань (z=c)
    struct LoadedFace {
        int elemId;      // номер елемента
        int faceId;      // номер грані (5 - верхня грань)
        double pressure; // величина навантаження
    };
    std::vector<LoadedFace> loadedFaces;

    // Конструктор
    Mesh(double length_a, double width_b, double height_c, 
         int divs_x, int divs_y, int divs_z);

    // Метод для генерації сітки
    void generateParallelepiped();

    // Обчислення півширини стрічки матриці
    int computeBandwidth() const;

    // Виведення базової інформації
    void printInfo() const;

private:
    // Допоміжна функція: обчислення глобального індексу вузла за його позицією в сітці
    inline int nodeIndex(int i, int j, int k) const {
        return i + j * (2 * nx + 1) + k * (2 * nx + 1) * (2 * ny + 1);
    }
};
