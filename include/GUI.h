#pragma once
#include <GLFW/glfw3.h>
#include <vector>
#include <array>
#include <string>
#include "Mesh.h"
#include "Assembly.h"

struct AppState {
    // Параметри задачі
    double a = 10.0, b = 5.0, c = 2.0;   // геометрія
    int nx = 2, ny = 2, nz = 1;           // сітка
    double P = 100.0;                      // навантаження
    double E = 2.1e11;                     // модуль Юнга
    double nu = 0.3;                       // коеф. Пуассона

    // Стан програми
    bool meshGenerated = false;
    bool femComputed = false;
    bool showDeformed = false;
    double progress = 0.0;
    int currentStep = 0;
    std::string statusMsg = "Готовий до роботи";

    // Дані
    Mesh* mesh = nullptr;
    Assembly* assembly = nullptr;
    std::vector<std::array<double, 6>>* sigma = nullptr;
    std::vector<double>* vonMises = nullptr;
    std::vector<double>* sigma1 = nullptr;
    std::vector<double>* sigma2 = nullptr;
    std::vector<double>* sigma3 = nullptr;
    std::array<std::array<std::array<double, 20>, 3>, 27>* DFIABG = nullptr;
    std::array<std::array<std::array<double, 8>, 2>, 9>* DPSITE = nullptr;
};

class GUI {
public:
    GUI();
    ~GUI();

    bool init(GLFWwindow* window, const char* glslVersion = "#version 130");
    void render(AppState& state);
    void shutdown();

private:
    void renderInputPanel(AppState& state);
    void renderMeshPanel(AppState& state);
    void renderResultsPanel(AppState& state);
    void renderStressPanel(AppState& state);
    void renderControlsPanel(AppState& state);

    // Вкладки
    int activeTab = 0;
    int activeResultTab = 0;

    // Таблиці
    char filterBuffer[256] = "";
    int sortColumn = -1;
    bool sortAscending = true;
};