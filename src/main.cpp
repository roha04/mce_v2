#include <iostream>
#include <Eigen/Dense>
#include <GLFW/glfw3.h>
#include "../include/Mesh.h"
#include "../include/Viewer.h"
#include "../include/Gauss.h"
#include "../include/Element.h"
#include "../include/Assembly.h"
#include "../include/Stress.h"
#include "../include/GUI.h"

#ifdef _WIN32
#include <windows.h>
#pragma comment(lib, "opengl32.lib")
#endif

int main() {
    #ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    #endif

    std::cout << "=== Програма МСЕ (Паралелепіпед) ===\n\n";

    // Ініціалізація Viewer (вікно OpenGL)
    Viewer viewer;
    if (!viewer.init(1280, 800, "МСЕ - Аналіз напружено-деформованого стану")) {
        std::cerr << "Помилка створення вікна!\n";
        return 1;
    }

    // Ініціалізація GUI
    GUI gui;
    GLFWwindow* window = glfwGetCurrentContext();
    if (!gui.init(window)) {
        std::cerr << "Помилка ініціалізації GUI!\n";
        return 1;
    }

    // Стан додатку
    AppState state;

    // Головний цикл
    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glLoadIdentity();

        // Відображення 3D сцени
        viewer.setupProjection();
        viewer.setupCamera();
        
        // Малюємо сітку, якщо є
        if (state.mesh && state.meshGenerated) {
            viewer.renderMesh(*state.mesh);
            
            // Малюємо деформовану форму, якщо є переміщення
            if (state.femComputed && state.assembly) {
                viewer.renderDeformed(*state.mesh, state.assembly->U);
                viewer.renderStresses(*state.mesh, *state.vonMises);
            }
        }

        // Малюємо GUI поверх 3D
        gui.render(state);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // Очищення
    gui.shutdown();
    
    delete state.mesh;
    delete state.assembly;
    delete state.sigma;
    delete state.vonMises;
    delete state.sigma1;
    delete state.sigma2;
    delete state.sigma3;
    delete state.DFIABG;
    delete state.DPSITE;

    return 0;
}