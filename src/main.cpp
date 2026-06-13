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
    #define NOMINMAX
    #include <windows.h>
    #undef min
    #undef max
    #pragma comment(lib, "opengl32.lib")
    #endif

int main() {
    #ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    #endif

    // Ініціалізація Viewer
    Viewer viewer;
    if (!viewer.init(1280, 800, "MCE - Analiz napruzeno-deformovanogo stanu")) {
        std::cerr << "Pomylka stvorennia vikna!\n";
        return 1;
    }

    // Ініціалізація GUI
    GUI gui;
    GLFWwindow* window = glfwGetCurrentContext();
    if (!gui.init(window)) {
        std::cerr << "Pomylka inicializacii GUI!\n";
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
            
            // Малюємо деформовану форму, якщо showDeformed
            if (state.showDeformed && state.femComputed && state.assembly) {
                // Автоматичний масштаб: щоб макс. переміщення було ~10% від розміру моделі
                double maxDisp = 0;
                for (int i = 0; i < state.mesh->nqp; ++i) {
                    double mag = state.assembly->getDisplacement(i).norm();
                    if (mag > maxDisp) maxDisp = mag;
                }
                if (maxDisp > 1e-30) {
                    double maxDim = std::max({state.mesh->a, state.mesh->b, state.mesh->c});
                    viewer.setDeformationScale(maxDim * 0.15 / maxDisp);
                }
                viewer.renderDeformed(*state.mesh, state.assembly->U);
            }
            
            // Малюємо кольорову карту напружень
            if (state.femComputed && state.vonMises) {
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