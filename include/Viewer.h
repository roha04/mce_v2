#pragma once
#include <GLFW/glfw3.h>
#include <Eigen/Dense>
#include "Mesh.h"
#include <functional>

class Viewer {
public:
    Viewer();
    ~Viewer();

    // Ініціалізація вікна
    bool init(int width, int height, const char* title);

    // Запуск головного циклу візуалізації
    void run(const Mesh& mesh, 
             const std::vector<double>* displacements = nullptr,
             const std::vector<double>* stresses = nullptr);

    // Перевірка, чи вікно відкрите
    bool isRunning() const { return !glfwWindowShouldClose(window); }

    // Встановлення параметрів відображення
    void setShowNodes(bool show) { showNodes = show; }
    void setShowElements(bool show) { showElements = show; }
    void setShowDeformed(bool show) { showDeformed = show; }
    void setDeformationScale(double scale) { defScale = scale; }

private:
    GLFWwindow* window;
    int width, height;

    // Параметри відображення
    bool showNodes;
    bool showElements;
    bool showDeformed;
    double defScale;

    // Стан камери
    double azimuth, elevation, distance;
    double centerX, centerY, centerZ;

    // Стан миші
    double lastMouseX, lastMouseY;
    bool mousePressed;

    // Колірні схеми для напружень
    struct Color {
        float r, g, b;
    };
    static Color stressToColor(double value, double minVal, double maxVal);

public:
    // Рендеринг (публічні для використання з GUI)
    void setupProjection();
    void setupCamera();
    void renderMesh(const Mesh& mesh);
    void renderNodes(const Mesh& mesh);
    void renderElements(const Mesh& mesh);
    void renderDeformed(const Mesh& mesh, const std::vector<double>& displacements);
    void renderStresses(const Mesh& mesh, const std::vector<double>& stresses);
    
    // Додаткові функції візуалізації
    void renderLoadArrows(const Mesh& mesh, double pressure);
    void renderNodeLabels(const Mesh& mesh);
    void renderLoadInfo(const Mesh& mesh, double pressure);

private:
    // Обробка вводу
    static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
    static void cursorPosCallback(GLFWwindow* window, double xpos, double ypos);
    static void scrollCallback(GLFWwindow* window, double xoffset, double yoffset);
    static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
    void handleMouseButton(int button, int action, int mods);
    void handleCursorPos(double xpos, double ypos);
    void handleScroll(double xoffset, double yoffset);
    void handleKey(int key, int scancode, int action, int mods);

    // Допоміжна функція для малювання лінії
    void drawLine(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Color& color);
    void drawBoxElement(const std::array<Eigen::Vector3d, 20>& nodeCoords, const Color& color);
    void drawQuad(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2,
                  const Eigen::Vector3d& p3, const Eigen::Vector3d& p4, 
                  const Color& c1, const Color& c2, const Color& c3, const Color& c4);
};