#include "Viewer.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <GLFW/glfw3.h>

#define _USE_MATH_DEFINES
#include <math.h>

// GLU функції через glfw
#include <GLFW/glfw3native.h>

// Допоміжна функція для OpenGL
static void glVertex(const Eigen::Vector3d& v) {
    glVertex3d(v.x(), v.y(), v.z());
}

// Реалізація gluPerspective через glFrustum
static void gluPerspective(double fovY, double aspect, double zNear, double zFar) {
    double f = 1.0 / tan(fovY * M_PI / 360.0);
    glFrustum(-zNear / f * aspect, zNear / f * aspect,
              -zNear / f, zNear / f,
              zNear, zFar);
}

// Реалізація gluLookAt
static void gluLookAt(double eyeX, double eyeY, double eyeZ,
                      double centerX, double centerY, double centerZ,
                      double upX, double upY, double upZ) {
    Eigen::Vector3d F(centerX - eyeX, centerY - eyeY, centerZ - eyeZ);
    F.normalize();
    Eigen::Vector3d Up(upX, upY, upZ);
    Eigen::Vector3d s = F.cross(Up);
    s.normalize();
    Eigen::Vector3d u = s.cross(F);

    GLdouble m[16] = {
         s.x(),  u.x(), -F.x(), 0.0,
         s.y(),  u.y(), -F.y(), 0.0,
         s.z(),  u.z(), -F.z(), 0.0,
         0.0,    0.0,    0.0,   1.0
    };

    glMultMatrixd(m);
    glTranslated(-eyeX, -eyeY, -eyeZ);
}

Viewer::Viewer()
    : window(nullptr), width(800), height(600),
      showNodes(true), showElements(true), showDeformed(false),
      defScale(0.1),
      azimuth(-45.0), elevation(30.0), distance(0.0),
      centerX(0.0), centerY(0.0), centerZ(0.0),
      lastMouseX(0.0), lastMouseY(0.0), mousePressed(false) {}

Viewer::~Viewer() {
    if (window) {
        glfwDestroyWindow(window);
    }
    glfwTerminate();
}

bool Viewer::init(int w, int h, const char* title) {
    width = w;
    height = h;

    if (!glfwInit()) {
        std::cerr << "Помилка ініціалізації GLFW!\n";
        return false;
    }

    // Встановлюємо версію OpenGL 3.3 (сумісність)
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);

    window = glfwCreateWindow(width, height, title, nullptr, nullptr);
    if (!window) {
        std::cerr << "Помилка створення вікна GLFW!\n";
        glfwTerminate();
        return false;
    }

    glfwMakeContextCurrent(window);

    // Встановлюємо функції зворотного виклику
    glfwSetWindowUserPointer(window, this);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetCursorPosCallback(window, cursorPosCallback);
    glfwSetScrollCallback(window, scrollCallback);
    glfwSetKeyCallback(window, keyCallback);

    // Налаштування OpenGL
    glClearColor(0.1f, 0.1f, 0.15f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Налаштування освітлення (базове)
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

    float lightPos[] = { 1.0f, 1.0f, 1.0f, 0.0f };
    float lightAmb[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    float lightDiff[] = { 0.8f, 0.8f, 0.8f, 1.0f };
    float lightSpec[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
    glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmb);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiff);
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpec);

    return true;
}

void Viewer::run(const Mesh& mesh, 
                  const std::vector<double>* displacements,
                  const std::vector<double>* stresses) {
    if (!window) return;

    // Обчислюємо центр сцени
    centerX = mesh.a / 2.0;
    centerY = mesh.b / 2.0;
    centerZ = mesh.c / 2.0;

    // Встановлюємо відстань камери на основі розміру моделі
    double maxDim = std::max({mesh.a, mesh.b, mesh.c});
    distance = maxDim * 2.5;

    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glLoadIdentity();

        // Налаштування проекції та камери
        setupProjection();
        setupCamera();

        // Малюнки
        renderMesh(mesh);

        // Якщо є переміщення, малюємо деформовану форму
        if (displacements && showDeformed) {
            renderDeformed(mesh, *displacements);
        }

        // Якщо є напруження, малюємо їх розподіл
        if (stresses) {
            renderStresses(mesh, *stresses);
        }

        // Інформація на екрані
        // (у простому варіанті - просто в консолі)
        // У майбутньому можна додати HUD

        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}

void Viewer::setupProjection() {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    double aspect = (double)width / (double)height;
    double maxDim = distance;
    double near = maxDim * 0.01;
    double far = maxDim * 100.0;
    gluPerspective(45.0, aspect, near, far);
    glMatrixMode(GL_MODELVIEW);
}

void Viewer::setupCamera() {
    double radAz = azimuth * M_PI / 180.0;
    double radEl = elevation * M_PI / 180.0;

    double x = centerX + distance * cos(radEl) * sin(radAz);
    double y = centerY + distance * cos(radEl) * cos(radAz);
    double z = centerZ + distance * sin(radEl);

    gluLookAt(x, y, z,
              centerX, centerY, centerZ,
              0.0, 0.0, 1.0);
}

void Viewer::renderMesh(const Mesh& mesh) {
    if (showElements) {
        // Малюємо ребра елементів сірим кольором
        Color edgeColor = {0.6f, 0.6f, 0.6f};
        
        for (const auto& elem : mesh.elements) {
            std::array<Eigen::Vector3d, 20> nodeCoords;
            for (int i = 0; i < 20; ++i) {
                nodeCoords[i] = mesh.nodes[elem[i]];
            }
            drawBoxElement(nodeCoords, edgeColor);
        }
    }

    if (showNodes) {
        // Малюємо вузли червоними точками
        glPointSize(5.0f);
        glDisable(GL_LIGHTING);
        glBegin(GL_POINTS);
        glColor3f(1.0f, 0.2f, 0.2f);
        for (const auto& node : mesh.nodes) {
            glVertex(node);
        }
        glEnd();
        glEnable(GL_LIGHTING);

        // Позначаємо закріплені вузли синім
        glPointSize(7.0f);
        glDisable(GL_LIGHTING);
        glBegin(GL_POINTS);
        glColor3f(0.2f, 0.2f, 1.0f);
        for (int fixedIdx : mesh.fixedNodes) {
            glVertex(mesh.nodes[fixedIdx]);
        }
        glEnd();
        glEnable(GL_LIGHTING);

        // Позначаємо навантажені вузли зеленим
        // (якщо визначені навантажені грані)
        if (!mesh.loadedFaces.empty()) {
            glPointSize(7.0f);
            glDisable(GL_LIGHTING);
            glBegin(GL_POINTS);
            glColor3f(0.2f, 1.0f, 0.2f);
            // Вузли на навантажених гранях верхньої грані
            // Для спрощення - всі вузли з z = c
            for (const auto& node : mesh.nodes) {
                if (std::abs(node.z() - mesh.c) < 1e-10) {
                    glVertex(node);
                }
            }
            glEnd();
            glEnable(GL_LIGHTING);
        }
    }
}

void Viewer::renderNodes(const Mesh& mesh) {
    // Вже включено в renderMesh
}

void Viewer::renderElements(const Mesh& mesh) {
    // Вже включено в renderMesh
}

void Viewer::renderDeformed(const Mesh& mesh, const std::vector<double>& displacements) {
    if (displacements.size() < 3 * mesh.nodes.size()) return;

    Color defColor = {1.0f, 0.8f, 0.2f}; // жовтий для деформованої форми
    glLineWidth(2.0f);
    glDisable(GL_LIGHTING);

    for (const auto& elem : mesh.elements) {
        std::array<Eigen::Vector3d, 20> nodeCoords;
        for (int i = 0; i < 20; ++i) {
            int idx = elem[i];
            Eigen::Vector3d def(
                defScale * displacements[3 * idx],
                defScale * displacements[3 * idx + 1],
                defScale * displacements[3 * idx + 2]
            );
            nodeCoords[i] = mesh.nodes[idx] + def;
        }
        drawBoxElement(nodeCoords, defColor);
    }

    glEnable(GL_LIGHTING);
    glLineWidth(1.0f);
}

void Viewer::renderStresses(const Mesh& mesh, const std::vector<double>& stresses) {
    if (stresses.size() < mesh.nodes.size()) return;

    double minStress = *std::min_element(stresses.begin(), stresses.end());
    double maxStress = *std::max_element(stresses.begin(), stresses.end());

    if (std::abs(maxStress - minStress) < 1e-10) {
        maxStress = minStress + 1.0;
    }

    // Малюємо зафарбовані грані елементів
    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);

    for (const auto& elem : mesh.elements) {
        // Збираємо координати вузлів
        std::array<Eigen::Vector3d, 20> nodeCoords;
        std::array<double, 20> elemStresses;
        for (int i = 0; i < 20; ++i) {
            nodeCoords[i] = mesh.nodes[elem[i]];
            elemStresses[i] = stresses[elem[i]];
        }

        // Малюємо 6 граней (кожна як 4 трикутники)
        // Нижня грань (вузли 0,1,2,3)
        Color c0 = stressToColor(elemStresses[0], minStress, maxStress);
        Color c1 = stressToColor(elemStresses[1], minStress, maxStress);
        Color c2 = stressToColor(elemStresses[2], minStress, maxStress);
        Color c3 = stressToColor(elemStresses[3], minStress, maxStress);
        
        glBegin(GL_TRIANGLES);
        drawQuad(nodeCoords[0], nodeCoords[1], nodeCoords[2], nodeCoords[3], c0, c1, c2, c3);
        glEnd();
    }

    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
}

Viewer::Color Viewer::stressToColor(double value, double minVal, double maxVal) {
    // Перетворення значення напруження в колір (синій -> зелений -> жовтий -> червоний)
    double t = (value - minVal) / (maxVal - minVal);
    t = std::max(0.0, std::min(1.0, t));

    Color c;
    if (t < 0.25) {
        // синій -> блакитний
        double s = t / 0.25;
        c = {0.0f, (float)s, 1.0f};
    } else if (t < 0.5) {
        // блакитний -> зелений
        double s = (t - 0.25) / 0.25;
        c = {0.0f, 1.0f, (float)(1.0 - s)};
    } else if (t < 0.75) {
        // зелений -> жовтий
        double s = (t - 0.5) / 0.25;
        c = {(float)s, 1.0f, 0.0f};
    } else {
        // жовтий -> червоний
        double s = (t - 0.75) / 0.25;
        c = {1.0f, (float)(1.0 - s), 0.0f};
    }
    return c;
}

void Viewer::drawLine(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Color& color) {
    glColor3f(color.r, color.g, color.b);
    glBegin(GL_LINES);
    glVertex(p1);
    glVertex(p2);
    glEnd();
}

void Viewer::drawQuad(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2,
                      const Eigen::Vector3d& p3, const Eigen::Vector3d& p4, 
                      const Color& c1, const Color& c2, const Color& c3, const Color& c4) {
    // Перший трикутник
    glColor3f(c1.r, c1.g, c1.b); glVertex(p1);
    glColor3f(c2.r, c2.g, c2.b); glVertex(p2);
    glColor3f(c4.r, c4.g, c4.b); glVertex(p4);
    // Другий трикутник
    glColor3f(c2.r, c2.g, c2.b); glVertex(p2);
    glColor3f(c3.r, c3.g, c3.b); glVertex(p3);
    glColor3f(c4.r, c4.g, c4.b); glVertex(p4);
}

void Viewer::drawBoxElement(const std::array<Eigen::Vector3d, 20>& nodes, const Color& color) {
    glColor3f(color.r, color.g, color.b);
    glLineWidth(1.5f);

    // Нижня грань (z = -1) - вузли 0,1,2,3
    glBegin(GL_LINE_LOOP);
    glVertex(nodes[0]); glVertex(nodes[1]);
    glVertex(nodes[2]); glVertex(nodes[3]);
    glEnd();

    // Верхня грань (z = 1) - вузли 4,5,6,7
    glBegin(GL_LINE_LOOP);
    glVertex(nodes[4]); glVertex(nodes[5]);
    glVertex(nodes[6]); glVertex(nodes[7]);
    glEnd();

    // Вертикальні ребра
    glBegin(GL_LINES);
    glVertex(nodes[0]); glVertex(nodes[4]);
    glVertex(nodes[1]); glVertex(nodes[5]);
    glVertex(nodes[2]); glVertex(nodes[6]);
    glVertex(nodes[3]); glVertex(nodes[7]);
    glEnd();

    // Серединні вузли на нижній грані
    glBegin(GL_LINES);
    glVertex(nodes[0]); glVertex(nodes[8]);  glVertex(nodes[8]);  glVertex(nodes[1]);
    glVertex(nodes[1]); glVertex(nodes[9]);  glVertex(nodes[9]);  glVertex(nodes[2]);
    glVertex(nodes[2]); glVertex(nodes[10]); glVertex(nodes[10]); glVertex(nodes[3]);
    glVertex(nodes[3]); glVertex(nodes[11]); glVertex(nodes[11]); glVertex(nodes[0]);

    // Серединні вузли на верхній грані
    glVertex(nodes[4]); glVertex(nodes[12]); glVertex(nodes[12]); glVertex(nodes[5]);
    glVertex(nodes[5]); glVertex(nodes[13]); glVertex(nodes[13]); glVertex(nodes[6]);
    glVertex(nodes[6]); glVertex(nodes[14]); glVertex(nodes[14]); glVertex(nodes[7]);
    glVertex(nodes[7]); glVertex(nodes[15]); glVertex(nodes[15]); glVertex(nodes[4]);

    // Вертикальні серединні ребра
    glVertex(nodes[0]); glVertex(nodes[16]); glVertex(nodes[16]); glVertex(nodes[4]);
    glVertex(nodes[1]); glVertex(nodes[17]); glVertex(nodes[17]); glVertex(nodes[5]);
    glVertex(nodes[2]); glVertex(nodes[18]); glVertex(nodes[18]); glVertex(nodes[6]);
    glVertex(nodes[3]); glVertex(nodes[19]); glVertex(nodes[19]); glVertex(nodes[7]);
    glEnd();
}

// ========== ОБРОБНИКИ ПОДІЙ ==========

void Viewer::mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
    auto* viewer = static_cast<Viewer*>(glfwGetWindowUserPointer(window));
    if (viewer) viewer->handleMouseButton(button, action, mods);
}

void Viewer::cursorPosCallback(GLFWwindow* window, double xpos, double ypos) {
    auto* viewer = static_cast<Viewer*>(glfwGetWindowUserPointer(window));
    if (viewer) viewer->handleCursorPos(xpos, ypos);
}

void Viewer::scrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
    auto* viewer = static_cast<Viewer*>(glfwGetWindowUserPointer(window));
    if (viewer) viewer->handleScroll(xoffset, yoffset);
}

void Viewer::keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    auto* viewer = static_cast<Viewer*>(glfwGetWindowUserPointer(window));
    if (viewer) viewer->handleKey(key, scancode, action, mods);
}

void Viewer::handleMouseButton(int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        mousePressed = (action == GLFW_PRESS);
        if (mousePressed) {
            glfwGetCursorPos(window, &lastMouseX, &lastMouseY);
        }
    }
}

void Viewer::handleCursorPos(double xpos, double ypos) {
    if (mousePressed) {
        double dx = xpos - lastMouseX;
        double dy = ypos - lastMouseY;
        azimuth += dx * 0.5;
        elevation += dy * 0.5;
        elevation = std::max(-89.0, std::min(89.0, elevation));
        lastMouseX = xpos;
        lastMouseY = ypos;
    }
}

void Viewer::handleScroll(double xoffset, double yoffset) {
    distance *= (1.0 - yoffset * 0.1);
    distance = std::max(0.1, distance);
}

void Viewer::handleKey(int key, int scancode, int action, int mods) {
    if (action != GLFW_PRESS) return;

    switch (key) {
        case GLFW_KEY_ESCAPE:
            glfwSetWindowShouldClose(window, true);
            break;
        case GLFW_KEY_N:
            showNodes = !showNodes;
            break;
        case GLFW_KEY_E:
            showElements = !showElements;
            break;
        case GLFW_KEY_D:
            showDeformed = !showDeformed;
            break;
        case GLFW_KEY_R:
            azimuth = -45.0;
            elevation = 30.0;
            break;
        case GLFW_KEY_UP:
            defScale *= 1.2;
            break;
        case GLFW_KEY_DOWN:
            defScale *= 0.8;
            break;
    }
}