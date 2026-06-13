#include "Viewer.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <GLFW/glfw3.h>
#include <cstdio>

#define _USE_MATH_DEFINES
#include <math.h>

// Допоміжна функція для OpenGL
static void glVertex(const Eigen::Vector3d& v) {
    glVertex3d(v.x(), v.y(), v.z());
}

// Реалізація gluPerspective
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

// Допоміжна функція для малювання сфери
static void drawSphere(double radius, int segments = 12) {
    for (int i = 0; i < segments; ++i) {
        double lat0 = M_PI * (-0.5 + (double)(i) / segments);
        double z0 = sin(lat0);
        double zr0 = cos(lat0);
        double lat1 = M_PI * (-0.5 + (double)(i + 1) / segments);
        double z1 = sin(lat1);
        double zr1 = cos(lat1);

        glBegin(GL_QUAD_STRIP);
        for (int j = 0; j <= segments; ++j) {
            double lng = 2 * M_PI * (double)(j) / segments;
            double x = cos(lng);
            double y = sin(lng);
            glNormal3d(x * zr0, y * zr0, z0);
            glVertex3d(x * zr0 * radius, y * zr0 * radius, z0 * radius);
            glNormal3d(x * zr1, y * zr1, z1);
            glVertex3d(x * zr1 * radius, y * zr1 * radius, z1 * radius);
        }
        glEnd();
    }
}

// Малювання стрілки
static void drawArrow(const Eigen::Vector3d& from, const Eigen::Vector3d& dir, 
                       double length, double headLength = 0.2, double headRadius = 0.08) {
    Eigen::Vector3d to = from + dir * length;
    Eigen::Vector3d normDir = dir.normalized();
    
    // Тіло стрілки
    glBegin(GL_LINES);
    glVertex(from);
    glVertex(to);
    glEnd();
    
    // Головка стрілки (конус)
    Eigen::Vector3d headBase = to - normDir * headLength;
    Eigen::Vector3d perp1 = normDir.cross(Eigen::Vector3d(0, 0, 1));
    if (perp1.norm() < 0.01) perp1 = normDir.cross(Eigen::Vector3d(0, 1, 0));
    perp1.normalize();
    Eigen::Vector3d perp2 = normDir.cross(perp1);
    perp2.normalize();
    
    int nSeg = 8;
    glBegin(GL_TRIANGLE_FAN);
    glVertex(to);
    for (int i = 0; i <= nSeg; ++i) {
        double angle = 2.0 * M_PI * i / nSeg;
        Eigen::Vector3d p = headBase + (perp1 * cos(angle) + perp2 * sin(angle)) * headRadius;
        glVertex(p);
    }
    glEnd();
}

Viewer::Viewer()
    : window(nullptr), width(800), height(600),
      showNodes(true), showElements(true), showDeformed(false),
      defScale(0.01),
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

    glfwSetWindowUserPointer(window, this);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetCursorPosCallback(window, cursorPosCallback);
    glfwSetScrollCallback(window, scrollCallback);
    glfwSetKeyCallback(window, keyCallback);

    glClearColor(0.12f, 0.12f, 0.17f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

    float lightPos[] = { 2.0f, 2.0f, 3.0f, 0.0f };
    float lightAmb[] = { 0.3f, 0.3f, 0.3f, 1.0f };
    float lightDiff[] = { 0.9f, 0.9f, 0.9f, 1.0f };
    float lightSpec[] = { 0.5f, 0.5f, 0.5f, 1.0f };
    glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
    glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmb);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiff);
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpec);

    return true;
}

void Viewer::renderLoadArrows(const Mesh& mesh, double pressure) {
    if (pressure == 0.0) return;
    
    glDisable(GL_LIGHTING);
    
    // Малюємо стрілки на верхній грані (z = c)
    // Створюємо сітку стрілок на поверхні
    int nx = mesh.nx;
    int ny = mesh.ny;
    
    double stepX = mesh.a / nx;
    double stepY = mesh.b / ny;
    double z = mesh.c;
    
    // Визначаємо масштаб: довжина стрілки пропорційна тиску
    double maxDim = std::max({mesh.a, mesh.b, mesh.c});
    double arrowLength = maxDim * 0.15 * (pressure / 1e6); // масштабуємо
    
    // Червоний колір для навантаження
    glColor3f(1.0f, 0.1f, 0.1f);
    
    // Малюємо стрілки в центрах елементів верхньої грані
    for (int ey = 0; ey < ny; ++ey) {
        for (int ex = 0; ex < nx; ++ex) {
            double cx = (ex + 0.5) * stepX;
            double cy = (ey + 0.5) * stepY;
            
            // Прозора червона стрілка вниз (навантаження притискає)
            drawArrow(Eigen::Vector3d(cx, cy, z + arrowLength * 0.1),
                     Eigen::Vector3d(0, 0, -1), arrowLength);
        }
    }
    
    // Малюємо пунктирну лінію по контуру верхньої грані
    glEnable(GL_LINE_STIPPLE);
    glLineStipple(1, 0x00FF);
    glColor3f(1.0f, 0.3f, 0.3f);
    glLineWidth(2.0f);
    glBegin(GL_LINE_LOOP);
    glVertex3d(0, 0, z);
    glVertex3d(mesh.a, 0, z);
    glVertex3d(mesh.a, mesh.b, z);
    glVertex3d(0, mesh.b, z);
    glEnd();
    glDisable(GL_LINE_STIPPLE);
    glLineWidth(1.0f);
    
    // Надпис "P"
    // (у спрощеному варіанті - через бібліотеку)
    
    glEnable(GL_LIGHTING);
}

void Viewer::renderNodeLabels(const Mesh& mesh) {
    if (!showNodes) return;
    
    glDisable(GL_LIGHTING);
    glPointSize(8.0f);
    
    // Малюємо сферичні вузли з кольоровим кодуванням
    for (size_t i = 0; i < mesh.nodes.size(); ++i) {
        const auto& node = mesh.nodes[i];
        
        // Визначаємо колір вузла залежно від типу
        if (std::abs(node.z() - mesh.c) < 1e-10) {
            // Верхня грань (навантажена) - червоні
            glColor3f(1.0f, 0.3f, 0.3f);
        } else if (std::abs(node.z()) < 1e-10) {
            // Нижня грань (закріплена) - сині
            glColor3f(0.3f, 0.3f, 1.0f);
        } else {
            // Внутрішні вузли - зелені
            glColor3f(0.3f, 1.0f, 0.3f);
        }
        
        // Малюємо вузол як сферу (тільки для показових вузлів)
        bool isCorner = false;
        for (int ei = 0; ei < mesh.nel; ++ei) {
            for (int ni = 0; ni < 20; ++ni) {
                if (mesh.elements[ei][ni] == (int)i && ni < 8) {
                    isCorner = true;
                    break;
                }
            }
            if (isCorner) break;
        }
        
        if (isCorner || i < 20) {
            // Кутові вузли - сфери
            glPushMatrix();
            glTranslated(node.x(), node.y(), node.z());
            drawSphere(0.08 * std::max({mesh.a, mesh.b, mesh.c}) / 10.0);
            glPopMatrix();
        } else {
            // Серединні вузли - точки
            glBegin(GL_POINTS);
            glVertex(node);
            glEnd();
        }
    }
    
    glEnable(GL_LIGHTING);
}

void Viewer::renderLoadInfo(const Mesh& mesh, double pressure) {
    // Малюємо інформаційний текст у 3D просторі
    // (через бібліотеку glut - не використовуємо,
    //  ця інформація відображається в GUI)
}

void Viewer::run(const Mesh& mesh, 
                  const std::vector<double>* displacements,
                  const std::vector<double>* stresses) {
    if (!window) return;

    centerX = mesh.a / 2.0;
    centerY = mesh.b / 2.0;
    centerZ = mesh.c / 2.0;

    double maxDim = std::max({mesh.a, mesh.b, mesh.c});
    distance = maxDim * 2.5;

    // Не запускаємо свій цикл - main.cpp керує рендерингом
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
    // Малюємо елементи сірим кольором
    if (showElements) {
        Color edgeColor = {0.5f, 0.5f, 0.6f};
        
        for (const auto& elem : mesh.elements) {
            std::array<Eigen::Vector3d, 20> nodeCoords;
            for (int i = 0; i < 20; ++i) {
                nodeCoords[i] = mesh.nodes[elem[i]];
            }
            drawBoxElement(nodeCoords, edgeColor);
        }
    }

    // Малюємо вузли
    renderNodeLabels(mesh);
    
    // Малюємо стрілки навантаження
    if (!mesh.loadedFaces.empty()) {
        double pressure = mesh.loadedFaces[0].pressure;
        renderLoadArrows(mesh, pressure);
    }
}

void Viewer::renderNodes(const Mesh& mesh) {
    renderNodeLabels(mesh);
}

void Viewer::renderElements(const Mesh& mesh) {
    if (showElements) {
        Color edgeColor = {0.5f, 0.5f, 0.6f};
        for (const auto& elem : mesh.elements) {
            std::array<Eigen::Vector3d, 20> nodeCoords;
            for (int i = 0; i < 20; ++i) {
                nodeCoords[i] = mesh.nodes[elem[i]];
            }
            drawBoxElement(nodeCoords, edgeColor);
        }
    }
}

void Viewer::renderDeformed(const Mesh& mesh, const std::vector<double>& displacements) {
    if (displacements.size() < 3 * mesh.nodes.size()) return;

    Color defColor = {1.0f, 0.8f, 0.2f};
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

    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);

    // Малюємо всі 6 граней кожного елемента з кольором напружень
    for (const auto& elem : mesh.elements) {
        std::array<Eigen::Vector3d, 20> nodeCoords;
        std::array<double, 20> elemStresses;
        for (int i = 0; i < 20; ++i) {
            nodeCoords[i] = mesh.nodes[elem[i]];
            elemStresses[i] = stresses[elem[i]];
        }

        // Визначаємо грані (кожна з 4 вузлів)
        const int faces[6][4] = {
            {0, 1, 2, 3},  // нижня
            {4, 5, 6, 7},  // верхня
            {0, 4, 7, 3},  // ліва
            {1, 5, 6, 2},  // права
            {0, 1, 5, 4},  // передня
            {2, 3, 7, 6}   // задня
        };

        for (int f = 0; f < 6; ++f) {
            Color c0 = stressToColor(elemStresses[faces[f][0]], minStress, maxStress);
            Color c1 = stressToColor(elemStresses[faces[f][1]], minStress, maxStress);
            Color c2 = stressToColor(elemStresses[faces[f][2]], minStress, maxStress);
            Color c3 = stressToColor(elemStresses[faces[f][3]], minStress, maxStress);
            
            glBegin(GL_TRIANGLES);
            drawQuad(nodeCoords[faces[f][0]], nodeCoords[faces[f][1]], 
                     nodeCoords[faces[f][2]], nodeCoords[faces[f][3]], 
                     c0, c1, c2, c3);
            glEnd();
        }
    }

    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
}

Viewer::Color Viewer::stressToColor(double value, double minVal, double maxVal) {
    double t = (value - minVal) / (maxVal - minVal);
    t = std::max(0.0, std::min(1.0, t));

    Color c;
    if (t < 0.25) {
        double s = t / 0.25;
        c = {0.0f, (float)s, 1.0f};
    } else if (t < 0.5) {
        double s = (t - 0.25) / 0.25;
        c = {0.0f, 1.0f, (float)(1.0 - s)};
    } else if (t < 0.75) {
        double s = (t - 0.5) / 0.25;
        c = {(float)s, 1.0f, 0.0f};
    } else {
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
    glColor3f(c1.r, c1.g, c1.b); glVertex(p1);
    glColor3f(c2.r, c2.g, c2.b); glVertex(p2);
    glColor3f(c4.r, c4.g, c4.b); glVertex(p4);
    
    glColor3f(c2.r, c2.g, c2.b); glVertex(p2);
    glColor3f(c3.r, c3.g, c3.b); glVertex(p3);
    glColor3f(c4.r, c4.g, c4.b); glVertex(p4);
}

void Viewer::drawBoxElement(const std::array<Eigen::Vector3d, 20>& nodes, const Color& color) {
    glColor3f(color.r, color.g, color.b);
    glLineWidth(1.5f);

    // Нижня грань
    glBegin(GL_LINE_LOOP);
    glVertex(nodes[0]); glVertex(nodes[1]);
    glVertex(nodes[2]); glVertex(nodes[3]);
    glEnd();

    // Верхня грань
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

    // Серединні ребра нижньої грані
    glBegin(GL_LINES);
    glVertex(nodes[0]); glVertex(nodes[8]);  glVertex(nodes[8]);  glVertex(nodes[1]);
    glVertex(nodes[1]); glVertex(nodes[9]);  glVertex(nodes[9]);  glVertex(nodes[2]);
    glVertex(nodes[2]); glVertex(nodes[10]); glVertex(nodes[10]); glVertex(nodes[3]);
    glVertex(nodes[3]); glVertex(nodes[11]); glVertex(nodes[11]); glVertex(nodes[0]);

    // Серединні ребра верхньої грані
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
            defScale *= 1.5;
            break;
        case GLFW_KEY_DOWN:
            defScale *= 0.67;
            break;
    }
}