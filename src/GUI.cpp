#include "GUI.h"
#include "Gauss.h"
#include "Element.h"
#include "Stress.h"
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <iomanip>

GUI::GUI() : activeTab(0), activeResultTab(0), sortColumn(-1), sortAscending(true) {
    memset(filterBuffer, 0, sizeof(filterBuffer));
}

GUI::~GUI() {}

bool GUI::init(GLFWwindow* window, const char* glslVersion) {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.IniFilename = nullptr;

    ImGui::StyleColorsDark();

    // Try to load system font with Cyrillic support
    // First try common Windows font paths
    const char* fontPaths[] = {
        "C:\\Windows\\Fonts\\Arial.ttf",
        "C:\\Windows\\Fonts\\arial.ttf",
        "C:\\Windows\\Fonts\\SegoeUI.ttf",
        "C:\\Windows\\Fonts\\segoeui.ttf",
        "C:\\Windows\\Fonts\\Tahoma.ttf",
        "C:\\Windows\\Fonts\\tahoma.ttf",
        "C:\\Windows\\Fonts\\Verdana.ttf",
        "C:\\Windows\\Fonts\\verdana.ttf",
        nullptr
    };

    // Cyrillic + Latin character ranges
    static const ImWchar cryllicRanges[] = {
        0x0020, 0x00FF,  // Basic Latin
        0x0100, 0x017F,  // Latin Extended-A
        0x0400, 0x04FF,  // Cyrillic
        0x0500, 0x052F,  // Cyrillic Supplement
        0x20AC, 0x20AC,  // Euro
        0
    };

    ImFont* font = nullptr;
    for (int i = 0; fontPaths[i] != nullptr; i++) {
        font = io.Fonts->AddFontFromFileTTF(fontPaths[i], 15.0f, nullptr, cryllicRanges);
        if (font) {
            break;
        }
    }
    
    if (!font) {
        // Fallback to default font (Latin only)
        font = io.Fonts->AddFontDefault();
    }
    
    if (font) {
        io.FontDefault = font;
    }

    if (!ImGui_ImplGlfw_InitForOpenGL(window, true) ||
        !ImGui_ImplOpenGL3_Init(glslVersion)) {
        return false;
    }
    return true;
}

void GUI::shutdown() {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}

void GUI::render(AppState& state) {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    // Головне меню
    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("Файл")) {
            if (ImGui::MenuItem("Вийти", "ESC")) {
                glfwSetWindowShouldClose(glfwGetCurrentContext(), true);
            }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Довідка")) {
            ImGui::Text("МСЕ - аналіз паралелепіпеда");
            ImGui::Text("20-вузлові серендипові елементи");
            ImGui::Separator();
            ImGui::Text("Управління:");
            ImGui::BulletText("ЛКМ+миша: обертання");
            ImGui::BulletText("Scroll: масштабування");
            ImGui::BulletText("ПКМ+миша: зміщення");
            ImGui::EndMenu();
        }
        ImGui::EndMainMenuBar();
    }

    // Бічна панель (ліва)
    ImGui::SetNextWindowPos(ImVec2(0, 20));
    ImGui::SetNextWindowSize(ImVec2(420, ImGui::GetIO().DisplaySize.y - 20));
    ImGui::Begin("Панель керування", nullptr, 
                 ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);

    // Вкладки (використовуємо radio buttons)
    ImGui::RadioButton("Введення даних", &activeTab, 0); ImGui::SameLine();
    ImGui::RadioButton("Сітка", &activeTab, 1); ImGui::SameLine();
    ImGui::RadioButton("Результати", &activeTab, 2); ImGui::SameLine();
    ImGui::RadioButton("Напруження", &activeTab, 3);
    ImGui::Separator();

    switch (activeTab) {
        case 0: renderInputPanel(state); break;
        case 1: renderMeshPanel(state); break;
        case 2: renderResultsPanel(state); break;
        case 3: renderStressPanel(state); break;
    }

    // Статус
    ImGui::Separator();
    if (state.currentStep > 0 && state.progress > 0) {
        ImGui::ProgressBar(state.progress / 100.0, ImVec2(-1, 20), 
                          (std::to_string((int)state.progress) + "%").c_str());
    }
    ImGui::Text("Статус: %s", state.statusMsg.c_str());

    ImGui::End();

    // Права панель - управління камерою
    renderControlsPanel(state);

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

void GUI::renderInputPanel(AppState& state) {
    ImGui::Text("Геометрія паралелепіпеда");
    ImGui::Separator();

    ImGui::InputDouble("Довжина a (X)", &state.a, 0.1, 1.0, "%.1f");
    ImGui::InputDouble("Ширина b (Y)", &state.b, 0.1, 1.0, "%.1f");
    ImGui::InputDouble("Висота c (Z)", &state.c, 0.1, 1.0, "%.1f");

    ImGui::Separator();
    ImGui::Text("Параметри сітки");
    ImGui::InputInt("Розбиття nx", &state.nx, 1, 2);
    ImGui::InputInt("Розбиття ny", &state.ny, 1, 2);
    ImGui::InputInt("Розбиття nz", &state.nz, 1, 1);
    
    if (state.nx < 1) state.nx = 1;
    if (state.ny < 1) state.ny = 1;
    if (state.nz < 1) state.nz = 1;

    ImGui::Separator();
    ImGui::Text("Характеристики матеріалу");
    ImGui::InputDouble("Модуль Юнга E", &state.E, 1e9, 1e11, "%.3e");
    ImGui::InputDouble("Коеф. Пуассона nu", &state.nu, 0.01, 0.1, "%.2f");

    ImGui::Separator();
    ImGui::Text("Навантаження");
    ImGui::InputDouble("Тиск P (верхня грань)", &state.P, 10.0, 100.0, "%.1f");

    ImGui::Separator();
    if (ImGui::Button("Згенерувати сітку", ImVec2(-1, 30))) {
        state.statusMsg = "Генерація сітки...";
        state.progress = 10.0;
        state.currentStep = 1;

        // Оновлюємо Mesh
        if (state.mesh) delete state.mesh;
        state.mesh = new Mesh(state.a, state.b, state.c, state.nx, state.ny, state.nz);
        state.mesh->generateParallelepiped();
        
        // Встановлюємо навантаження
        for (auto& face : state.mesh->loadedFaces) {
            face.pressure = state.P;
        }

        state.meshGenerated = true;
        state.femComputed = false;
        state.statusMsg = "Сітка згенерована. Натисніть 'Розрахувати'";
        state.progress = 20.0;
    }

    if (state.meshGenerated) {
        ImGui::Text("Вузлів: %d, Елементів: %d", state.mesh->nqp, state.mesh->nel);
        ImGui::Text("Півширина стрічки: %d", state.mesh->ng);
    }

    ImGui::Separator();
    if (ImGui::Button("Повний розрахунок МСЕ", ImVec2(-1, 40)) && state.meshGenerated) {
        state.statusMsg = "Обчислення базисних функцій...";
        state.progress = 30.0;
        state.currentStep = 2;

        // Крок 2: Базисні функції
        if (state.DFIABG) delete state.DFIABG;
        if (state.DPSITE) delete state.DPSITE;
        state.DFIABG = new std::array<std::array<std::array<double, 20>, 3>, 27>();
        state.DPSITE = new std::array<std::array<std::array<double, 8>, 2>, 9>();
        Gauss::computeDFIABG(*state.DFIABG);
        Gauss::computeDPSITE(*state.DPSITE);

        state.progress = 50.0;
        state.statusMsg = "Збірка глобальної системи...";
        state.currentStep = 3;

        // Крок 3-4: Збірка та крайові умови
        if (state.assembly) delete state.assembly;
        state.assembly = new Assembly();
        state.assembly->initialize(state.mesh->nqp, state.mesh->ng);
        state.assembly->assemble(*state.mesh, *state.DFIABG, *state.DPSITE);
        state.assembly->applyBoundaryConditions(state.mesh->fixedNodes);

        state.progress = 70.0;
        state.statusMsg = "Розв'язування СЛАР...";
        state.currentStep = 5;

        // Крок 5: Розв'язування
        state.assembly->solve();

        state.progress = 85.0;
        state.statusMsg = "Обчислення напружень...";
        state.currentStep = 6;

        // Крок 6: Напруження
        if (state.sigma) delete state.sigma;
        if (state.vonMises) delete state.vonMises;
        if (state.sigma1) delete state.sigma1;
        if (state.sigma2) delete state.sigma2;
        if (state.sigma3) delete state.sigma3;
        
        state.sigma = new std::vector<std::array<double, 6>>();
        state.vonMises = new std::vector<double>();
        state.sigma1 = new std::vector<double>();
        state.sigma2 = new std::vector<double>();
        state.sigma3 = new std::vector<double>();

        Stress::computeStresses(*state.mesh, state.assembly->U, *state.DFIABG, *state.sigma);
        Stress::computeVonMises(*state.sigma, *state.vonMises);
        Stress::computePrincipalStresses(*state.sigma, *state.sigma1, *state.sigma2, *state.sigma3);

        state.femComputed = true;
        state.progress = 100.0;
        state.currentStep = 7;
        state.statusMsg = "Розрахунок завершено!";
    }
}

void GUI::renderMeshPanel(AppState& state) {
    if (!state.mesh || !state.meshGenerated) {
        ImGui::Text("Спочатку згенеруйте сітку на вкладці 'Введення'");
        return;
    }

    ImGui::Text("Інформація про сітку");
    ImGui::Separator();

    ImGui::Text("Розміри: %.1f x %.1f x %.1f", state.mesh->a, state.mesh->b, state.mesh->c);
    ImGui::Text("Вузлів (nqp): %d", state.mesh->nqp);
    ImGui::Text("Елементів (nel): %d", state.mesh->nel);
    ImGui::Text("Півширина стрічки (ng): %d", state.mesh->ng);
    ImGui::Text("Закріплених вузлів: %zu", state.mesh->fixedNodes.size());
    ImGui::Text("Навантажених граней: %zu", state.mesh->loadedFaces.size());

    ImGui::Separator();
    ImGui::Text("Вузли сітки");

    // Фільтр
    ImGui::InputText("Фільтр", filterBuffer, sizeof(filterBuffer));
    
    // Таблиця вузлів
    static ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | 
                                    ImGuiTableFlags_ScrollY | ImGuiTableFlags_Sortable;
    
    if (ImGui::BeginTable("nodes_table", 4, flags, ImVec2(0, 300))) {
        ImGui::TableSetupColumn("№", ImGuiTableColumnFlags_DefaultSort);
        ImGui::TableSetupColumn("X");
        ImGui::TableSetupColumn("Y");
        ImGui::TableSetupColumn("Z");
        ImGui::TableHeadersRow();

        int maxShow = std::min(100, state.mesh->nqp);
        for (int i = 0; i < maxShow; ++i) {
            auto& node = state.mesh->nodes[i];
            
            // Фільтр
            std::string rowStr = std::to_string(i) + " " + 
                std::to_string(node.x()) + " " + 
                std::to_string(node.y()) + " " + 
                std::to_string(node.z());
            if (strlen(filterBuffer) > 0) {
                std::string filter(filterBuffer);
                std::transform(filter.begin(), filter.end(), filter.begin(), ::tolower);
                std::string lowerRow = rowStr;
                std::transform(lowerRow.begin(), lowerRow.end(), lowerRow.begin(), ::tolower);
                if (lowerRow.find(filter) == std::string::npos) continue;
            }

            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::Text("%d", i);
            ImGui::TableNextColumn(); ImGui::Text("%.3f", node.x());
            ImGui::TableNextColumn(); ImGui::Text("%.3f", node.y());
            ImGui::TableNextColumn(); ImGui::Text("%.3f", node.z());
        }
        ImGui::EndTable();
    }
    ImGui::Text("Показано перші %d вузлів", std::min(100, state.mesh->nqp));
    if (state.mesh->nqp > 100) {
        ImGui::Text("... та ще %d вузлів", state.mesh->nqp - 100);
    }

    ImGui::Separator();
    if (ImGui::CollapsingHeader("Масив зв'язності NT")) {
        if (ImGui::BeginTable("nt_table", 21, ImGuiTableFlags_Borders | ImGuiTableFlags_ScrollY, ImVec2(0, 200))) {
            ImGui::TableSetupColumn("Елемент");
            for (int j = 0; j < 20; ++j) {
                ImGui::TableSetupColumn(std::to_string(j+1).c_str());
            }
            ImGui::TableHeadersRow();

            int maxElem = std::min(10, state.mesh->nel);
            for (int e = 0; e < maxElem; ++e) {
                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::Text("%d", e);
                for (int j = 0; j < 20; ++j) {
                    ImGui::TableNextColumn(); ImGui::Text("%d", state.mesh->elements[e][j]);
                }
            }
            ImGui::EndTable();
        }
    }
}

void GUI::renderResultsPanel(AppState& state) {
    if (!state.femComputed || !state.assembly) {
        ImGui::Text("Спочатку виконайте розрахунок на вкладці 'Введення'");
        return;
    }

    ImGui::Text("Результати розрахунку");
    ImGui::Separator();

    ImGui::RadioButton("Переміщення", &activeResultTab, 0); ImGui::SameLine();
    ImGui::RadioButton("Деформація", &activeResultTab, 1);

    if (activeResultTab == 0) {
        // Таблиця переміщень
        ImGui::InputText("Фільтр", filterBuffer, sizeof(filterBuffer));
        
        if (ImGui::BeginTable("disp_table", 5, 
            ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY | ImGuiTableFlags_Sortable,
            ImVec2(0, 400))) {
            
            ImGui::TableSetupColumn("Вузол");
            ImGui::TableSetupColumn("ux");
            ImGui::TableSetupColumn("uy");
            ImGui::TableSetupColumn("uz");
            ImGui::TableSetupColumn("|u|");
            ImGui::TableHeadersRow();

            int maxShow = std::min(200, state.mesh->nqp);
            for (int i = 0; i < maxShow; ++i) {
                auto u = state.assembly->getDisplacement(i);
                double mag = u.norm();

                // Фільтр
                if (strlen(filterBuffer) > 0) {
                    std::string filter(filterBuffer);
                    if (std::to_string(i).find(filter) == std::string::npos) continue;
                }

                // Пропускаємо нульові (нев'язані вузли)
                if (mag < 1e-20) continue;

                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::Text("%d", i);
                ImGui::TableNextColumn(); ImGui::Text("%.4e", u.x());
                ImGui::TableNextColumn(); ImGui::Text("%.4e", u.y());
                ImGui::TableNextColumn(); ImGui::Text("%.4e", u.z());
                ImGui::TableNextColumn(); ImGui::Text("%.4e", mag);
            }
            ImGui::EndTable();
        }
    } else {
        // Візуалізація деформованої форми (відображається в OpenGL)
        ImGui::Text("Деформована форма:");
        ImGui::Text("Натисніть 'D' на 3D сцені для показу деформації");
        ImGui::Separator();
        
        if (state.mesh->nodes.size() > 0) {
            // Знаходимо максимальне переміщення
            double maxDisp = 0;
            for (int i = 0; i < state.mesh->nqp; ++i) {
                double mag = state.assembly->getDisplacement(i).norm();
                if (mag > maxDisp) maxDisp = mag;
            }
            ImGui::Text("Макс. переміщення: %.4e", maxDisp);
            ImGui::Text("Розмір моделі: ~%.1f", std::max({state.mesh->a, state.mesh->b, state.mesh->c}));
            double relDisp = maxDisp / std::max({state.mesh->a, state.mesh->b, state.mesh->c});
            ImGui::Text("Відносна деформація: %.4e", relDisp);
        }
    }
}

void GUI::renderStressPanel(AppState& state) {
    if (!state.femComputed || !state.sigma) {
        ImGui::Text("Спочатку виконайте розрахунок на вкладці 'Введення'");
        return;
    }

    ImGui::Text("Напруження у вузлах");
    ImGui::Separator();

    // Статистика
    double max_vM = *std::max_element(state.vonMises->begin(), state.vonMises->end());
    double min_vM = *std::min_element(state.vonMises->begin(), state.vonMises->end());
    double max_s1 = *std::max_element(state.sigma1->begin(), state.sigma1->end());
    double min_s3 = *std::min_element(state.sigma3->begin(), state.sigma3->end());

    ImGui::Text("Макс. напр. за Мізесом: %.4f", max_vM);
    ImGui::Text("Макс. головне σ₁: %.4f", max_s1);
    ImGui::Text("Мін. головне σ₃: %.4f", min_s3);
    ImGui::Separator();

    ImGui::RadioButton("Компоненти", &activeResultTab, 0); ImGui::SameLine();
    ImGui::RadioButton("Головні", &activeResultTab, 1); ImGui::SameLine();
    ImGui::RadioButton("Мізес", &activeResultTab, 2);
    ImGui::Separator();

    ImGui::InputText("Фільтр", filterBuffer, sizeof(filterBuffer));

    if (ImGui::BeginTable("stress_table", 8,
        ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY,
        ImVec2(0, 400))) {

        if (activeResultTab == 0) {
            ImGui::TableSetupColumn("Вузол");
            ImGui::TableSetupColumn("σ_xx");
            ImGui::TableSetupColumn("σ_yy");
            ImGui::TableSetupColumn("σ_zz");
            ImGui::TableSetupColumn("σ_xy");
            ImGui::TableSetupColumn("σ_yz");
            ImGui::TableSetupColumn("σ_xz");
            ImGui::TableSetupColumn("σ_vM");
        } else if (activeResultTab == 1) {
            ImGui::TableSetupColumn("Вузол");
            ImGui::TableSetupColumn("σ₁");
            ImGui::TableSetupColumn("σ₂");
            ImGui::TableSetupColumn("σ₃");
        } else {
            ImGui::TableSetupColumn("Вузол");
            ImGui::TableSetupColumn("σ_vM");
        }
        ImGui::TableHeadersRow();

        int maxShow = std::min(200, state.mesh->nqp);
        for (int i = 0; i < maxShow; ++i) {
            // Фільтр
            if (strlen(filterBuffer) > 0) {
                std::string filter(filterBuffer);
                if (std::to_string(i).find(filter) == std::string::npos) continue;
            }

            // Пропускаємо нульові
            if ((*state.vonMises)[i] < 1e-10) continue;

            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::Text("%d", i);

            if (activeResultTab == 0) {
                ImGui::TableNextColumn(); ImGui::Text("%.4f", (*state.sigma)[i][0]);
                ImGui::TableNextColumn(); ImGui::Text("%.4f", (*state.sigma)[i][1]);
                ImGui::TableNextColumn(); ImGui::Text("%.4f", (*state.sigma)[i][2]);
                ImGui::TableNextColumn(); ImGui::Text("%.4f", (*state.sigma)[i][3]);
                ImGui::TableNextColumn(); ImGui::Text("%.4f", (*state.sigma)[i][4]);
                ImGui::TableNextColumn(); ImGui::Text("%.4f", (*state.sigma)[i][5]);
                ImGui::TableNextColumn(); ImGui::Text("%.4f", (*state.vonMises)[i]);
            } else if (activeResultTab == 1) {
                ImGui::TableNextColumn(); ImGui::Text("%.4f", (*state.sigma1)[i]);
                ImGui::TableNextColumn(); ImGui::Text("%.4f", (*state.sigma2)[i]);
                ImGui::TableNextColumn(); ImGui::Text("%.4f", (*state.sigma3)[i]);
            } else {
                ImGui::TableNextColumn(); ImGui::Text("%.4f", (*state.vonMises)[i]);
            }
        }
        ImGui::EndTable();
    }
}

void GUI::renderControlsPanel(AppState& state) {
    ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x - 250, 20));
    ImGui::SetNextWindowSize(ImVec2(250, 180));
    ImGui::Begin("Управління 3D", nullptr, 
                 ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize);

    ImGui::Text("Відображення:");
    ImGui::Checkbox("Вузли [N]", &state.meshGenerated); // dummy - handled by viewer
    ImGui::Checkbox("Елементи [E]", &state.meshGenerated);
    ImGui::Checkbox("Деформація [D]", &state.femComputed);
    
    ImGui::Separator();
    ImGui::Text("Масштаб деформації:");
    static float defScale = 0.01f;
    if (ImGui::SliderFloat("##scale", &defScale, 0.0001f, 10.0f, "%.4f")) {
        // Viewer отримає це через callback
    }
    
    ImGui::Separator();
    if (ImGui::Button("Скинути камеру [R]", ImVec2(-1, 25))) {
        // Viewer отримає через ключ R
    }

    ImGui::End();
}