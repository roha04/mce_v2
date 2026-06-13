cmake_minimum_required(VERSION 3.14)
project(MCE_Project VERSION 1.0.0 LANGUAGES CXX C)

# Встановлюємо стандарт C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Включаємо модуль для автоматичного завантаження залежностей
include(FetchContent)

message(STATUS "Завантаження Eigen (математика)...")
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG        3.4.0
)
FetchContent_MakeAvailable(eigen)

message(STATUS "Завантаження GLFW (створення вікна)...")
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
  glfw
  GIT_REPOSITORY https://github.com/glfw/glfw.git
  GIT_TAG        3.3.8
)
FetchContent_MakeAvailable(glfw)

# Завантаження ImGui
message(STATUS "Завантаження ImGui (GUI)...")
FetchContent_Declare(
  imgui
  GIT_REPOSITORY https://github.com/ocornut/imgui.git
  GIT_TAG        v1.91.0
)
FetchContent_MakeAvailable(imgui)

# Створюємо виконуваний файл
add_executable(MCE_App 
    src/main.cpp
    src/FEM/Mesh.cpp
    src/FEM/Gauss.cpp
    src/FEM/Element.cpp
    src/FEM/Assembly.cpp
    src/FEM/Stress.cpp
    src/Viewer.cpp
    src/GUI.cpp
)

# Вказуємо, де шукати заголовкові файли (.h)
target_include_directories(MCE_App PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
    ${imgui_SOURCE_DIR}
    ${imgui_SOURCE_DIR}/backends
)

# Підключаємо бібліотеки до нашої програми
target_link_libraries(MCE_App PRIVATE
    glfw
    Eigen3::Eigen
)

# Додаємо ImGui файли
target_sources(MCE_App PRIVATE
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
)
