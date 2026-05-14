#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include "engine/PIMDSystem.h"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

int main() {
    if (!glfwInit()) return -1;
    
    GLFWwindow* window = glfwCreateWindow(800, 800, "Thermodynamics: Dynamic Phase Transition (Melting)", NULL, NULL);
    if (!window) { glfwTerminate(); return -1; }
    glfwMakeContextCurrent(window);
    
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) return -1;

    // --- Start in a Deep Freeze ---
    // Start at T = 0.05f to form the perfect crystal lattice
    // Mass = 20.0f (Good inertia), Temp = 0.05f (Cold, but NOT zero!)
    // 108 atoms for a perfect 3x3x3 FCC grid
    PIMDSystem system(108, 32, 15.0f, 0.05f);
    
    for (int i = 0; i < system.numAtoms; ++i) {
        system.atoms[i].mass = 1.0f; 
        system.atoms[i].charge = 0.0f; 
    }

    glEnable(GL_DEPTH_TEST);

    while (!glfwWindowShouldClose(window)) {
        // Dynamic Heating
        if (system.temperature < 100.0f) {
            system.temperature += 0.1f; 
            system.omegaP = std::sqrt((float)system.numBeads) * system.temperature;
        }

        system.update(0.001f);

        // --- NEW REAL-TIME UI ---
        std::string title = "PIMD Engine | Temp: " + std::to_string(system.temperature).substr(0, 5); 
        glfwSetWindowTitle(window, title.c_str());
        // ------------------------
        
        // ... (rest of your camera and rendering code) ... 

        glClearColor(0.05f, 0.05f, 0.05f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        float ratio = width / (float)height;
        glViewport(0, 0, width, height);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glm::mat4 proj = glm::perspective(glm::radians(45.0f), ratio, 0.1f, 100.0f);
        glLoadMatrixf(glm::value_ptr(proj));

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        // Keep the camera close enough to see the crystal, but far enough to see the gas
        float time = glfwGetTime();
        float camX = std::sin(time * 0.2f) * 18.0f; 
        float camZ = std::cos(time * 0.2f) * 18.0f;
        
        float center = system.L / 2.0f;
        glm::mat4 view = glm::lookAt(
            glm::vec3(camX + center, center + 3.0f, camZ + center), 
            glm::vec3(center, center, center),                      
            glm::vec3(0.0f, 1.0f, 0.0f)                             
        );
        glLoadMatrixf(glm::value_ptr(view));

        // Draw Boundary Box
        glColor3f(0.0f, 0.3f, 0.0f);
        glBegin(GL_LINES);
            float L = system.L;
            glVertex3f(0,0,0); glVertex3f(L,0,0); glVertex3f(L,0,0); glVertex3f(L,0,L);
            glVertex3f(L,0,L); glVertex3f(0,0,L); glVertex3f(0,0,L); glVertex3f(0,0,0);
            glVertex3f(0,L,0); glVertex3f(L,L,0); glVertex3f(L,L,0); glVertex3f(L,L,L);
            glVertex3f(L,L,L); glVertex3f(0,L,L); glVertex3f(0,L,L); glVertex3f(0,L,0);
            glVertex3f(0,0,0); glVertex3f(0,L,0); glVertex3f(L,0,0); glVertex3f(L,L,0);
            glVertex3f(L,0,L); glVertex3f(L,L,L); glVertex3f(0,0,L); glVertex3f(0,L,L);
        glEnd();

        // Draw Atoms
        for (int i = 0; i < system.numAtoms; ++i) {
            float r = (float)(i % 5) / 5.0f + 0.2f;
            float g = (float)(i % 7) / 7.0f + 0.2f;
            float b = (float)(i % 3) / 3.0f + 0.4f;
            glColor3f(r, g, b);
            
            glBegin(GL_LINES); 
            for (int k = 0; k < system.numBeads; ++k) {
                glm::vec3 p1 = system.atoms[i].beads[k].pos;
                glm::vec3 p2 = system.atoms[i].beads[(k + 1) % system.numBeads].pos;
                glVertex3f(p1.x, p1.y, p1.z);
                glVertex3f(p2.x, p2.y, p2.z);
            }
            glEnd();

            glPointSize(4.0f); 
            glBegin(GL_POINTS);
            for (auto& bead : system.atoms[i].beads) glVertex3f(bead.pos.x, bead.pos.y, bead.pos.z);
            glEnd();
        }

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}