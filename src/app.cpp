#include <iostream>
#include <GLFW/glfw3.h>
#include "mujoco/mujoco.h"
#include <dtMath/dtMath.h>

// mouse interaction
char err[1000];
mjModel *m;
mjData *d;
mjVFS vfs;
mjvCamera cam;
mjvOption opt;
mjvScene scn;
mjrContext ctx;

// mouse interaction
bool button_left = false;
bool button_middle = false;
bool button_right = false;
double lastx = 0;
double lasty = 0;

// keyboard callback
void keyboard(GLFWwindow *window, int key, int scancode, int act, int mods)
{
    // backspace: reset simulation
    if (act == GLFW_PRESS && key == GLFW_KEY_BACKSPACE)
    {
        mj_resetData(m, d);
        mj_forward(m, d);
    }
}

// mouse button callback
void mouse_button(GLFWwindow *window, int button, int act, int mods)
{
    // update button state
    button_left = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS);
    button_middle = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS);
    button_right = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS);

    // update mouse position
    glfwGetCursorPos(window, &lastx, &lasty);
}

// mouse move callback
void mouse_move(GLFWwindow *window, double xpos, double ypos)
{
    // no buttons down: nothing to do
    if (!button_left && !button_middle && !button_right)
    {
        return;
    }

    // compute mouse displacement, save
    double dx = xpos - lastx;
    double dy = ypos - lasty;
    lastx = xpos;
    lasty = ypos;

    // get current window size
    int width, height;
    glfwGetWindowSize(window, &width, &height);

    // get shift key state
    bool mod_shift = (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
                      glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS);

    // determine action based on mouse button
    mjtMouse action;
    if (button_right)
    {
        action = mod_shift ? mjMOUSE_MOVE_H : mjMOUSE_MOVE_V;
    }
    else if (button_left)
    {
        action = mod_shift ? mjMOUSE_ROTATE_H : mjMOUSE_ROTATE_V;
    }
    else
    {
        action = mjMOUSE_ZOOM;
    }

    // move camera
    mjv_moveCamera(m, action, dx / height, dy / height, &scn, &cam);
}

// scroll callback
void scroll(GLFWwindow *window, double xoffset, double yoffset)
{
    // emulate vertical mouse motion = 5% of window height
    mjv_moveCamera(m, mjMOUSE_ZOOM, 0, -0.05 * yoffset, &scn, &cam);
}

int main(int argc, char **argv)
{

    // load a model
    // m = mj_loadXML("model/humanoid/humanoid.xml", NULL, err, 1000);
    m = mj_loadXML("model/quadip/QuadIP.xml", NULL, err, 1000);
    if (!m)
    {
        std::cout << err << std::endl;
        return -1;
    }

    // make data
    d = mj_makeData(m);

    // init GLFW
    glfwInit();
    GLFWwindow *window = glfwCreateWindow(1200, 900, "Demo", NULL, NULL);
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    mjv_defaultCamera(&cam);
    mjv_defaultOption(&opt);
    // mjv_defaultPerturb(&pert);
    mjv_defaultScene(&scn);
    mjr_defaultContext(&ctx);

    // create scene & context
    mjv_makeScene(m, &scn, 2000);
    mjr_makeContext(m, &ctx, mjFONTSCALE_150);

    // install GLFW mouse and keyboard callbacks
    glfwSetKeyCallback(window, keyboard);
    glfwSetCursorPosCallback(window, mouse_move);
    glfwSetMouseButtonCallback(window, mouse_button);
    glfwSetScrollCallback(window, scroll);

    // run simulation
    while (!glfwWindowShouldClose(window))
    {
        mjtNum simstart = d->time;

        while (d->time - simstart < 1.0 / 60.0)
        {
            std::cout << d->time << std::endl;
            double tau[12] = {0.0, -0.1, 0.0, 0.0, -0.1, 0.0, 0.0, -0.1, 0.0, 0.0, -0.1, 0.0};
            std::copy(&tau[0], &tau[0] + 12, d->ctrl);
            mj_step(m, d);
        }

        // get framebuffer viewport
        mjrRect viewport = {0, 0, 0, 0};
        glfwGetFramebufferSize(window, &viewport.width, &viewport.height);

        // update scene and render
        mjv_updateScene(m, d, &opt, NULL, &cam, mjCAT_ALL, &scn);
        mjr_render(viewport, &scn, &ctx);

        // swap OpenGL buffers
        glfwSwapBuffers(window);

        // process pending GUI events
        glfwPollEvents();
    }

    // free scene and context
    mjv_freeScene(&scn);
    mjr_freeContext(&ctx);

    // free model and data
    mj_deleteData(d);
    mj_deleteModel(m);

// terminate GLFW
#if defined(__APPLE__) || defined(_WIN32)
    glfwTerminate();
#endif

    std::cout << "Good-bye~" << std::endl;
    return 0;
}