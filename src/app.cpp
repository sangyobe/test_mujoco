#include <iostream>
#include <app.h>
#include "mujoco/mujoco.h"

int main(int argc, char **argv)
{
    char err[1000];
    mjModel *m;
    mjData *d;
    mjVFS vfs;

    std::cout << "Hello!" << std::endl;

    // load a model
    m = mj_loadXML("models/humanoid.xml", NULL, err, 1000);
    if (!m)
    {
        std::cout << err << std::endl;
        return -1;
    }

    // make data
    d = mj_makeData(m);

    // run simulation
    while (d->time < 10)
    {
        std::cout << d->time << std::endl;
        mj_step(m, d);
    }

    // free model and data
    mj_deleteData(d);
    mj_deleteModel(m);

    std::cout << "Good-bye~" << std::endl;
    return 0;
}