

#include <iostream>
using namespace std;

#include "objs.hpp"



int main() {
    chargedParticle cp{{-2,0,0}, {2,0,0}, 1.602e-19, 1.67e-27};
    magDipole m{{0,0,0},{10,0,0.1}};
    double time{0.8};
    double timestep{1e-5};
    int writeskip{100};

    for(int i = 0; i<int(time/timestep); i++) {
        cp.timestep(m, timestep);
        if(i % writeskip == 0) cout << cp.loc << "\n";
    }

}