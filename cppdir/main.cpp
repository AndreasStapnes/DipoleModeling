

#include <iostream>
#include <fstream>
using namespace std;

#include "objs.hpp"



int main() {
    ifstream file{};
    file.open("in.txt");
    double mass{}, charge{}, time{}, timestep{};
    int writeskip{};
    vec loc{{0,0,0}}, vel{{0,0,0}}, moment{{0,0,0}};
    string dump;
    file >> mass >> dump;
    file >> charge >> dump;
    file >> loc >> dump;
    file >> vel >> dump; 
    chargedParticle cp{loc, vel, charge, mass};
    file >> loc >> dump;
    file >> moment >> dump;

    magDipole m{loc,moment};
    cp.reacting_dipoles.push_back(&m);
    file >> time >> dump;
    file >> timestep >> dump;
    file >> writeskip >> dump;

    for(int i = 0; i<int(time/timestep); i++) {
        cp.timestep(timestep);
        if(i % writeskip == 0) cout << cp.loc << " " << cp.vel << "\n";
    }

}