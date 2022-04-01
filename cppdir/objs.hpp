#pragma once
#include <array>
#include <initializer_list>
#include <iostream>
#include <vector>

using namespace std;

#define mu0 4*3.141592*1e-7


struct vec {
    std::array<double, 3> loc;
    
    vec(double x, double y, double z);
    vec(const double(elems)[3]);
    vec(const std::array<double,3> elems);
    vec(const vec& rhs);

    vec& operator*=(const double &rhs);
    vec& operator/=(const double &rhs);
    vec operator*(const double &rhs) const;
    vec operator/(const double &rhs) const;
    vec& operator+=(const vec &rhs);
    vec& operator-=(const vec &rhs);
    vec operator+(const vec &rhs) const;
    vec operator-(const vec &rhs) const;

    double len() const;
    vec unit() const;

    double dot(const vec& rhs) const;
};

vec cross(const vec& a, const vec& b);

ostream& operator<<(ostream& cout, const vec& rhs);
pair<vec, vec> operator*(pair<vec, vec> elem, double rhs);
pair<vec, vec> operator*(double lhs, pair<vec, vec> elem);
istream& operator>>(istream& cin, vec& rhs);

class magDipole {
    public:
    vec loc;
    vec m;

    magDipole(vec loc, vec m);

    vec B(const vec &loc);

    private:

};

struct chargedParticle {

    public:

    vector<magDipole*> reacting_dipoles;
    vec loc;
    vec vel;
    double charge;
    double mass;
    chargedParticle(const chargedParticle& rhs);
    chargedParticle(const vec& loc, const vec& vel, double charge, double mass);
    
    chargedParticle& operator+=(pair<const vec&, const vec&> diff);
    chargedParticle& operator-=(pair<const vec&, const vec&> diff);
    chargedParticle operator+(pair<const vec&, const vec&> diff);
    chargedParticle operator-(pair<const vec&, const vec&> diff);

    pair<vec,vec> state();

    vec F(magDipole dipole);
    pair<vec,vec> timediff();
    void timestep(double h=1e-3);
};