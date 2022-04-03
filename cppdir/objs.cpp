#include "objs.hpp"
#include <initializer_list>
#include <numeric>
#include <math.h>
#include <iostream>
#include <exception>

vec::vec(double x, double y, double z): loc{x,y,z} {}
vec::vec(const double (elems)[3]) {
    for(int i = 0; i<3; i++) {
        this->loc[i] = *elems;
        elems ++;
    }
}
vec::vec(const std::array<double, 3> elems): loc{elems} {}
vec::vec(const vec &rhs):loc{rhs.loc} {}



vec& vec::operator+=(const vec& rhs) {
    for(int i = 0; i<3; i++) this->loc[i] += rhs.loc[i];
    return *this;
}
vec& vec::operator-=(const vec& rhs) {
    for(int i = 0; i<3; i++) this->loc[i] -= rhs.loc[i];
    return *this;
}
vec vec::operator+(const vec &rhs) const{
    vec result{*this};
    result += rhs;
    return result;
}
vec vec::operator-(const vec &rhs) const{
    vec result{*this};
    result -= rhs;
    return result;
}
vec& vec::operator*=(const double &rhs){
    for(int i = 0; i<3; i++) this->loc[i]*=rhs;
    return *this;
}
vec& vec::operator/=(const double &rhs){
    (*this)*=(1/rhs);
    return *this;
}
vec vec::operator*(const double &rhs) const {
    vec result{*this};
    result *= rhs;
    return result;
}
vec vec::operator/(const double &rhs) const {
    vec result{*this};
    result /= rhs;
    return result;
}

double vec::len() const {
    double sqsum{};
    for(int i = 0; i<3; i++) sqsum += (this->loc[i])*(this->loc[i]);
    return sqrt(sqsum);
}
vec vec::unit() const {
    return (*this) / this->len();
}

double vec::dot(const vec& rhs) const{
    double sum{0};
    for(int i = 0; i<3; i++) {
        sum += this->loc[i]*rhs.loc[i];
    }
    return sum;
}


vec cross(const vec& a, const vec& b) {
    double xa = a.loc[0]; double ya = a.loc[1]; double za = a.loc[2];
    double xb = b.loc[0]; double yb = b.loc[1]; double zb = b.loc[2];
    return vec{(ya*zb-za*yb), (za*xb-xa*zb), (xa*yb-ya*xb)};
}


ostream& operator<<(ostream& ostr, const vec& rhs) {
    ostr << "(";
    for(int i = 0; i<2; i++) ostr << rhs.loc[i] << ",";
    ostr << rhs.loc[2] << ")";
    return ostr;
}

pair<vec, vec> operator*(pair<vec, vec> elem, double rhs) {
    return pair<vec, vec>{elem.first*rhs, elem.second*rhs};
}
pair<vec, vec> operator*(double lhs, pair<vec, vec> elem) {
    return elem*lhs;
}

istream& operator>>(istream& cin, vec& rhs) {
    char l_p, r_p, c1, c2;
    double x,y,z;
    cin >> l_p >> x >> c1 >> y >> c2 >> z >> r_p;
    if(!cin || l_p != '(' || c1 != ',' || c2 != ',' || r_p != ')') {
     throw invalid_argument{"text did not match vector-format (..,..,..)"};}
    rhs.loc[0] = x;
    rhs.loc[1] = y;
    rhs.loc[2] = z;
    return cin;
}


magDipole::magDipole(vec loc, vec m):m{m}, loc{loc} {}

vec magDipole::B(const vec &loc) {
    //mu0/(4*np.pi*r_siz*r_siz*r_siz) * (3*r.dir()*np.dot(r.dir(), self)
    vec r {loc - this->loc};
    //cout << "r=" << r << "\n";
    double r_siz{r.len()};
    vec rhat{r.unit()};
    //cout << "rhat=" << rhat << "\n";
    vec returnval = (rhat*3*(rhat.dot(this->m)) - this->m) * mu0/(4*M_PI*r_siz*r_siz*r_siz);
    //cout << "returnval=" << returnval << "\n";
    return returnval;
}

chargedParticle::chargedParticle(const chargedParticle& rhs) : loc{rhs.loc}, vel{rhs.vel}, charge{rhs.charge}, mass{rhs.mass}, reacting_dipoles{rhs.reacting_dipoles} {}
chargedParticle::chargedParticle(const vec& loc, const vec& vel, double charge, double mass): loc{loc}, vel{vel}, charge{charge}, mass{mass}, reacting_dipoles{} {}
chargedParticle& chargedParticle::operator+=(pair<const vec&, const vec&> diff) {
    this->loc += diff.first;
    this->vel += diff.second;
    //cout << "!vel="<<(this->vel)<<" !";
    return *this;
}
chargedParticle& chargedParticle::operator-=(pair<const vec&, const vec&> diff) {
    this->loc -= diff.first;
    this->vel -= diff.second;
    return *this;
}
chargedParticle chargedParticle::operator+(pair<const vec&, const vec&> diff) {
    chargedParticle result{*this};
    result += diff;
    return result;
}
chargedParticle chargedParticle::operator-(pair<const vec&, const vec&> diff) {
    chargedParticle result{*this};
    result -= diff;
    return result;
}

pair<vec,vec> chargedParticle::state() {
    return pair<vec,vec>{this->loc, this->vel};
}

vec chargedParticle::F(magDipole dipole) {
    return cross(this->vel, dipole.B(this->loc))*this->charge;
}

pair<vec,vec> chargedParticle::timediff() {
    vec F{0,0,0};
    for(magDipole* dp : this->reacting_dipoles) {F += this->F(*dp);}
    return pair<vec,vec>{this->vel, F/(this->mass)};
}

void chargedParticle::timestep(double h) {
    auto f1 = (*this).timediff();
    auto f2 = ((*this) + f1*(h/2)).timediff();
    
    auto f3 = ((*this) + f2*(h/2)).timediff();
    auto f4 = ((*this) + f3*h).timediff();
    //cout << "f2="<<f2.first << "," << f2.second <<" ";
    this->loc += ((f1.first + f2.first*2 + f3.first*2 + f4.first)*h/6);
    this->vel += ((f1.second + f2.second*2 + f3.second*2 + f4.second)*h/6);
}