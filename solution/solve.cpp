#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;
namespace fs = std::filesystem;
using C = complex<double>;
const double PI0 = acos(-1.0);

struct Station { string id; double angle; C gain; double pole_tau, zero_tau, cal_rms; };
struct Obs { int id, si; double f; C raw, z; double la, ph; };
struct Query { int id, si; double f; };
struct ModelCfg {
    int n = 3, disp = 1;
    vector<double> init, lo, hi;
};
struct P { double rho[3], h[2], m, tau, c; };
struct CalSample { double f; C ref, raw; };
struct ResponseFit { C gain; double pole_tau, zero_tau, rms; };

static double clamp01(double x) { return max(0.0, min(1.0, x)); }
static double wrapdeg(double x) {
    x = fmod(x + 180.0, 360.0);
    if (x < 0) x += 360.0;
    return x - 180.0;
}
static double modphase(double a, double b) { return wrapdeg(a - b); }

static vector<string> toks(const string& line) {
    istringstream in(line);
    vector<string> out; string x;
    while (in >> x) out.push_back(x);
    return out;
}

static P decode(const vector<double>& y, const ModelCfg& c) {
    P p{};
    auto logmap = [&](int k) { return exp(log(c.lo[k]) + y[k] * (log(c.hi[k]) - log(c.lo[k]))); };
    p.rho[0] = logmap(0); p.rho[1] = logmap(1); p.rho[2] = logmap(2);
    p.h[0] = logmap(3); p.h[1] = logmap(4);
    p.m = c.lo[5] + y[5] * (c.hi[5] - c.lo[5]);
    p.tau = logmap(6);
    p.c = c.lo[7] + y[7] * (c.hi[7] - c.lo[7]);
    return p;
}

static C forward_z(double f, const P& p, int disp_layer, double mu0, double pi) {
    C iw(0.0, 2.0 * pi * f);
    C sig[3];
    for (int j = 0; j < 3; ++j) sig[j] = 1.0 / p.rho[j];
    C disp = pow(iw * p.tau, p.c);
    if (0 <= disp_layer && disp_layer < 3) sig[disp_layer] *= 1.0 + p.m / (1.0 + disp);
    C z = sqrt(iw * mu0 / sig[2]);
    for (int j = 1; j >= 0; --j) {
        C gamma = sqrt(iw * mu0 * sig[j]);
        C eta = sqrt(iw * mu0 / sig[j]);
        C t = tanh(gamma * p.h[j]);
        z = eta * (z + eta * t) / (eta + z * t);
    }
    return z;
}

static double huber(double n, double d) {
    return n <= d ? 0.5 * n * n : d * n - 0.5 * d * d;
}

struct Fit {
    double obj;
    vector<double> y;
    vector<double> residual;
};

static double objective(const vector<double>& y, const ModelCfg& cfg, const vector<Obs>& obs,
                        double mu0, double pi, double delta, vector<double>* rr = nullptr) {
    P p = decode(y, cfg);
    double out = 0.0;
    vector<double> r;
    if (rr) r.reserve(obs.size() * 2);
    for (const Obs& o : obs) {
        C z = forward_z(o.f, p, cfg.disp, mu0, pi);
        double a = log10(abs(z)) - o.la;
        double ph = modphase(atan2(z.imag(), z.real()) * 180.0 / pi, o.ph) * pi / 180.0;
        double n = hypot(a, ph);
        out += huber(n, delta);
        if (rr) { r.push_back(a); r.push_back(ph); }
    }
    if (rr) *rr = move(r);
    return out;
}

static bool solve_linear(vector<vector<double>> a, vector<double> b, vector<double>& x) {
    int n = static_cast<int>(b.size());
    for (int c = 0; c < n; ++c) {
        int piv = c;
        for (int r = c + 1; r < n; ++r) if (fabs(a[r][c]) > fabs(a[piv][c])) piv = r;
        if (fabs(a[piv][c]) < 1e-14) return false;
        swap(a[piv], a[c]); swap(b[piv], b[c]);
        double q = a[c][c];
        for (int j = c; j < n; ++j) a[c][j] /= q;
        b[c] /= q;
        for (int r = 0; r < n; ++r) if (r != c) {
            double z = a[r][c];
            if (z == 0.0) continue;
            for (int j = c; j < n; ++j) a[r][j] -= z * a[c][j];
            b[r] -= z * b[c];
        }
    }
    x = b; return true;
}

static double clock_value(double x, const vector<double>& c) {
    double y = 0.0, p = 1.0;
    for (double a : c) { y += a * p; p *= x; }
    return y;
}

static double clock_slope(double x, const vector<double>& c) {
    double y = 0.0, p = 1.0;
    for (size_t k = 1; k < c.size(); ++k) { y += k * c[k] * p; p *= x; }
    return y;
}

static double frequency_from_tick(double tick, const vector<double>& c, double beacon_lo, double beacon_hi) {
    vector<double> roots;
    double lo = max(-20.0, beacon_lo - 1.5), hi = min(20.0, beacon_hi + 1.5);
    int scans = 2400;
    double prev_x = lo, prev_y = clock_value(prev_x, c) - tick;
    for (int i = 1; i <= scans; ++i) {
        double x = lo + (hi - lo) * i / scans;
        double y = clock_value(x, c) - tick;
        if (prev_y == 0.0 || y == 0.0 || (prev_y < 0.0) != (y < 0.0)) {
            double a = prev_x, b = x, fa = prev_y;
            for (int it = 0; it < 80; ++it) {
                double m = 0.5 * (a + b), fm = clock_value(m, c) - tick;
                if ((fa < 0.0) == (fm < 0.0)) { a = m; fa = fm; } else b = m;
            }
            roots.push_back(0.5 * (a + b));
        }
        prev_x = x; prev_y = y;
    }
    double chosen = roots.empty() ? 0.5 * (beacon_lo + beacon_hi) : roots.front();
    double best = numeric_limits<double>::infinity();
    double center = 0.5 * (beacon_lo + beacon_hi);
    for (double x : roots) {
        if (clock_slope(x, c) <= 0.0) continue;
        double penalty = (x < beacon_lo - 1e-9 || x > beacon_hi + 1e-9) ? 10.0 : 0.0;
        double score = penalty + fabs(x - center) * 1e-3 + fabs(clock_value(x, c) - tick);
        if (score < best) { best = score; chosen = x; }
    }
    return exp(chosen);
}

static ResponseFit fit_response(const vector<CalSample>& rows) {
    auto score = [&](double pole, double zero, C* gain_out) {
        C num(0.0, 0.0); double den = 0.0;
        for (const auto& r : rows) {
            C corrected = r.raw * C(1.0, 2.0 * PI0 * r.f * pole) / C(1.0, 2.0 * PI0 * r.f * zero);
            num += corrected * conj(r.ref); den += norm(r.ref);
        }
        C g = num / den; if (gain_out) *gain_out = g;
        double e = 0.0;
        for (const auto& r : rows) {
            C pred = g * C(1.0, 2.0 * PI0 * r.f * zero) / C(1.0, 2.0 * PI0 * r.f * pole) * r.ref;
            e += norm(pred - r.raw);
        }
        return e;
    };
    double best_pole = 0.0, best_zero = 0.0, best_score = numeric_limits<double>::infinity();
    for (int i = 0; i <= 60; ++i) {
        double pole = 0.02 * i / 60.0;
        for (int j = 0; j <= 48; ++j) {
            double zero = 0.012 * j / 48.0;
            double s = score(pole, zero, nullptr);
            if (s < best_score) { best_score = s; best_pole = pole; best_zero = zero; }
        }
    }
    double sp = 0.0005, sz = 0.00035;
    for (int it = 0; it < 90; ++it) {
        bool improved = false;
        for (int dim = 0; dim < 2; ++dim) {
            for (int dir : {-1, 1}) {
                double pole = best_pole + (dim == 0 ? dir * sp : 0.0);
                double zero = best_zero + (dim == 1 ? dir * sz : 0.0);
                if (pole < 0.0 || pole > 0.02 || zero < 0.0 || zero > 0.012) continue;
                double s = score(pole, zero, nullptr);
                if (s < best_score) { best_score = s; best_pole = pole; best_zero = zero; improved = true; }
            }
        }
        if (!improved) { sp *= 0.55; sz *= 0.55; }
    }
    C g; double e = score(best_pole, best_zero, &g);
    return {g, best_pole, best_zero, sqrt(e / max<size_t>(1, rows.size()))};
}

static Fit optimize(vector<double> start, const ModelCfg& cfg, const vector<Obs>& obs,
                    double mu0, double pi, double delta) {
    const int n = 8;
    for (double& x : start) x = clamp01(x);
    double cur = objective(start, cfg, obs, mu0, pi, delta);
    double lambda = 0.01;
    for (int iter = 0; iter < 180; ++iter) {
        vector<double> r;
        objective(start, cfg, obs, mu0, pi, delta, &r);
        int rows = static_cast<int>(r.size());
        vector<vector<double>> J(rows, vector<double>(n));
        for (int k = 0; k < n; ++k) {
            vector<double> yp = start;
            double step = 2e-5;
            yp[k] = clamp01(yp[k] + step);
            vector<double> rp;
            objective(yp, cfg, obs, mu0, pi, delta, &rp);
            double den = yp[k] - start[k];
            if (den == 0.0) den = step;
            for (int i = 0; i < rows; ++i) J[i][k] = (rp[i] - r[i]) / den;
        }
        vector<vector<double>> A(n, vector<double>(n, 0.0));
        vector<double> b(n, 0.0);
        for (int i = 0; i < rows; i += 2) {
            double norm = hypot(r[i], r[i + 1]);
            double w = norm <= delta ? 1.0 : delta / max(norm, 1e-15);
            for (int rr = 0; rr < 2; ++rr) {
                double wi = sqrt(w), ri = r[i + rr] * wi;
                for (int u = 0; u < n; ++u) {
                    double ju = J[i + rr][u] * wi;
                    b[u] -= ju * ri;
                    for (int v = 0; v < n; ++v) A[u][v] += ju * J[i + rr][v] * wi;
                }
            }
        }
        for (int k = 0; k < n; ++k) A[k][k] += lambda * (A[k][k] + 1e-8);
        vector<double> step;
        if (!solve_linear(A, b, step)) break;
        double sn = 0.0; for (double z : step) sn += z * z;
        vector<double> trial = start;
        for (int k = 0; k < n; ++k) trial[k] = clamp01(trial[k] + step[k]);
        double nxt = objective(trial, cfg, obs, mu0, pi, delta);
        if (nxt < cur) {
            start = move(trial); cur = nxt; lambda = max(1e-9, lambda * 0.45);
            if (sn < 1e-12) break;
        } else {
            lambda = min(1e8, lambda * 4.0);
        }
    }
    vector<double> rr;
    objective(start, cfg, obs, mu0, pi, delta, &rr);
    return {objective(start, cfg, obs, mu0, pi, delta), move(start), move(rr)};
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    fs::path manifest(argv[1]);
    fs::path dir = manifest.parent_path();
    map<string,string> files;
    {
        ifstream f(manifest);
        string line;
        while (getline(f, line)) {
            auto t = toks(line); if (t.size() >= 2 && t[0][0] != '#') files[t[0]] = t[1];
        }
    }
    int clock_order = 2;
    double mu0 = 4.0 * PI0 * 1e-7, pi = PI0, delta = 0.018;
    {
        ifstream f(dir / files["constants"]); string line;
        while (getline(f,line)) { auto t=toks(line); if(t.size()<2) continue; if(t[0]=="mu0") mu0=stod(t[1]); if(t[0]=="pi") pi=stod(t[1]); if(t[0]=="huber_delta") delta=stod(t[1]); if(t[0]=="clock_order") clock_order=stoi(t[1]); }
    }
    vector<double> clock(clock_order + 1, 0.0);
    double beacon_lo = 20.0, beacon_hi = -20.0;
    {
        int m = clock_order + 1;
        vector<vector<double>> normal(m, vector<double>(m, 0.0));
        vector<double> rhs(m, 0.0); string line; ifstream f(dir / files["clock"]);
        while (getline(f, line)) {
            auto t = toks(line); if (t.size() < 3 || t[0][0] == '#') continue;
            double tick = stod(t[1]), x = log(stod(t[2]));
            beacon_lo = min(beacon_lo, x); beacon_hi = max(beacon_hi, x);
            vector<double> b(m, 1.0);
            for (int i = 1; i < m; ++i) b[i] = b[i - 1] * x;
            for (int i = 0; i < m; ++i) { rhs[i] += b[i] * tick; for (int j = 0; j < m; ++j) normal[i][j] += b[i] * b[j]; }
        }
        solve_linear(normal, rhs, clock);
    }
    ModelCfg cfg;
    {
        ifstream f(dir / files["model"]); string line;
        while (getline(f,line)) {
            auto t=toks(line); if(t.empty() || t[0][0]=='#') continue;
            if(t[0]=="layer_count") cfg.n=stoi(t[1]);
            else if(t[0]=="dispersive_layer") cfg.disp=stoi(t[1]);
            else if(t[0]=="initial_rho") for(int i=0;i<3;++i) { if((int)cfg.init.size()<8) cfg.init.resize(8); cfg.init[i]=stod(t[i+1]); }
            else if(t[0]=="initial_thickness") { if(cfg.init.size()<8) cfg.init.resize(8); cfg.init[3]=stod(t[1]); cfg.init[4]=stod(t[2]); }
            else if(t[0]=="initial_m") { if(cfg.init.size()<8) cfg.init.resize(8); cfg.init[5]=stod(t[1]); }
            else if(t[0]=="initial_tau") { if(cfg.init.size()<8) cfg.init.resize(8); cfg.init[6]=stod(t[1]); }
            else if(t[0]=="initial_c") { if(cfg.init.size()<8) cfg.init.resize(8); cfg.init[7]=stod(t[1]); }
            else if(t[0]=="lower_rho") { cfg.lo.resize(8); for(int i=0;i<3;++i) cfg.lo[i]=stod(t[i+1]); }
            else if(t[0]=="upper_rho") { cfg.hi.resize(8); for(int i=0;i<3;++i) cfg.hi[i]=stod(t[i+1]); }
            else if(t[0]=="lower_thickness") { for(int i=0;i<2;++i) cfg.lo[3+i]=stod(t[i+1]); }
            else if(t[0]=="upper_thickness") { for(int i=0;i<2;++i) cfg.hi[3+i]=stod(t[i+1]); }
            else if(t[0]=="lower_m") cfg.lo[5]=stod(t[1]); else if(t[0]=="upper_m") cfg.hi[5]=stod(t[1]);
            else if(t[0]=="lower_tau") cfg.lo[6]=stod(t[1]); else if(t[0]=="upper_tau") cfg.hi[6]=stod(t[1]);
            else if(t[0]=="lower_c") cfg.lo[7]=stod(t[1]); else if(t[0]=="upper_c") cfg.hi[7]=stod(t[1]);
        }
    }
    vector<Station> stations; unordered_map<string,int> si;
    {
        ifstream f(dir / files["stations"]); string line;
        while(getline(f,line)){ auto t=toks(line); if(t.size()<2 || t[0][0]=='#') continue; si[t[0]]=stations.size(); stations.push_back({t[0],stod(t[1]),C(1,0),0,0,0}); }
    }
    {
        ifstream f(dir / files["calibration"]); string line; vector<vector<CalSample>> samples(stations.size());
        while(getline(f,line)){ auto t=toks(line); if(t.size()<6 || t[0][0]=='#') continue; int k=si[t[0]]; samples[k].push_back({stod(t[1]),C(stod(t[2]),stod(t[3])),C(stod(t[4]),stod(t[5]))}); }
        for(size_t k=0;k<stations.size();++k){ auto fit=fit_response(samples[k]); stations[k].gain=fit.gain; stations[k].pole_tau=fit.pole_tau; stations[k].zero_tau=fit.zero_tau; stations[k].cal_rms=fit.rms; }
    }
    vector<Obs> obs;
    {
        ifstream f(dir / files["observations"]); string line;
        while(getline(f,line)){ auto t=toks(line); if(t.size()<5 || t[0][0]=='#') continue; int id=stoi(t[0]), k=si[t[1]]; double freq=frequency_from_tick(stod(t[2]),clock,beacon_lo,beacon_hi); C raw(stod(t[3]),stod(t[4])); C z=raw*C(1.0,2.0*pi*freq*stations[k].pole_tau)/(stations[k].gain*C(1.0,2.0*pi*freq*stations[k].zero_tau))/exp(C(0,stations[k].angle*pi/180.0)); obs.push_back({id,k,freq,raw,z,log10(abs(z)),atan2(z.imag(),z.real())*180.0/pi}); }
    }
    vector<Query> queries;
    { ifstream f(dir / files["queries"]); string line; while(getline(f,line)){ auto t=toks(line); if(t.size()<3||t[0][0]=='#') continue; queries.push_back({stoi(t[0]),si[t[1]],frequency_from_tick(stod(t[2]),clock,beacon_lo,beacon_hi)}); } }
    vector<double> y0(8); for(int k=0;k<8;++k){ double a=cfg.init[k]; if(k==5||k==7) y0[k]=(a-cfg.lo[k])/(cfg.hi[k]-cfg.lo[k]); else y0[k]=(log(a)-log(cfg.lo[k]))/(log(cfg.hi[k])-log(cfg.lo[k])); y0[k]=clamp01(y0[k]); }
    Fit best{numeric_limits<double>::infinity(),{}, {}};
    for(int s=0;s<12;++s){ vector<double> st=y0; if(s==1) for(double& x:st)x=0.5; else if(s>1) for(int k=0;k<8;++k) st[k]=clamp01(y0[k]+0.14*sin((s+2)*(k+1)*1.37)); Fit z=optimize(st,cfg,obs,mu0,pi,delta); if(z.obj<best.obj) best=move(z); }
    P fit=decode(best.y,cfg);
    fs::create_directories("state");
    { ofstream f("state/clock.tsv"); f<<setprecision(17); for(size_t i=0;i<clock.size();++i) f<<"a"<<i<<' '<<clock[i]<<'\n'; }
    { ofstream f("state/gains.tsv"); f<<setprecision(17); for(auto&s:stations) f<<s.id<<'\t'<<s.gain.real()<<'\t'<<s.gain.imag()<<'\t'<<s.pole_tau<<'\t'<<s.zero_tau<<'\t'<<s.cal_rms<<'\n'; }
    { ofstream f("state/model.tsv"); f<<setprecision(17); f<<"rho0 "<<fit.rho[0]<<"\nrho1 "<<fit.rho[1]<<"\nrho2 "<<fit.rho[2]<<"\nh0 "<<fit.h[0]<<"\nh1 "<<fit.h[1]<<"\nm "<<fit.m<<"\ntau "<<fit.tau<<"\nc "<<fit.c<<"\nobjective "<<best.obj<<'\n'; }
    { ofstream f("state/predictions.tsv"); f<<setprecision(17); for(auto&o:obs){C z=forward_z(o.f,fit,cfg.disp,mu0,pi); double la=log10(abs(z)), ph=atan2(z.imag(),z.real())*180/pi; f<<o.id<<'\t'<<la<<'\t'<<ph<<'\t'<<la-o.la<<'\t'<<modphase(ph,o.ph)<<'\n';} }
    { ofstream f("answers.tsv"); f<<setprecision(17); for(auto&q:queries){C z=forward_z(q.f,fit,cfg.disp,mu0,pi); f<<q.id<<'\t'<<log10(abs(z))<<'\t'<<atan2(z.imag(),z.real())*180/pi<<'\n';} }
    return 0;
}
