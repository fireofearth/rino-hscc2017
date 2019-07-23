/* ============================================================================
 File   : ode_def.h
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 The place to declare the systems of ODEs or DDEs on which he/she wants to perform reachability
 Setup of dimension of the system (sysdim), initial conditions and parameters is in ode_def.cpp
 ============================================================================ */
#ifndef ODE_DEF_H
#define ODE_DEF_H

#include <exception>
#include "filib_interval.h"
#include "tadiff.h" 
#include "fadiff.h"
#include "fadbad_aa.h"

using namespace std;

#define innerapprox 1

enum Problem {NN, BRUSSELATOR, VANDERPOLL, CIRCLE, BOX, LORENZATTR, BALLISTIC, LOTKAVOLTERRA, FITHUGHNAGUMO};

// dimension of the system of ODE/DDE to analyze:
extern int sysdim;
// dimension of uncertain input
extern int jacdim; // Jacobian will be dimension sysdim * jacdim
extern vector<AAF> inputs;   // uncertain inputs and parameters : some will be used in initial condition, some as uncertain parameters
extern vector<AAF> center_inputs;
extern vector<interval> eps;

// for subdivisions of the initial domain to refine precision
extern int nb_subdiv_init; // number of subdivisiions
extern double recovering; // percentage of recovering between subdivisions
extern vector<vector<vector<interval>>> Xouter_print, Xouter_robust_print, Xouter_minimal_print, Xinner_print,
    Xinner_robust_print, Xinner_minimal_print, Xexact_print; // store results of subdivision
extern vector<double> t_print; // times where results are stored
extern int current_subdiv;
extern int current_iteration;

// for robust inner-approximations
extern int uncontrolled;  // number of uncontrolled parameters (forall params)
extern int controlled;  // number of controlled parameters (forall params)
extern vector<bool> is_uncontrolled; // for each input, uncontrolled or controlled (for robust inner-approx)
extern vector<bool> is_initialcondition; // for each input, initial condition or parameter (for robust inner-approx)
extern int variable;  // number of non constant parameters
extern vector<bool> is_variable; // for each parameter, constant or variable

struct OdeInit {
    int sysdim;
    int jacdim;
    int nb_subdiv_init;
    double tau;
    double t_begin;
    double t_end;
    int order;
    vector<AAF> inputs;
};

class OdeFunc {
private:
    Problem p;
public:
    OdeFunc(Problem p) {
        this-> p = p;
    }

    OdeFunc() {
        this->p = NN;
    }

    struct OdeInit get_initial_variables() {
        /**
         * Colin: problem specifications to initialize
         */
        switch(this->p) {
            case BRUSSELATOR:
                return { 2, 2, 1, 0.05, 0, 10., 4,
                         { interval(0.9, 1), interval(0, 0.1) }
                };
            case VANDERPOLL:
                return { 2, 2, 1, 0.05, 0., 10., 5,
                        { interval(2., 2.1), interval(0., 0.1) }
                };
            case CIRCLE:
                return { 2, 2, 1, 0.05, 0., 10., 4,
                        { interval(-1.6, -1.4), interval(1.4, 1.6) }
                };
            case BOX:
                return { 2, 2, 1, 0.05, 0., 10., 4,
                    { interval(0.10, 0.2), interval(0.10,0.2) }
                };
            case LORENZATTR:
                return { 3, 3, 1, 0.05, 0., 10., 4,
                    { interval(0.0,0.1), interval(0.0,0.1), interval(0.0,0.1) }
                };
            case BALLISTIC:
                return {4, 4, 1, 0.1, 0., 4., 3,
                        { interval(181.,185.), 3.14159/180*interval(2.5,3.5), interval(0.0,0.01), interval(0.0,0.01)}
                };
            case LOTKAVOLTERRA:
                return {2, 2, 1, 0.05, 0., 20., 4,
                    {interval(0.25, 0.28), interval(0.16, 0.19)}
                };
            case FITHUGHNAGUMO:
                return { 2, 2, 1, 0.05, 0., 25, 4,
                    {interval(1., 1.2), interval(2.25, 2.45)}
                };
            default:
                throw logic_error("problem is not specified");
        }
    }

    template <class C>
    void operator()(vector<C> &yp, vector<C> y) {
        double g = 9.81; // gravity in m/s^2
        double rho = 1204.4; // air density in g/m3
        double a = 0.000126677; // cross section du bullet d=12.7mm, cross section=pi R^2
        double d = 0.45; // drag coefficient
        double m = 14.3; // mass du bullet in g

        /**
         * Colin: problem formulas to switch to
         */
        switch(this->p) {
            case BRUSSELATOR: // OK
                yp[0] = 1 - (1.5 + 1) * y[0] + y[0] * y[0] * y[1];
                yp[1] = 1.5 * y[0] - y[0] * y[0] * y[1];
                break;
            case VANDERPOLL: // meh. Loses accuracy after t > 4.5
                yp[0] = 1.0 * (1 - y[1] * y[1]) * y[0] - y[1];
                yp[1] = y[0];
                break;
            case CIRCLE: // OK
                yp[0] = -y[1] + y[0] * (1 - y[0] * y[0] - y[1] * y[1]);
                yp[1] =  y[0] + y[1] * (1 - y[0] * y[0] - y[1] * y[1]);
                break;
            case BOX: // OK
                yp[0] = (y[1] + 0.2 * y[0]) * (1 - y[0] * y[0]);
                yp[1] = - y[0] * (1 - y[1] * y[1]);
                break;
            case LORENZATTR: // doesn't give anything
                yp[0] = 10. * (y[1] - y[0]);
                yp[1] = y[0] * (28. - y[2]) - y[1];
                yp[2] = y[0] * y[1] - 2.667 * y[2];
                break;
            case BALLISTIC: // from Goubault/Putot
                yp[0] = - g*sin(y[1])-rho*y[0]*y[0]*a*d/(2.0*m); // velocity v
                yp[1] = - g*cos(y[1])/y[0]; // angle gamma with respect to the x axis
                yp[2] = y[0]*cos(y[1]); // position x
                yp[3] = y[0]*sin(y[1]); // position y
                break;
            case LOTKAVOLTERRA: // OK
                yp[0] = (0.5 - y[1]) * y[0];
                yp[1] = (-0.5 + y[0]) * y[1];
                break;
            case FITHUGHNAGUMO:
                yp[0] = y[0] - 0.3333 * y[0] * y[0] * y[0] - y[1] + 0.875;
                yp[1] = 0.08 * (y[0] + 0.7 - 0.8 * y[1]);
                break;
            default:
                throw logic_error("problem is not specified");
        }
    }
};

// for ODEs : initialize the state variable (and center for inner-approximation)
void set_initialconditions(vector<AAF> &x, vector<AAF> &xcenter, vector<vector<AAF>> &J);

// for ODEs and DDEs: define bounds for parameters and inputs, value of delay d0 if any, and parameters of integration (timestep, order of TM)
void init_system(OdeFunc &odef, double &t_begin, double &t_end, double &tau, double &d0, int &nb_subdiv, int &order);

// specific to subdivisions
void init_subdiv(int current_subdiv, vector<AAF> inputs_save, int param_to_subdivide);

#endif