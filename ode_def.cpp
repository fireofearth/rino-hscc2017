
/* ============================================================================
 File   : ode_def.cpp
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
The place to define initial conditions and parameters the systems of ODEs or DDEs on which to perform reachability
============================================================================ */

#include <assert.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>

#include "filib_interval.h"
#include "tadiff.h" 
#include "fadiff.h"
#include "fadbad_aa.h"
#include "ode_def.h"
#include "matrix.h"

using namespace std;

int sysdim; // dimension of system of ODE/DDE
int jacdim;  //  Jacobian will be dimension sysdim * jacdim

// parameters of the system of the ODEs
vector<AAF> inputs; // uncertain inputs and parameters : some will be used in initial condition, some as uncertain parameters
vector<AAF> center_inputs;
vector<interval> eps;

// for subdivisions of the initial domain to refine precision
int nb_subdiv_init; // number of subdivisiions
double recovering; // percentage of recovering between subdivisions
vector<vector<vector<interval>>> Xouter_print, Xouter_robust_print, Xouter_minimal_print, Xinner_print, Xinner_robust_print, Xinner_minimal_print, Xexact_print; // store results of subdivision
vector<double> t_print; // times where results are stored
int current_subdiv;
int current_iteration;

// for robust inner-approximations
int uncontrolled; // number of uncontrolled parameters (forall params)
int controlled; // number of controlled parameters (forall params)
vector<bool> is_uncontrolled; // for each input, uncontrolled or controlled (for robust inner-approx)
vector<bool> is_initialcondition; // for each input, initial condition or parameter (for robust outer-approx)
int variable; // number of non constant parameters
vector<bool> is_variable;  // for each parameter, constant or variable

// for ODEs : initialize the state variable (and center for inner-approximation)
void set_initialconditions(vector<AAF> &x, vector<AAF> &xcenter, vector<vector<AAF>> &J)
{
    // par défaut
    for (int i=0 ; i<sysdim ; i++) {
        x[i] = inputs[i];
        xcenter[i] = center_inputs[i];
    }
    setId(J);
}

// for ODEs and DDEs: define bounds for parameters and inputs, value of delay d0 if any, and parameters of integration (timestep, order of TM)
void init_system(OdeFunc &odef, double &t_begin, double &t_end, double &tau, double &d0, int &nb_subdiv, int &order) {

    /**
     * Colin: problem dimensions
     */
    struct OdeInit odeinit = odef.get_initial_variables();
    sysdim = odeinit.sysdim;
    jacdim = odeinit.jacdim;
    nb_subdiv_init = odeinit.nb_subdiv_init;
    t_begin = odeinit.t_begin;
    t_end = odeinit.t_end;
    tau = odeinit.tau;
    order = odeinit.order;
    inputs = odeinit.inputs;

    interval temp;
    int nb_points;
    
    uncontrolled = 0;
    controlled = 0;
    is_uncontrolled = vector<bool>(jacdim);
    variable = 0;
    is_variable = vector<bool>(jacdim);
    is_initialcondition = vector<bool>(jacdim);
    for (int i=0 ; i<jacdim; i++) {
        is_uncontrolled[i] = false;
        is_variable[i] = false;
        is_initialcondition[i] = false; // by definition, initial conditions are controlled
    }

    nb_points = (t_end-t_begin)/tau+1;
    
    // common to EDO and DDE
    /**
     *
     */
    center_inputs = vector<AAF>(jacdim);
    eps = vector<interval>(jacdim);
    for (int i=0 ; i<jacdim ; i++)
    {
        if (is_uncontrolled[i])
            uncontrolled ++;
        if (!is_uncontrolled[i] && !is_initialcondition[i])
            controlled++;
        temp = inputs[i].convert_int();
        center_inputs[i] = mid(temp);
        eps[i] = temp-mid(temp);
    }
    
    cout << "controlled=" << controlled  << " uncontrolled=" << uncontrolled << endl;
    
    // for saving results
  //  cout << "(t_end-t_begin)*nb_subdiv/d0+1=" << ((t_end-t_begin)/d0+1)*(nb_subdiv+1) << endl;
    Xouter_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    Xouter_robust_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    Xouter_minimal_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    Xinner_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    Xinner_robust_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    Xinner_minimal_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    Xexact_print = vector<vector<vector<interval>>>(nb_subdiv_init+1,vector<vector<interval>>(nb_points, vector<interval>(sysdim)));
    t_print = vector<double>(nb_points);
        
}

void init_subdiv(int current_subdiv, vector<AAF> inputs_save, int param_to_subdivide)
{
    center_inputs = vector<AAF>(jacdim);
    eps = vector<interval>(jacdim);
    
    interval save = inputs_save[param_to_subdivide].convert_int();
    double delta = (save.sup()-save.inf())/nb_subdiv_init;
    if ((current_subdiv > 1) && (current_subdiv < nb_subdiv_init))
        inputs[param_to_subdivide] = interval(
                save.inf()+delta*(current_subdiv-1-recovering),
                save.inf()+delta*(current_subdiv+recovering
                ));

    else if (current_subdiv == 1)
        inputs[param_to_subdivide] = interval(
                save.inf()+delta*(current_subdiv-1),
                save.inf()+delta*(current_subdiv+recovering
                ));
    else if (current_subdiv == nb_subdiv_init)
        inputs[param_to_subdivide] = interval(
                save.inf()+delta*(current_subdiv-1-recovering),
                save.inf()+delta*(current_subdiv
                ));
    cout << "inputs[param_to_subdivide] " << inputs[param_to_subdivide] << endl;
    
   
     interval   temp = inputs[param_to_subdivide].convert_int();
        center_inputs[param_to_subdivide] = mid(temp);
        eps[param_to_subdivide] = temp-mid(temp);
    
}
