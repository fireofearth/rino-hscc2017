/* ============================================================================
 File   : main.cpp
 Author : Sylvie Putot, Ecole Polytechnique (France)
 
 Part of the RINO package for Inner and Outer Reachability Analysis.
 
 This is the main file:
    - the system to analyze is chosen (among those defined in file ode_def.h/cpp)
    - subdivision of the input domain may be asked
    - the actual reachability analyis on the system of ODE or DDE is called
    - results are output (output files + gnuplot script)
 ============================================================================ */

#include "filib_interval.h"
#include "tadiff.h" 
#include "fadiff.h"
#include "fadbad_aa.h"
#include "utils.h"
#include "ode_def.h"
#include "matrix.h"
#include "inner.h"
#include "ode_integr.h"
#include <iostream>
#include <ostream>
#include <fstream>
#include <ctime>
#include <array>

#include <stdlib.h>

using namespace std;


vector<ofstream> outFile_outer_minimal;   //  minimal outer-approximated range for each variable of the system
vector<ofstream> outFile_outer_robust;  // robust outer-approximated range for each variable of the system
vector<ofstream> outFile_outer;   //  maximal outer-approximated range for each variable of the system
vector<ofstream> outFile_inner_minimal;   //  minimal inner-approximated range for each variable of the system
vector<ofstream> outFile_inner;   //  maximal inner-approximated range for each variable of the system
vector<ofstream> outFile_inner_robust;   // robust inner-approximated range for each variable of the system
vector<ofstream> outFile_center;

ofstream outFile_width_ratio;     //  min on xi ( width of inner-approx (xi) / width of outer-approx (xi) )
ofstream outFile_meanerror_outer; // mean on xi of error between outer-approx and analytical solution if any
ofstream outFile_meanerror_inner; // mean on xi of error between inner-approx and analytical solution if any
ofstream outFile_meanerror_diff;  // mean on xi of error between outer-approx and inner-approx
ofstream outFile_relmeanerror_outer; // mean on xi of error between outer-approx and analytical solution if any, over width of exact tube
ofstream outFile_relmeanerror_inner; // mean on xi of error between inner-approx and analytical solution if any, over width of exact tube
ofstream outFile_relmeanerror_diff;  // mean on xi of error between outer-approx and inner-approx, over width of over-approx tube

void open_outputfiles();
void print_finalstats(clock_t begin);
void generate_gnuplot_script();

void print_initstats(vector<AAF> &x);

void print_ErrorMeasures(int current_iteration, vector<AAF> inputs_save, double d0);
void print_finalsolution(vector<AAF> inputs_save, int max_it, double d0);

int main(int argc, char* argv[]) {
    // all these parameters in initialization functions defined in ode_def.cpp
    double tn;    // current time
    double tau;   // integration time step (fixed step for now)
    double t_begin; // starting time of initialization
    double t_end; // ending time of integration
    int order;    // order of Taylor expansion
    
    double d0; // = 1;   // delay in DDE
    int nb_subdiv; // = 10;   // number of Taylor models on [0,d0]
    OdeFunc odef;

    if (argc == 2) {
        int syschoice = atoi(argv[1]);
        odef = OdeFunc(static_cast<Problem>(syschoice));
    }
    else {
        cout << "choose number...\n";
        return 0;
    }

    clock_t begin = clock();

    init_system(odef, t_begin, t_end, tau, d0, nb_subdiv, order);

    open_outputfiles();

    // Colin: seems to be used to preserve inputs
    vector<AAF> inputs_save(jacdim);
    for (int i=0 ; i<jacdim; i++)
        inputs_save[i] = inputs[i];

    vector<vector<AAF>> J(sysdim, vector<AAF>(jacdim));
    vector<AAF> x(sysdim);
    vector<AAF> xcenter(sysdim);

    for (current_subdiv=1; current_subdiv <= nb_subdiv_init; current_subdiv++) {
        if (nb_subdiv_init > 1) init_subdiv(current_subdiv, inputs_save, 0);
        current_iteration = 0;

        cout << "center_inputs[0]" << center_inputs[0] << endl;
        cout << "center_inputs[1]" <<  center_inputs[1] << endl;
        cout << inputs[0] << endl;
        cout << inputs[1] << endl;

        set_initialconditions(x,xcenter,J);  //            setId(J0);

        tn = t_begin;
        print_initstats(inputs);

        HybridStep_ode cur_step = init_ode(odef,xcenter,x,J,tn,tau,order);

        while (cur_step.tn < t_end) {
            // build Taylor model for Value and Jacobian and deduce guards for each active mode
            cur_step.TM_build();
            cur_step.TM_evalandprint_solutionstep(eps,cur_step.tn+tau);
            cur_step.init_nextstep(tau);
        }
    }
    print_finalsolution(inputs_save, (t_end-t_begin)/tau, d0);
    
    print_finalstats(begin);

    generate_gnuplot_script();
}

// printig solution from all stored results of subdivisiions
void print_finalsolution(vector<AAF> inputs_save, int max_it, double d0)
{
    // print final solution + output some stats / error information
    bool no_hole = true;
    
    for (current_iteration = 0; current_iteration <= max_it; current_iteration++)
    {
        for (int i=0 ; i<sysdim ; i++)
        {
            Xouter_print[0][current_iteration][i] = Xouter_print[1][current_iteration][i];
            Xouter_robust_print[0][current_iteration][i] = Xouter_robust_print[1][current_iteration][i];
            Xinner_print[0][current_iteration][i] = Xinner_print[1][current_iteration][i];
            Xinner_robust_print[0][current_iteration][i] = Xinner_robust_print[1][current_iteration][i];
            for (int j=2; j <=nb_subdiv_init; j++)
            {
                Xouter_print[0][current_iteration][i] = hull(Xouter_print[0][current_iteration][i],Xouter_print[j][current_iteration][i]);
                Xouter_robust_print[0][current_iteration][i] = hull(Xouter_robust_print[0][current_iteration][i],Xouter_robust_print[j][current_iteration][i]);
                // verifying there is no hole in the inner-approx
                if (!((Xinner_print[j][current_iteration][i].sup() >= Xinner_print[0][current_iteration][i].inf()) &&
                      (Xinner_print[0][current_iteration][i].sup() >= Xinner_print[j][current_iteration][i].inf())))
                {
                    no_hole = false;
                }
                // Warning: joining tubes is correct for inner approx only if no hole
                Xinner_print[0][current_iteration][i] = hull(Xinner_print[0][current_iteration][i],Xinner_print[j][current_iteration][i]);
                Xinner_robust_print[0][current_iteration][i] = hull(Xinner_robust_print[0][current_iteration][i],Xinner_robust_print[j][current_iteration][i]);
            }
            if (nb_subdiv_init > 1)
                outFile_outer[i] << t_print[current_iteration] << "\t" << inf(Xouter_print[0][current_iteration][i]) << "\t" << sup(Xouter_print[0][current_iteration][i]) << endl;
        }
        
        print_ErrorMeasures(current_iteration,inputs_save,d0);
    }
    
    if (no_hole)
        cout << "NO HOLE when joining the inner-approx tubes";
}

void print_ErrorMeasures(int current_iteration, vector<AAF> inputs_save, double d0) {
    double aux, minwidth_ratio, sum, rel_sum;
    vector<interval> Xexact(sysdim);
    
    
    // correct only of no hole
    // min on xi of ratio of width of inner / width of outer
    minwidth_ratio = (sup(Xinner_print[0][current_iteration][0])-inf(Xinner_print[0][current_iteration][0]))/(sup(Xouter_print[0][current_iteration][0])-inf(Xouter_print[0][current_iteration][0]));
    for (int i=1 ; i<sysdim ; i++) {
        aux = (sup(Xinner_print[0][current_iteration][i])-inf(Xinner_print[0][current_iteration][i]))/(sup(Xouter_print[0][current_iteration][i])-inf(Xouter_print[0][current_iteration][i]));
        if (minwidth_ratio > aux)
            minwidth_ratio = aux;
    }
    if (t_print[current_iteration] != 0)
    outFile_width_ratio << t_print[current_iteration] << "\t" << minwidth_ratio << endl;
    
    // mean over the xi of the error between over-approx and inner-approx
    sum = 0;
    rel_sum = 0;
    for (int i=0 ; i<sysdim ; i++)
    {
        aux = max(sup(Xouter_print[0][current_iteration][i])-sup(Xinner_print[0][current_iteration][i]),inf(Xinner_print[0][current_iteration][i])-inf(Xouter_print[0][current_iteration][i]));
        sum += aux;
        rel_sum += aux / (sup(Xouter_print[0][current_iteration][i])-inf(Xouter_print[0][current_iteration][i]));
    }
    sum = sum/sysdim;
    rel_sum = rel_sum/sysdim;
    outFile_meanerror_diff << t_print[current_iteration] << "\t" << sum << endl;
    outFile_relmeanerror_diff << t_print[current_iteration] << "\t" << rel_sum << endl;
}

void read_parameters(const char * params_filename, double &tau, double &t_end, int &order, char * sys_name, char *initial_condition)
{
    const int LINESZ = 2048;
    char buff[LINESZ];
    char output_variables[1000];
    
    cout << "****** Reading system parameter from cfg file ******" << endl;
    sprintf(sys_name,"sys"); // default name of main system to read from the xml file
    //  ifstream inputFileStream("examples/ex_Hyst/brusselator/brusselator.cfg");
    
    FILE *params_file = fopen(params_filename,"r");
    if (params_file == NULL)
        cout << "Error reading " << params_filename << ": file not found" << endl;
    while (fgets(buff,LINESZ,params_file)) {
        sscanf(buff, "system = %s\n", sys_name);
        sscanf(buff, "initially = %[^\n]\n", initial_condition);   // tell separator is newline, otherwise by default it is white space
        sscanf(buff, "time-horizon = %lf\n", &t_end);
        sscanf(buff, "sampling-time = %lf\n", &tau);
        sscanf(buff, "output-variables = %[^\n]\n", output_variables);
        sscanf(buff, "order = %d\n", &order);
    }
    fclose(params_file);
    cout << "system name = " << sys_name << endl;
    cout << "initial condition = " << initial_condition << endl;
    cout << "tau = " << tau << endl;
    cout << "t_end = " << t_end << endl;
    cout << "output_variables = " << output_variables << endl;
    cout << "order = " << order << endl;
    cout << "****** End parameter reading ******" << endl << endl;
}

void generate_gnuplot_script()
{
    ofstream gnuplot_script;
    gnuplot_script.open("output/gnuplot_script.gp");
    
    // ---------- plotting variables ----------
    int nb_lignes, nb_colonnes;
    nb_lignes = sqrt(sysdim);
    nb_colonnes = nb_lignes;
    if (nb_lignes * nb_colonnes < sysdim)
        nb_lignes++;
    if (nb_lignes * nb_colonnes < sysdim)
        nb_colonnes++;
    
    // plotting indiviudually each solution to a file
    gnuplot_script << "set term pngcairo font \"Helvetica,18\"" << endl;
    
    for (int i=0 ; i<sysdim ; i++) {
        if (sysdim == 1)
            gnuplot_script << "set output \"x.png\"" << endl;
        else
            gnuplot_script << "set output \"x"<<i+1<<".png\"" << endl;

        if (nb_subdiv_init > 1)
            gnuplot_script << "set title '"<<nb_subdiv_init<<" subdivisions'" << endl;
        gnuplot_script << "set xlabel 't (seconds)'" << endl;

        gnuplot_script << "set ylabel 'x(t)'" << endl;

        if ((uncontrolled > 0)) // && (controlled > 0))
        {
            gnuplot_script << "set style fill noborder"<<endl;
            gnuplot_script << "plot 'x"<<i+1<<"outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  title \"maximal outer flowpipe\", 'x"<<i+1<<"outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, ";
            gnuplot_script << "'x"<<i+1<<"inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  title \"maximal inner flowpipe\", ";
            gnuplot_script << "'x"<<i+1<<"outer_robust.out' using 1:2 w l lt 5  lc rgb \"red\"  lw 2 title \"robust outer flowpipe\", 'x"<<i+1<<"outer_robust.out' using 1:3 w l lt 5  lc rgb \"red\"  lw 2notitle,";
            gnuplot_script << "'x"<<i+1<<"inner_robust.out' using 1:2:3 w filledcu lc rgb \"red\" title \"robust inner flowpipe\", ";
            gnuplot_script << "'x"<<i+1<<"outer_minimal.out' using 1:2 w l lt 7 lc rgb '#7e2f8e' lw 2 title \"minimal outer flowpipe\", 'x"<<i+1<<"outer_minimal.out' using 1:3 w l lt 7 lc rgb '#7e2f8e' lw 2 notitle,";
            gnuplot_script << "'x"<<i+1<<"inner_minimal.out' using 1:2:3 w filledcu lc rgb '#7e2f8e' title \"minimal inner flowpipe\" " << endl;
        }
        else if (controlled > 0)
        {
            gnuplot_script << "set style fill noborder"<<endl;
            gnuplot_script << "plot 'x"<<i+1<<"outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  title \"maximal outer flowpipe\", 'x"<<i+1<<"outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, ";
            gnuplot_script << "'x"<<i+1<<"inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  title \"maximal inner flowpipe\", ";
            gnuplot_script << "'x"<<i+1<<"outer_minimal.out' using 1:2 w l lt 7 lc rgb '#7e2f8e' lw 2 title \"minimal outer flowpipe\", 'x"<<i+1<<"outer_minimal.out' using 1:3 w l lt 7 lc rgb '#7e2f8e' lw 2 notitle,";
            gnuplot_script << "'x"<<i+1<<"inner_minimal.out' using 1:2:3 w filledcu lc rgb '#7e2f8e' title \"minimal inner flowpipe\" " << endl;
            
        }
        else // no uncertain parameters, only innitial conditions
        {
            gnuplot_script << "set style fill noborder"<<endl;
            gnuplot_script << "plot 'x"<<i+1<<"outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  title \"maximal outer flowpipe\", 'x"<<i+1<<"outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, ";
            gnuplot_script << "'x"<<i+1<<"inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  title \"maximal inner flowpipe\" " << endl;
            
        }
        gnuplot_script << "unset output" << endl;
    }
    
    gnuplot_script << "set term pngcairo font \"Helvetica,18\"" << endl;
    gnuplot_script << "set output \"xi.png\"" << endl;
    if (nb_subdiv_init > 1)
        gnuplot_script << "set title '"<<nb_subdiv_init<<" subdivisions'" << endl;

    gnuplot_script << "set xlabel 't (seconds)'" << endl;
    
    if (uncontrolled > 0)
    {
        gnuplot_script << "set style fill noborder"<<endl;
        
        gnuplot_script << "plot 'x1outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  title \"maximal outer flowpipe\", 'x1outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, ";
        gnuplot_script << "'x1inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  title \"maximal inner flowpipe\", ";
        gnuplot_script << "'x1outer_robust.out' using 1:2 w l lt 5  lc rgb \"red\"  lw 2 title \"robust outer flowpipe\", 'x1outer_robust.out' using 1:3 w l lt 5  lc rgb \"red\"  lw 2 notitle,";
        gnuplot_script << "'x1inner_robust.out' using 1:2:3 w filledcu lc rgb \"red\" title \"robust inner flowpipe \", ";
        gnuplot_script << "'x1outer_minimal.out' using 1:2 w l lt 7 lc rgb '#7e2f8e' lw 2 title \"minimal outer flowpipe\", 'x1outer_minimal.out' using 1:3 w l lt 7 lc rgb '#7e2f8e' lw 2 notitle,";
        gnuplot_script << "'x1inner_minimal.out' using 1:2:3 w filledcu lc rgb '#7e2f8e' title \"minimal inner flowpipe\", ";
        
        gnuplot_script << " 'x2outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, 'x2outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, ";
        gnuplot_script << "'x2inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  notitle , ";
        gnuplot_script << "'x2outer_robust.out' using 1:2 w l lt 5  lc rgb \"red\"  lw 2 notitle, 'x2outer_robust.out' using 1:3 w l lt 5  lc rgb \"red\"  lw 2 notitle,";
        gnuplot_script << "'x2inner_robust.out' using 1:2:3 w filledcu lc rgb \"red\" notitle, ";
        gnuplot_script << "'x2outer_minimal.out' using 1:2 w l lt 7 lc rgb '#7e2f8e' lw 2 notitle, 'x2outer_minimal.out' using 1:3 w l lt 7 lc rgb '#7e2f8e' lw 2 notitle,";
        gnuplot_script << "'x2inner_minimal.out' using 1:2:3 w filledcu lc rgb '#7e2f8e' notitle " << endl;
    }
    else if (controlled > 0)
    {
        gnuplot_script << "set style fill noborder"<<endl;
        
        gnuplot_script << "plot 'x1outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  title \"maximal outer flowpipe \", 'x1outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, ";
        gnuplot_script << "'x1inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  title \"maximal inner flowpipe \", ";
        gnuplot_script << "'x1outer_minimal.out' using 1:2 w l lt 7 lc rgb '#7e2f8e' lw 2 title \"minimal outer flowpipe \", 'x1outer_minimal.out' using 1:3 w l lt 7 lc rgb '#7e2f8e' lw 2 notitle,";
        gnuplot_script << "'x1inner_minimal.out' using 1:2:3 w filledcu lc rgb '#7e2f8e' title \"minimal inner flowpipe \", ";
        
        gnuplot_script << " 'x2outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, 'x2outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, ";
        gnuplot_script << "'x2inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  notitle, ";
        gnuplot_script << "'x2outer_minimal.out' using 1:2 w l lt 7 lc rgb '#7e2f8e' lw 2 notitle, 'x2outer_minimal.out' using 1:3 w l lt 7 lc rgb '#7e2f8e' lw 2 notitle,";
        gnuplot_script << "'x2inner_minimal.out' using 1:2:3 w filledcu lc rgb '#7e2f8e' notitle " << endl;
    }
    else // only initial conditions => only maximal flowpipes
    {
        gnuplot_script << "set style fill noborder"<<endl;
        
        gnuplot_script << "plot 'x1outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  title \"maximal outer flowpipe \", 'x1outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, ";
        gnuplot_script << "'x1inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  title \"maximal inner flowpipe \", ";
        
        gnuplot_script << " 'x2outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, 'x2outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, ";
        gnuplot_script << "'x2inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  notitle " << endl;
    }
    gnuplot_script << "unset output" << endl;
    
    
    // plotting the width ration
    gnuplot_script << "set term pngcairo" << endl;
    gnuplot_script << "set output \"width_ratio.png\"" << endl;
    gnuplot_script << "set xlabel 't (seconds)'" << endl;
    gnuplot_script << "set ylabel 'min over x_i of width ratios'" << endl;
    gnuplot_script << "plot 'width_ratio.out' w l title \"width(inner-approx) / width (outer-approx)" << endl;
    gnuplot_script << "unset output" << endl;
    
    // plotting  mean on xi of error between outer-approx and analytical solution if any
    gnuplot_script << "set term pngcairo" << endl;
    gnuplot_script << "set output \"meanerror.png\"" << endl;
    gnuplot_script << "set xlabel 't (seconds)'" << endl;
    gnuplot_script << "set ylabel 'mean over x_i of max error'" << endl;
    if (nb_subdiv_init > 1)
        gnuplot_script << "set title '"<<nb_subdiv_init<<" subdivisions'" << endl;
    gnuplot_script << "plot 'meanerror_outer.out' w l lt 1 title \"outer-approx error\", 'meanerror_inner.out' w l lt 1 dashtype 2 title \"inner-approx error\", 'meanerror_diff.out' w l lt 2 title \"max distance between inner and outer-approx\"" << endl;
    gnuplot_script << "unset output" << endl;
    
    // plotting  mean on xi of relative error between outer-approx and analytical solution if any (relative because divided by width of exact or over-approx solution)
    gnuplot_script << "set term pngcairo" << endl;
    gnuplot_script << "set output \"relmeanerror.png\"" << endl;
    gnuplot_script << "set xlabel 't (seconds)'" << endl;
    gnuplot_script << "set ylabel 'mean over x_i of max relative error'" << endl;
    if (nb_subdiv_init > 1)
        gnuplot_script << "set title '"<<nb_subdiv_init<<" subdivisions'" << endl;
    gnuplot_script << "plot 'relmeanerror_outer.out' w l lt 1 title \"outer-approx error\", 'relmeanerror_inner.out' w l lt 1 dashtype 2 title \"inner-approx error\", 'relmeanerror_diff.out' w l lt 2 title \"max distance between inner and outer-approx\"" << endl;
    gnuplot_script << "unset output" << endl;
    
    gnuplot_script << "set terminal aqua" << endl;
    
    gnuplot_script.close();
    
    cout << "......" << endl;
    cout << "Result files are in the output directory" << endl;
    cout << "Visualize them with gnuplot by : cd output; gnuplot -p 'gnuplot_script.gp' ; cd .." << endl;
    cout << "The gnuplot script will also save files as png files (in same directory)" << endl;
    system("cd output; gnuplot -p 'gnuplot_script.gp' ; cd ..");
}



void open_outputfiles() {
    system("rm -r output");
    system("mkdir output");
    
    outFile_outer = vector<ofstream>(sysdim);   // output outer-approximated range for each variable of the system
    outFile_outer_robust = vector<ofstream>(sysdim);
    outFile_outer_minimal = vector<ofstream>(sysdim);
    outFile_inner = vector<ofstream>(sysdim); // output inner-approximated range for each variable of the system
    outFile_inner_robust = vector<ofstream>(sysdim);
    outFile_inner_minimal = vector<ofstream>(sysdim);
    outFile_center  = vector<ofstream>(sysdim);
    
    stringstream file_name;
   
    for (int i=0 ; i<sysdim ; i++) {
        file_name.str("");
        file_name << "output/x" << i+1 << "outer.out";
        outFile_outer[i].open(file_name.str().c_str());
        file_name.str("");
        file_name << "output/x" << i+1 << "outer_robust.out";
        outFile_outer_robust[i].open(file_name.str().c_str());
        file_name.str("");
        file_name << "output/x" << i+1 << "outer_minimal.out";
        outFile_outer_minimal[i].open(file_name.str().c_str());
        file_name.str("");
        file_name << "output/x" << i+1 << "center.out";
        outFile_center[i].open(file_name.str().c_str());
        file_name.str("");
        file_name << "output/x" << i+1 << "inner.out";
        outFile_inner[i].open(file_name.str().c_str());
        file_name.str("");
        file_name << "output/x" << i+1 << "inner_robust.out";
        outFile_inner_robust[i].open(file_name.str().c_str());
        file_name.str("");
        file_name << "output/x" << i+1 << "inner_minimal.out";
        outFile_inner_minimal[i].open(file_name.str().c_str());
    }

    outFile_width_ratio.open("output/width_ratio.out");
    outFile_meanerror_outer.open("output/meanerror_outer.out");
    outFile_meanerror_inner.open("output/meanerror_inner.out");
    outFile_meanerror_diff.open("output/meanerror_diff.out");
    outFile_relmeanerror_outer.open("output/relmeanerror_outer.out");
    outFile_relmeanerror_inner.open("output/relmeanerror_inner.out");
    outFile_relmeanerror_diff.open("output/relmeanerror_diff.out");
}



void print_initstats(vector<AAF> &x)
{
   
    // print initial conditions of the ODE
    for (int i=0 ; i<sysdim ; i++) {
        // print in exit files
        outFile_outer[i] << 0 << "\t" << inf(x[i].convert_int()) << "\t" << sup(x[i].convert_int()) << endl;
        outFile_outer_robust[i] << 0 << "\t" << inf(x[i].convert_int()) << "\t" << sup(x[i].convert_int()) << endl;
        outFile_outer_minimal[i] << 0 << "\t" << inf(x[i].convert_int()) << "\t" << sup(x[i].convert_int()) << endl;
        outFile_center[i] << 0<< "\t" << mid(x[i].convert_int()) << "\t" << mid(x[i].convert_int()) << endl;
        outFile_inner[i] << 0 << "\t" << inf(x[i].convert_int()) << "\t" << sup(x[i].convert_int()) << endl;
        outFile_inner_robust[i] << 0 << "\t" << inf(x[i].convert_int()) << "\t" << sup(x[i].convert_int()) << endl;
        outFile_inner_minimal[i] << 0 << "\t" << inf(x[i].convert_int()) << "\t" << sup(x[i].convert_int()) << endl;
    }
    outFile_width_ratio << 0 << "\t" << 1.0 << endl;
    
    cout << "printing at t=0 : x=" << endl;
    for (int i=0 ; i<sysdim ; i++)
        cout << "x[" << i <<"]=" << x[i] << "\t";
    cout << endl;
    cout << "x0=" << endl;
    for (int i=0 ; i<sysdim ; i++)
        cout << "x0[" << i <<"]=" << mid(x[i].convert_int()) << "\t";
    cout << endl;
}

// print after the end of the analysis
void print_finalstats(clock_t begin)
{
    clock_t end = clock();
    // double end_time = getTime ( );
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "elpased time (sec) =" << elapsed_secs << endl;
    
    
    for (int i=0 ; i<sysdim ; i++) {
        outFile_outer[i].close();
        outFile_outer_robust[i].close();
        outFile_outer_minimal[i].close();
        outFile_inner[i].close();
        outFile_inner_robust[i].close();
        outFile_inner_minimal[i].close();
        outFile_center[i].close();
    }
    outFile_width_ratio.close();
    outFile_meanerror_outer.close();
    outFile_meanerror_inner.close();
    outFile_meanerror_diff.close();
    outFile_relmeanerror_outer.close();
    outFile_relmeanerror_inner.close();
    outFile_relmeanerror_diff.close();
}

