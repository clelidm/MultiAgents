#include <cstdlib>	//for string
#include <string>
using namespace std;   // for string

//parameters
const int N = 1000;
const double eta=0.1; //, srml; /* srml = sqrt(m^2 - lambda^2) */
const double c=-1.; //cost, has to be <0
const double epsi=0.;// epsi \in [0.,1.]
const string dir= "./191-195_N"; //"./../now_Nav1_N";
const string graph_name = "graph.gexf";

//const double qq_a=1.;
//const double qq_b=1.;

//number of actions
const int min_i_Nit=1;  // minimum 1, maximum max_Nit
const int max_i_Nit=9;  // minimum 1, maximum 18 : nb_actions=10**(max_i_Nit)  --> for the moment, not larger than 10 for other reasons...
const int64_t max_Nit = 1e6;	//max value 9e18 -- see INT64_MAX
//start registrering
const int i_register=6; // it was 6

//average
const int N_average=100;

//histo
const double histo_x_min=0;
const double histo_x_max=100;
const int histo_N_steps=10*(int)(histo_x_max-histo_x_min);
const double histo_Dx=(histo_x_max - histo_x_min)*1./(histo_N_steps);

/******************************************************************************/
/************************ AGENTS  AND FOOD ************************************/
/******************************************************************************/
typedef struct
{
    int name;           //m = name of the ressource
    bool free;          //true if free
    double u;           //payoff of the ressource
    int nb_visit;       //total nb of visits of the ressource (m was empty)
    int nb_try;         //total nb of try of the ressource
    int A_name;         //name of last agent visiting
    int P_nb_elements;  //nb of different visits for computing the mean persistence
    double arrival_time;          //time at which I arrive at the resource
    double leaving_time;
    double T_occ;       //the time the resource has spent being occupied
} food;

typedef struct
{
    int name;
    bool home;      // true if agent is at home
    int *pos;        //food position (m), if out
                     //last food visited, if home
    double *Rew;     //total reward since the beginning/or since register activated
    double *reward; //reward since the beginning for each ressource i    
    int *nb_visit;  //nb of visit for each ressource m - the N-eme value is the total nb of visits
    int *nb_try;
    //double *reward_one_strategy;
    double *reward_diff_strategy;
} agent;

/******************************************************************************/
/********************************* RESULTS ************************************/
/******************************************************************************/
typedef struct
{
    int tries;
    int visits;
} system_;

typedef struct
{
    int N_av;
    int *nb_try; //nb for each spot m
    int *nb_visit; //nb for each spot m
    int *eff_nb_try; //nb for each spot m
    int *eff_nb_visit; //nb for each spot m
    double *e;
    double *e_real;
    double *Pi;
    double *P;
    double *Q;
    //double *nb_A;
}results;

typedef struct
{
    int N_samples;
    int *histo_try; //histo effective nb of try
    int *histo_visit;
    int *histo_try_m; //histo effective nb of try
    int *histo_visit_m;
    double mean_try;
    double mean_visit;
    double m2_visit;
    double m2_try;
    double mean_try_m;
    double mean_visit_m;
    double m2_visit_m;
    double m2_try_m;
}histo;
/******************************************************************************/
/********************************* GRAPH **************************************/
/******************************************************************************/
typedef struct
{
    int id=0;
    //bool oriented=false;
    int source=N;
    int target=N;
    double weight=0.;
} edge;
