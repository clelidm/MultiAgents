#include <vector>
#include <string>
#include "data.h"
using namespace std;

/******************************************************************************/
/*********************************  PRINTS ************************************/
/******************************************************************************/
void print_data_food(food*, int, int, int);
void print_nb_visit(agent*, int);
string file_name2(int Nit, string name);
/******************************************************************************/
/********************************  STORE DATA  ********************************/
/******************************************************************************/
void store_data(results*, histo *, food*, agent *, system_, double);
void print_stored_data(results, histo, food*, int);
void print_stored_data_tab(results*, histo*, food*, int);
void print_kn_i_m(agent *, int Nit, int i);
void print_e_m(food *, int Nit, double t);
/******************************************************************************/
/*********************************  DISTRIB  **********************************/
/******************************************************************************/
void remplir(double *tab_distrib, double xmin, double xmax, double Dx, double xj);
/********************************************************************/
/**************************    GENERATOR    *************************/
/********************************************************************/
void initialise_rand();
double uni(double b=1, double a=0);
double rand_exp(double);
int rand_int(int b=100, int a=0);
/******************************************************************************/
/************************* TOOLS CONTAINERS ***********************************/
/******************************************************************************/
agent remove_agent(int, vector<agent>& vec);
/******************************************************************************/
/********************** INITIALISATION ****************************************/
/******************************************************************************/
double *init_tab(int);
food *init_food(int, double payoff_fn(int));
agent *init_agent(int, double payoff_fn(int));
results init_results(int);
histo init_histo(int);
results* init_results_tab(int N, int Nit);
histo* init_histo_tab(int N, int Nit);
void re_init_all(agent *A, food *F, system_ *total_nb, double payoff_fn(int));
/******************************************************************************/
/*************************** PERSISTENCE **************************************/
/******************************************************************************/
void persistence_up(food *, int, int, int);
/******************************************************************************/
/*********************  ACTIONS -- OUT & HOME  ********************************/
/******************************************************************************/
void go_out(int *H, vector<agent>& li_home, vector<agent>& li_out, food *F,
            int search_strategy(agent, food*, system_ *, double*), system_ *total_nb, double t, list<int>&);
void go_home(int *H, vector<agent>& li_home, vector<agent>& li_out, food *F, double t);
/******************************************************************************/
/***************************  GRAPH FUNCTIONS  ********************************/
/******************************************************************************/
void create_edges(list< pair<int, double> >& li_A, int nb_eff, list<edge>& li_E);
void create_graph(list<edge> &li_E, int);




