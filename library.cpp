#include <cstdlib>	//permet d'utiliser alea (et pas que)
#include <iostream> //cout
#include <cmath>
#include <ctime>
#include <vector>
#include <list>
#include "data.h"
using namespace std;

/********************************************************************/
/**************************    GENERATOR    *************************/
/********************************************************************/
void initialise_rand()  { srand48((unsigned)time(NULL)); }
double uni(double b=1, double a=0) { return a+(b-a)*drand48(); }
double rand_exp(double lambda) { return -log(uni())/lambda;}
int rand_int(int b=100, int a=0) { return (int) (uni()*(b-a));}
/******************************************************************************/
/******************************** PRIOR ***************************************/
/******************************************************************************/
/**** settlers ****/
const double q_try=10.;
double q_a_set_N1000(int i=1, int m=1) //N=1000
{
    if(i>=196 && i<=207)
      {
	if(m==i){return q_try;}
	else {return 1;}
      }
     else {return q_try;} 
}
double q_a_set_N1000_right(int i=1, int m=1) //N=1000
{
    if(i>=191 && i<=195)
    {
        if(m==i){return q_try;}
        else {return 1;}
    }
    else {return 1;}
}
double q_a_set_bef(int i=1, int m=1) //test N=500 before
{
    if(i>=45 && i<=50)
    {
        if(m==i){return q_try;}
        else {return 1;}
    }
    else {return q_try;}
}
double q_a_set_aft(int i=1, int m=1) //test N=500 after
{
    if(i==145)
    {
        if(m==i){return q_try;}
        else {return 1;}
    }
    else {return q_try;}  //no settlers
}
double q_b_set(int i=1, int m=1)
{
    return q_try+1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//choosen prior
double q_a(int i=1, int m=1) // i=agent_name, m=spot   
//has to be a>=1
{
    return q_a_set_N1000_right(i, m); //10;
}
double q_b(int i=1, int m=1)
//has to be a+b>=1
{
    return 11; //q_b_set(i, m);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/******************************************************************************/
/************************* TOOLS CONTAINERS ***********************************/
/******************************************************************************/
agent remove_agent(int i, vector<agent>& vec)
{
    agent Ai = *(vec.begin() + i);
    vec.erase(vec.begin() + i);
    return Ai;
}
/******************************************************************************/
/************************* TOOLS TABLES ***************************************/
/******************************************************************************/
//// fill tables
void init_rew(double *tab, double payoff_fn(int), int A_name)
{    for(int m=0; m<N; m++)  {tab[m]=payoff_fn(m) + c*q_b(A_name, m)/q_a(A_name, m);}   }

void zero_tab_int(int* tab, int N_size)
{   for(int i=0; i<N_size; i++)  {tab[i]=0;}   }

void zero_tab_double(double* tab)
{   for(int i=0; i<N; i++)  {tab[i]=0;}   }

//// create tables
double *init_tab(int N)
{
    double *u = (double *)malloc(N*sizeof(double));
    for(int i=0; i<N; i++)  {u[i]=0;}
    return u;
}
int *init_tab_int(int N)
{
    int *u = (int *)malloc(N*sizeof(int));
    for(int i=0; i<N; i++)  {u[i]=0;}
    return u;
}
 
double *init_tab_rew(double payoff_fn(int), int A_name)
{
    double *u = (double *)malloc(N*sizeof(double));    init_rew(u, payoff_fn, A_name);
    return u;
}
/******************************************************************************/
/********************** INITIALISATION ****************************************/
/******************************************************************************/
food *init_food(int N, double payoff_fn(int))
{
    food *F = (food *)malloc(N*sizeof(food));
    for(int i=0; i<N; i++)
    {
        F[i].name=i;
        F[i].free=true;
        F[i].u=payoff_fn(i);
        F[i].nb_visit=0;
        F[i].nb_try=0;
        F[i].A_name=N; //no agent have been there already
        F[i].P_nb_elements=0;
        F[i].arrival_time=0.;
        F[i].leaving_time=0.;
        F[i].T_occ=0.;	//total_time spot is occupied since the start
    }
    return F;
}


agent *init_agent(int N, double payoff_fn(int))
{
    agent *A=(agent *)malloc(N*sizeof(agent));
    for(int i=0; i<N; i++)
    {
        A[i].name=i;
        A[i].home=true;
        A[i].pos=(int *)malloc(sizeof(int));
            A[i].pos[0]=N;  //position of the last visited food, or current position //N=none of the food spot: 1rst time leaving home
        A[i].Rew=(double *)malloc(sizeof(double));
            A[i].Rew[0]=0;
        A[i].reward = init_tab(N);
        //A[i].reward_one_strategy = init_tab_rew(payoff_fn);
        A[i].reward_diff_strategy = init_tab_rew(payoff_fn, i);
        A[i].nb_visit = init_tab_int(N+1);
        A[i].nb_try = init_tab_int(N+1);
    }
    return A;
}

results init_results(int N)
{
    results R;
    R.N_av=0;
    R.nb_try=init_tab_int(N);
    R.nb_visit=init_tab_int(N);
    R.eff_nb_try=init_tab_int(N);   //effective nb of try at each patch m
    R.eff_nb_visit=init_tab_int(N); //effective nb of successfull try at each patch m
    R.e=init_tab(N);                //q(m) = probability that one agent going to m find the place empty
    R.e_real=init_tab(N);           //Real value of e(m) = pba that the site is empty
    R.Pi=init_tab(N);
    R.P=init_tab(N);
    R.Q=init_tab(N);    
    return R;
}

results* init_results_tab(int N, int Nit)
{
    results *R=(results *)malloc(Nit*sizeof(results));
    for(int i=0; i<Nit; i++)
    {
        R[i].N_av=0;
        R[i].nb_try=init_tab_int(N);
        R[i].nb_visit=init_tab_int(N);
        R[i].eff_nb_try=init_tab_int(N);
        R[i].eff_nb_visit=init_tab_int(N);
        R[i].e_real=init_tab(N);
        R[i].e=init_tab(N);
        R[i].Pi=init_tab(N);
        R[i].P=init_tab(N);
        R[i].Q=init_tab(N);
    }
    return R;
}

histo init_histo(int N)
{
    histo H;
    H.N_samples=0;
    //agent
    H.histo_try=init_tab_int(histo_N_steps); //histo effective nb of try
    H.histo_visit=init_tab_int(histo_N_steps);
    //food
    H.histo_try_m=init_tab_int(histo_N_steps); //histo effective nb of try
    H.histo_visit_m=init_tab_int(histo_N_steps);
    //agent
    H.mean_try=0;       //mean
    H.mean_visit=0;
    H.m2_visit=0;       //moment of order 2
    H.m2_try=0;
    //food
    H.mean_try_m=0;
    H.mean_visit_m=0;
    H.m2_visit_m=0;
    H.m2_try_m=0;
    return H;
}

histo* init_histo_tab(int N, int Nit)
{
    histo *H=(histo *)malloc(Nit*sizeof(histo));
    for(int i=0; i<Nit; i++)
    {
        H[i].N_samples=0;
        H[i].histo_try=init_tab_int(histo_N_steps); //histo effective nb of try
        H[i].histo_visit=init_tab_int(histo_N_steps);
        H[i].histo_try_m=init_tab_int(histo_N_steps); //histo effective nb of try
        H[i].histo_visit_m=init_tab_int(histo_N_steps);
        H[i].mean_try=0;
        H[i].mean_visit=0;
        H[i].m2_visit=0;
        H[i].m2_try=0;
    }
    return H;
}

void re_init_all(agent *A, food *F, system_ *total_nb, double payoff_fn(int))
{
    (*total_nb).visits=0;     (*total_nb).tries=0;
    for(int i=0; i<N; i++)
    {
        A[i].home=true;
        *(A[i].pos)=N;  //N=none of the food spot: 1rst time leaving home
        *(A[i].Rew)=0;
        zero_tab_double(A[i].reward);
        //init_rew(A[i].reward_one_strategy, payoff_fn);
        init_rew(A[i].reward_diff_strategy, payoff_fn, i);
        zero_tab_int(A[i].nb_visit, N+1);
        zero_tab_int(A[i].nb_try, N+1);
        
        F[i].free=true;
        F[i].nb_visit=0;
        F[i].nb_try=0;
        F[i].A_name=N; //no agent have been there already
        F[i].P_nb_elements=0;
        F[i].arrival_time=0.;
        F[i].leaving_time=0.;
        F[i].T_occ=0.;
    }
}
/******************************************************************************/
/*************************** PERSISTENCE **************************************/
/******************************************************************************/
void persistence_up(food *F, int m, int agent_name, int N)
{
    //first time m is visited
    if (F[m].A_name == N)
    { F[m].A_name=agent_name;   F[m].P_nb_elements = 1; }      // F[m].P.sum=1;   F[m].P.nb_elements=1;}
    //a different agent visits m
    else if (F[m].A_name != agent_name)
    {F[m].A_name=agent_name;   F[m].P_nb_elements+= 1; }
    //m is visited again by the same agent  //else if (F[m].A_name == agent_name)  ---> nothing to do
}
/******************************************************************************/
/****************************  SEARCH EMPTY FOOD  *****************************/
/******************************************************************************/
//// RANDOM CHOICE OF THE FOOD
int rand_food(agent Ai, food *F, system_ *total_nb, double *cost)
{
    int m=rand_int(N); //choose a random ressource
    while(!(F[m].free))
    {
        *cost+=c;          //paye a cost to change place
        Ai.reward[m]+=c;   //update the table of payoff of Ai
        Ai.nb_try[m]+=1;        Ai.nb_try[N]+=1;
        F[m].nb_try+=1;    //update the nb of try
        (*total_nb).tries++;
        
        //implement rewards:
        //Ai.reward_one_strategy[m] = F[m].u + c * (F[m].nb_try + q_b(Ai.name, m))/(F[m].nb_visit + q_a(Ai.name, m));
        Ai.reward_diff_strategy[m]+=c/(Ai.nb_visit[m] + q_a(Ai.name, m));
        
        m=rand_int(N);  //choose a new place
    }
    return m;
}

//// CHOICE OF THE FOOD WITH A CERTAIN PROBABILITY
double* rewait(double *reward_tab, double *total_weight)
{
    *total_weight=0;
    double *backup_tab=(double *)malloc(N*sizeof(double));
    double r_a;
    
    for(int i=0; i<N; i++)
    {
        r_a=reward_tab[i];
        if (r_a==0)       {backup_tab[i]=1.;} //{backup_tab[i]=0.;} //
        else if(r_a<0)    {backup_tab[i]=0;}  //1./(1-a);} //{backup_tab[i]=0.;}
        else            {backup_tab[i]=r_a;}
        *total_weight+=backup_tab[i];
    }
    return backup_tab;
}

int best(double *w, double total_weight, food* F)
{
    double xi=uni(total_weight);
    int m=0;
    double sum=w[0]; //weight for m=0
    
    while(sum<xi)
    { m++;   sum+=w[m]; }
    if(m>=N) {cout << "error reweight function" << endl; return N;}
    else return m;
}

int best_prob_food(agent Ai, food *F, system_ *total_nb, double *cost)
{
    double total_weight=0;
    //double *weight= rewait(Ai.reward, &total_weight);
    double *weight= rewait(Ai.reward_diff_strategy, &total_weight);
    
    int m;
    if (total_weight==0)
    {   m = rand_food(Ai, F, total_nb, cost);    }
    else
    {
        m=best(weight, total_weight, F);
        while(!(F[m].free))
        {
            *cost+=c;          //paye a cost to change place
            Ai.reward[m]+=c;   //update the table of payoff of Ai
            Ai.nb_try[m]+=1;        Ai.nb_try[N]+=1;
            F[m].nb_try+=1;    //update the nb of try
            (*total_nb).tries+=1;
            
            //implement rewards:
            //Ai.reward_one_strategy[m] = F[m].u + c * (F[m].nb_try + q_b(Ai.name, m))/(F[m].nb_visit + q_a(Ai.name, m));
            Ai.reward_diff_strategy[m]+=c/(Ai.nb_visit[m] + q_a(Ai.name, m));
            
            //choose a new place
            total_weight-=weight[m]; //change the total weight
            weight[m]=0;             //and put to 0 the proba to go to m right now
            m=best(weight, total_weight, F); //choose a new place
        }
    }
    
    free(weight);
    return m;
}

//// CHOICE OF THE FOOD BY ORDERING, RANGING FROM THE BEST TO THE WORST
bool best_first (const pair<int, double>& x, const pair<int, double>& y)
{    return (x.second>y.second);    }

int best_sort_food(agent Ai, food *F, system_ *total_nb, double *cost)
{
    list< pair<int, double> > li_food;
    
    for(int m=0; m<N; m++)
    {        li_food.push_back(make_pair(m, Ai.reward_diff_strategy[m]));    }
    
    li_food.sort(best_first);
    
    list< pair<int, double> >::const_iterator it = li_food.begin();
    int m=(*it).first;
    
    while(!(F[m].free))
    {
        //COST
        *cost+=c;          //paye a cost to change place
        Ai.reward[m]+=c;   //update the table of payoff of Ai
        Ai.nb_try[m]+=1;        Ai.nb_try[N]+=1;
        F[m].nb_try+=1;    //update the nb of try
        (*total_nb).tries+=1;
        
        //implement rewards:
        //Ai.reward_one_strategy[m] = F[m].u + c * (F[m].nb_try + q_b(Ai.name, m))/(F[m].nb_visit + q_a(Ai.name, m));
        Ai.reward_diff_strategy[m]+=c/(Ai.nb_visit[m] + q_a(Ai.name, m));
        
        //go to the next place
        ++it;   m=(*it).first;
        if (it==li_food.end())   {cout << "error, all foods occupied" << endl; }
    }
    
    li_food.clear();
    return m;
}

//// CHOICE OF THE FOOD BY TAKING THE BEST ONE
double search_best_food(double *rewardd)
{
    int m=0;
    double best_reward=rewardd[m];
    
    for (int i=1; i<N; i++)
    {        if(rewardd[i]>best_reward) {best_reward=rewardd[i]; m=i;}    }
    return m;
}

int best_max_food(agent Ai, food *F, system_ *total_nb, double *cost)
{
    double *backup_reward=(double *)malloc(N*sizeof(double));
    
    for(int m=0; m<N; m++)
    {        backup_reward[m]=Ai.reward_diff_strategy[m];  }
    
    int m = search_best_food(backup_reward);
    
    int i_try=1;    //nb of tries
    
    //cout << "try=";
    
    while(!(F[m].free))
    {
        //COST
        *cost+=c;          //paye a cost to change place
        Ai.reward[m]+=c;   //update the table of payoff of Ai
        Ai.nb_try[m]+=1;        Ai.nb_try[N]+=1;
        F[m].nb_try+=1;    //update the nb of try
        (*total_nb).tries+=1;
        
        //implement rewards:
        //Ai.reward_one_strategy[m] = F[m].u + c * (F[m].nb_try + q_b(Ai.name, m))/(F[m].nb_visit + q_a(Ai.name, m));
        Ai.reward_diff_strategy[m]+=c/(Ai.nb_visit[m] + q_a(Ai.name, m));
        
        //test
        //if (m==0) { cout << F[m].nb_try << ", "; }
        
        //go to the next place
        backup_reward[m]=-100;
        m = search_best_food(backup_reward);
        
        i_try++;
        if (backup_reward[m]==-100)   {cout << "error, all foods occupied:" << i_try  << " r=" << backup_reward[m] << endl; break;}
    }
    
    free(backup_reward);
    //cout << "m = " << m << endl;
    return m;
}
/******************************************************************************/
/*********************  ACTIONS -- OUT & HOME  ********************************/
/******************************************************************************/
void go_out(int *H, vector<agent>& li_home, vector<agent>& li_out, food *F, int search_strategy(agent, food*, system_*, double*), system_ *total_nb, double t, list<int>& li_order_out)
{
    //take out 1 random agent among the ones that are at home (H)
    int i=rand_int(*H); int m;
    agent Ai = remove_agent(i, li_home);
    Ai.home=false;
    
    //order exit of the agents
    if(*(Ai.pos)==N){li_order_out.push_back(Ai.name);}
    
    //find an empty place
    double cost=0.;       //total cost //current reward of the agent
    if(uni()<=epsi)
    { m = rand_food(Ai, F, total_nb, &cost); }
    else
    { m = search_strategy(Ai, F, total_nb, &cost); }   //{ m = best_food(Ai, F, &cost); }
    
    //m is occupied
    F[m].free=false;
    F[m].nb_visit++;   (*total_nb).visits++;    F[m].nb_try++;      (*total_nb).tries++;
    F[m].arrival_time=t;
    
    //Ai goes on F[m]
    *(Ai.pos)=m;
    Ai.reward[m]+=F[m].u;       //reward on the site m
    *(Ai.Rew)+= (F[m].u + cost);   //total reward
    Ai.nb_visit[m]+=1;         Ai.nb_visit[N]+=1;
    Ai.nb_try[m]+=1;        Ai.nb_try[N]+=1;
    
    //test
    //if (m==0) { cout << F[m].nb_visit << ", visit=" << F[m].nb_visit << ", "; }
    //cout << (*total_nb).tries << " , visit=" << (*total_nb).visits << endl;
    
    //implement rewards:
    //Ai.reward_one_strategy[m] = F[m].u + c * (F[m].nb_try + q_b(Ai.name, m))/(F[m].nb_visit + q_a(Ai.name, m));
    Ai.reward_diff_strategy[m] = F[m].u + c * (Ai.nb_try[m] + q_b(Ai.name, m))/(Ai.nb_visit[m] + q_a(Ai.name, m));
    
    //implement persistence
    persistence_up(F, m, Ai.name, N);
    
    //actualise the list of agent out and the number H
    li_out.push_back(Ai);
    *H = *H-1;
}

void go_home(int *H, vector<agent>& li_home, vector<agent>& li_out, food *F, double t)
{
    //take home 1 random agent among the ones that are outside (N-H)
    int i=rand_int(N- *H);
    agent Ai = remove_agent(i, li_out);
    Ai.home=true;
    
    //free the previous position of the agent
    int m=*(Ai.pos);
    F[m].free=true;
    F[m].leaving_time=t;
    F[m].T_occ+= F[m].leaving_time - F[m].arrival_time;
    
    //actualise the list and number H of agent at home
    li_home.push_back(Ai);
    *H = *H+1;
}
/******************************************************************************/
/**********************************  TEST  ************************************/
/******************************************************************************/
/*void test_sort_list(food* F)
{
    list< pair<int, double> > li;
    
    for(int m=0; m<N; m++)
    {
        li.push_back(make_pair(m, F[m].u));
    }
    
    list< pair<int, double> >::const_iterator it;
    cout << "list: " << endl;
    for(it=li.begin(); it!=li.end(); ++it)
    {
        cout << (*it).first << ' ' << (*it).second << endl;
    }
    cout << endl << endl;
    
    li.sort(); cout << "list sorted: " << endl;
    for(it=li.begin(); it!=li.end(); ++it)
    {
        cout << (*it).first << ' ' << (*it).second << endl;
    }
    cout << endl;
    
    li.sort(best_first); cout << "list sorted: " << endl;
    for(it=li.begin(); it!=li.end(); ++it)
    {
        cout << (*it).first << ' ' << (*it).second << endl;
    }
    cout << endl;
    
    li.clear(); 
    //cout << "list: ";
      //           for(it=li.begin(); it!=li.end(); ++it)
       //          {
       //          cout << *it << ' ';
       //          }
       //          cout << endl;
}
*/