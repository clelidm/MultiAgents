//g++ -std=c++11 library.cpp main.cpp print.cpp graph.cpp
// or g++ -std=c++0x library.cpp main.cpp print.cpp graph.cpp
//time ./a.out
/////-lm -lgsl -lgslcblas
#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <list>
#include <vector>
#include <string>
#include "library.h"
//#include "data.h"
//#include <cinttypes> //for type int64_t in C++11 -- not necessary. Necessary to call INT64_MAX

using namespace std;

/******************************************************************************/
/********************** INITIALISATION ****************************************/
/******************************************************************************/
double fn_payoff(int i, int order)
{
    //if (order%2==0) {    return 1+9*pow((i/(N-1.) - 1.), order); }
    //else            {    return 1+9*pow(1-(i/(N-1.)), order); }
    return 1 + 9*pow(1-(i/(N-1)), order);
}

double payoff(int i)
{   return 10.-9.*i/(N-1);    } //or {   return fn_payoff(i, 2);   }
/******************************************************************************/
/***************************** CREATE HISTO ***********************************/
/******************************************************************************/
double effective_nb_visits_for_Ai(agent *A, int Nit, int i)
{
    //double total_nb_try_Ai=0;
    double p_visit=0, p_try=0;
    double I_tried_Ai=0;    //sum(p_try**2)
    double I_visited_Ai=0;  //sum(p_visit**2)
    double eff_nb_tried_Ai=0;    // 1/I_tried_Ai
    double eff_nb_visited_Ai=0;  // 1/I_visited_Ai
    
    I_tried_Ai=0;    //sum(p_try**2)
    I_visited_Ai=0;  //sum(p_visit**2)
    
    for(int m=0; m<N; m++)
    {
        p_visit=((double)A[i].nb_visit[m])/A[i].nb_visit[N];
        p_try=((double)A[i].nb_try[m])/A[i].nb_try[N];
        
        I_tried_Ai+=(p_try*p_try);           I_visited_Ai+=(p_visit*p_visit);
    }
    eff_nb_tried_Ai=1./I_tried_Ai;           eff_nb_visited_Ai=1./I_visited_Ai;

    return eff_nb_visited_Ai;
}
/******************************************************************************/
void make_histo(int N_steps, double x_max, double x_min, int N_realizations, agent* A, int Nit)// int max_it
{
    double *t=(double *)malloc(N_steps*sizeof(double));
    double sum=0, norm=0, x;
    double delta_x = (x_max - x_min)*1./(N_steps);
    int k,i; //n_max=0;
    
    for(int k =0; k<N_steps; k++){t[k]=0;}
    
    //while(n_max < N_realizations){
    
    //for(int m=0; m < N_realizations; m++){
        
        //for each agent i:
        for(i=0;i<N;i++){
            x = effective_nb_visits_for_Ai(A,Nit,i);
            if (x>=x_min && x<=x_max)
            {
                k = (int)((x-x_min)/delta_x);
                t[k]++;
                sum = sum + 1;
            }
            else {
                cout << endl << "Increase x_max" <<endl;
            }
        }
    //}
    //n_max++;
    //}
    for(int m=0; m < N_steps; m++){
        t[m]=t[m]*1./(N*N_realizations*delta_x);
        norm+=t[m];
    }
    //check
    //cout << "The sum is normalized: " << norm << " " << (double) sum/N_realizations <<endl;
    
    
    fstream filee ("distrib-eta"+to_string(eta)+"-"+to_string(N_realizations)+"realiz-"+"Nit"+to_string(Nit)+".txt", ios::out);
    if(filee.is_open())
    {
        filee << "i \t\tcounts[i]"<<endl<<endl;
        for(int z=0; z<=N_steps; z++)
        {filee << x_min+z*delta_x << "\t\t" << t[z] <<endl;}
    }
    else
    {
        cout << "File could not be opened." << endl;
    }
    filee.close();
}

/******************************************************************************/
/***************************** PRINT_PAYOFF ***********************************/
/******************************************************************************/
void print_payoff(agent *A, int Nit, list<int>& li_order_out)
{
    if (li_order_out.size()!=N) {cout << "Some agents has not been out yet" << endl;}
    list<int>::iterator it;
    int i, compt=0;
    double rr=0;
    
    fstream file(file_name2(Nit, "payoff"), ios::out);
    file << "#order\t\t Ai.name \t\t mean_payoff of Ai Ui" << endl;
    for(it=li_order_out.begin(); it!=li_order_out.end(); ++it)
        {
            i =*it; rr=0;
            for (int j=0;j<N;j++)
                { rr+=A[i].reward[j]; }
            if (A[i].nb_visit[N]==0) {file << compt << " " << i << "\t\t" << 0 << "\t\t" << 0 << endl; }
            else
                { file << compt << " " << i << "\t\t" << *(A[i].Rew)/(A[i].nb_visit[N]) << "\t\t" << rr/(A[i].nb_visit[N]) << endl; }
            compt++;
        }
    file.close();
}

void print_payoff2(agent *A, int Nit, list<int>& li_order_out)
{
    double rr=0;
    
    fstream file(file_name2(Nit, "payoff"), ios::out);
    file << "#Ai \t\t mean_payoff of Ai Ui" << endl;
    for (int i=0; i<N; i++)
        {
        rr=0;
        for (int j=0;j<N;j++)
            { rr+=A[i].reward[j]; }
        if (A[i].nb_visit[N]==0) {file << i << "\t\t" << 0 << "\t\t" << 0 << endl; }
        else
        { file << i << "\t\t" << *(A[i].Rew)/(A[i].nb_visit[N]) << "\t\t" << rr/(A[i].nb_visit[N]) << endl; }
    }
    file.close();
}
/******************************************************************************/
/***************************** RUN ********************************************/
/******************************************************************************/
int run(agent* A, food *F, system_ *total_nb, int search_strategy(agent, food*, system_*, double*), bool print, int64_t N_act=1000, int Nit=7) //N_act = nb of actions
{
    // variables utiles
    int64_t i, i_Nit=10; 
    int pow_Nit=1;
    //agent AA;
    
    //initilisation
    int H=N; // nb of agent at home
    double t=0.;
    double Dt=0, Dt_h=0, Dt_o=0.;
    vector<agent> li_home;
    vector<agent> li_out;
    
    list<int> li_order_out; //order of the agents going out
    //int compt=0;

    li_out.clear();         //nobody out
    for(i=0; i<N; i++)      //all agents at home
        {   li_home.push_back(A[i]);   }

    //************************* RUN
    //Phome //fstream file("P_home_t.dat", ios::out);
    //Phome //file << 0 << ' ' << 1. << endl;
    for(i=1; i<=N_act; i++)
    {
        //cout << i <<": ";
        if(H==N)
        {
            Dt=rand_exp(eta*H);  //going out: rate eta*H
                    //cout << "GO. home size= " << li_home.size() << " --> ";
            go_out(&H, li_home, li_out, F, search_strategy, total_nb,t+Dt, li_order_out);
                    //cout << li_home.size() << " nb_visit_sites: ";
                    //AA = li_out.back();
                    //cout << AA.reward.size() << endl;
        }
        else if (H==0)
        {
            Dt=rand_exp(N);  //going home: rate mu*(N-H) with mu=1
            go_home(&H, li_home, li_out, F, t+Dt);
        }
        else
        {
            Dt_o=rand_exp(eta*H);  //going out: rate eta*H
            Dt_h=rand_exp(N-H);  //going home: rate mu*(N-H) with mu=1
            if (Dt_o < Dt_h)
            {
                        //cout << "GO. home size= " << li_home.size() << " --> ";
                Dt=Dt_o;  //going out
                go_out(&H, li_home, li_out, F, search_strategy, total_nb, t+Dt, li_order_out);
                        //cout << li_home.size() << " nb_visit_sites: ";
                        //AA = li_out.back();
                        //cout << AA.reward.size() << endl;
            }
            else
            {
                Dt=Dt_h;  //going out
                go_home(&H, li_home, li_out, F, t+Dt);
            }
        }
        t+=Dt;
        //Phome //file << t << ' ' << ((float)H)/N << endl;
        
        //#############  PRINT #############
        if(i==i_Nit && print)
        {   
            if(pow_Nit>=min_i_Nit)
            {
                print_data_food(F, pow_Nit, (*total_nb).tries, (*total_nb).visits);
                make_histo(histo_N_steps, histo_x_max, histo_x_min, 1, A, pow_Nit);   //(int N_steps, double x_max, double x_min, int N_realizations, agent* A, int Nit)
                print_e_m(F, pow_Nit, t);
                
                if(pow_Nit==Nit) { print_kn_i_m(A, pow_Nit, 1); }
            }
            i_Nit=i_Nit*10;
            pow_Nit++;
        }
    }
    return t;
    //#############  PRINT FIN #############
    //cout << endl << (li_home.front()).reward.size() << endl;
    //print_reward_for_agent_i(li_home.front());
}

void run_store(agent* A, food *F, results *R_av, histo *histo_av, system_ *total_nb, int search_strategy(agent, food*, system_*, double*), bool print, int64_t N_act=1000, int Nit=3) //N_act = nb of actions
{
    // variables utiles
    int64_t i, i_Nit=10;
    int pow_Nit=1;
    //agent AA;
    
    //initilisation
    int H=N; // nb of agent at home
    double t=0.;
    double Dt=0, Dt_h=0, Dt_o=0.;
    vector<agent> li_home;
    vector<agent> li_out;
    
    list<int> li_order_out;
    list<int>::iterator it;
    bool li_print=true;
    //int compt=0; //count agents going out, until they all went out ones
    
    //fstream file(file_name2(0, "position"), ios::out);
    bool record; //to record only after a certain time
    
    li_out.clear();         //nobody out
    for(i=0; i<N; i++)      //all agents at home
    {   li_home.push_back(A[i]);   }
    
    //************************* RUN
    for(i=1; i<=N_act; i++)
    {
        if(H==N)
        {
            Dt=rand_exp(eta*H);  //going out: rate eta*H
            go_out(&H, li_home, li_out, F, search_strategy, total_nb, t+Dt, li_order_out);
        }
        else if (H==0)
        {
            Dt=rand_exp(N);  //going home: rate mu*(N-H) with mu=1
            go_home(&H, li_home, li_out, F, t+Dt);
        }
        else
        {
            Dt_o=rand_exp(eta*H);  //going out: rate eta*H
            Dt_h=rand_exp(N-H);  //going home: rate mu*(N-H) with mu=1
            if (Dt_o < Dt_h)
            {
                Dt=Dt_o;  //going out
                go_out(&H, li_home, li_out, F, search_strategy, total_nb, t+Dt, li_order_out);
            }
            else
            {
                Dt=Dt_h;  //going out
                go_home(&H, li_home, li_out, F, t+Dt);
            }
        }
        t+=Dt;

        /*if (li_print) {
            cout << i << " " << li_order_out.size() << "\t";
            for(it=li_order_out.begin(); it!=li_order_out.end(); ++it)
                {cout << *it << " ";}
            cout << endl;
            if (li_order_out.size()==N) {li_print=false;}
        }*/
        //file << i << "\t";
        //for(int j=0; j<N; j++) {file << *(A[j].pos) << "\t";}
        //file << endl;
        //#############  STORE #############
        if(i==i_Nit)
        {
            //store_data(&R_av[i_Nit-1], F, total_nb);
            if(pow_Nit>=min_i_Nit)  {
                cout << "Nit10**=" << pow_Nit << " ";  store_data(R_av+(pow_Nit-1), histo_av+(pow_Nit-1), F, A, *total_nb, t);
                print_payoff(A, pow_Nit, li_order_out);
                }
            i_Nit=i_Nit*10;
            pow_Nit++;
        }
    }
    //file.close();
}
/******************************************************************************/
/***************************** RUN  AVERAGE ***********************************/
/******************************************************************************/
void run_average(int search_strategy(agent, food*, system_*, double*), int N_av=100, int64_t N_act=1000, int Nit=3)
{
    bool print=false; double t=1.;

    // INITIALISATION
    agent* A=init_agent(N, payoff);
    food* F= init_food(N, payoff);
    system_ total_nb;     total_nb.visits=0;     total_nb.tries=0;
    results R_av = init_results(N);
    histo H_av = init_histo(N);
    
    for(int i=0; i<N_av; i++)
    {
        cout << "Nb averages=" << i << " ";
        t=run(A, F, &total_nb, search_strategy, print, N_act, Nit);      //1e8 = 20 min for N=100 agents
        //print_data_food(F, 51, total_nb.tries, total_nb.visits);
        store_data(&R_av, &H_av, F, A, total_nb, t);
        print_kn_i_m(A, Nit, i);
        re_init_all(A, F, &total_nb, payoff);
    }
    print_stored_data(R_av, H_av, F, Nit);
}

void run_average2(int search_strategy(agent, food*, system_*, double*), int N_av=100, int64_t N_act=1000, int Nit=3)
{
    bool print=false;
    
    // INITIALISATION
    agent* A=init_agent(N, payoff);
    food* F= init_food(N, payoff);
    system_ total_nb;     total_nb.visits=0;     total_nb.tries=0;
    results *R_av = init_results_tab(N, Nit);  //create a tab of Nit results (for Nit=2 to Nit_max)
    histo *H_av = init_histo_tab(N, Nit);
    
    for(int i=0; i<N_av; i++)
    {
        cout << "Nb averages=" << i << endl;
        run_store(A, F, R_av, H_av, &total_nb, search_strategy, print, N_act, Nit);      //1e8 = 20 min for N=100 agents
        //print_data_food(F, 51, total_nb.tries, total_nb.visits);
        print_kn_i_m(A, Nit, i);
        re_init_all(A, F, &total_nb, payoff);
    }
    print_stored_data_tab(R_av, H_av, F, Nit);
    free(R_av);
    free(H_av);
}
/******************************************************************************/
/****************************  SEARCH EMPTY FOOD  *****************************/
/******************************************************************************/
int rand_food(agent Ai, food *F, system_ *total_nb, double *cost);
int best_prob_food(agent Ai, food *F, system_ *total_nb, double *cost);
int best_sort_food(agent Ai, food *F, system_ *total_nb, double *cost);
int best_max_food(agent Ai, food *F, system_ *total_nb, double *cost);

/******************************************************************************/
/************************** MAIN ***********************************************/
/******************************************************************************/
int main()
{
    // INITIALISATION
    initialise_rand();
    agent* A=init_agent(N, payoff);
    food* F= init_food(N, payoff);
    system_ T;     T.visits=0;     T.tries=0;
    
    // CREATE FOLDER "data_N"
    string str = "mkdir -p "+dir+to_string(N)+"_e"+to_string(eta);
    string str2 = "mkdir -p "+dir+to_string(N)+"_e"+to_string(eta)+"/record";
    const char *create_file = str.c_str();
    const char *create_file2 = str2.c_str();
    system(create_file); //windows : system("cmd /c tskill explorer");
    system(create_file2); //windows : system("cmd /c tskill explorer");
    
    // RUN
    bool print=true;
    //run(A, F, &T, best_max_food, print, max_Nit, max_i_Nit);      //1e8 = 20 min for N=100 agents
    
    // RUN AVERAGE
    //run_average(best_max_food, N_average,  max_Nit, max_i_Nit);
    run_average2(best_max_food, N_average, max_Nit, max_i_Nit);
    
    //make_histo(500,20,0,50000);
    
    // FREE MEMORY
    free(A);
    free(F);
    
  return 0;
}

//evolution of nb of try by one agent before doing one visit during time = evolution of proba of success during time for each agent
//how may tries in one spot m, from when it starts to be occupied, up to when it is empty again = evolution of proba of success during time on one spot
