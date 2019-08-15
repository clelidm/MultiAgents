#include <cstdlib>	//permet d'utiliser alea (et pas que)
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include "data.h"
using namespace std;

/******************************************************************************/
/*******************************  FILE NAMES **********************************/
/******************************************************************************/
string file_name(int Nit, string name)
{
    return dir + to_string(N) + "_e"+to_string(eta) + "/" + name + "_Nit" + to_string(Nit) + ".dat";
}
string file_name2(int Nit, string name)
{
    return dir + to_string(N) + "_e"+to_string(eta) + "/Nit" + to_string(Nit) + "_" + name + ".dat";
}
string file_name_rk(int Nit, string name, int rank)
{
    return dir + to_string(N) + "_e"+to_string(eta) + "/record/Nit" + to_string(Nit) + "_" + name + to_string(rank) +".dat";
}
/******************************************************************************/
/******************************** STORE DATA  *********************************/
/******************************************************************************/
void make_distrib(double x, int *tab)
{
    if (x>=histo_x_min && x<=histo_x_max)
    {
        int k = (int)((x-histo_x_min)/histo_Dx);
        tab[k]++;
    }
}

void effective_nb_Ai(agent *A, int i, double *eff_visit_Ai, double *eff_try_Ai)
{
    double p_visit=0, p_try=0;
    double I_tried_Ai=0;    //sum(p_try**2)
    double I_visited_Ai=0;  //sum(p_visit**2)
    
    for(int m=0; m<N; m++)
    {
        p_visit=((double)A[i].nb_visit[m])/A[i].nb_visit[N];
        p_try=((double)A[i].nb_try[m])/A[i].nb_try[N];
        
        I_tried_Ai+=(p_try*p_try);           I_visited_Ai+=(p_visit*p_visit);
    }
    *eff_try_Ai=1./I_tried_Ai;           *eff_visit_Ai=1./I_visited_Ai;
}

void effective_nb_at_m(agent *A, food F, double *eff_visit_Fm, double *eff_try_Fm)
{
    int m=F.name;
    double p_visit=0, p_try=0;
    double I_tried_m=0;    //sum(p_try**2)
    double I_visited_m=0;  //sum(p_visit**2)
    
    if(F.nb_try==0){*eff_try_Fm=0; *eff_visit_Fm=0;} //F has never been visited
    else
    {
        for(int i=0; i<N; i++)  //all agents visiting the food place F
        {
            p_try=((double)A[i].nb_try[m])/F.nb_try;
            p_visit=((double)A[i].nb_visit[m])/F.nb_visit;
        
            I_tried_m+=(p_try*p_try);           I_visited_m+=(p_visit*p_visit);
        }
    *eff_try_Fm=1./I_tried_m;        *eff_visit_Fm=1./I_visited_m;
    }
}

void store_data(results *R, histo* H, food *F, agent *A, system_ total_nb, double t)
{
    int nb_tot_try=0, nb_tot_visit=0;
    
    // Pi_m = #_try[m] / total_nb_try
    // Q_m = #_visit[m] / total_nb_visit    --> gives the most occupied sites
    // e_m = #_visit[m] / #_try[m] or 0 if #_try[m]=0  == pb m is empty   (visit = success try)
    // P_m = #_visit[m] / #_diff_successive_agent visiting m = persistence of the site m
    // nb_A = #_diff_agent = total nb of different agents visiting m
    
    (*R).N_av++; //one more element for computing the average

    for (int m = 0; m < N; m++)
    {
        //nb_visit and nb_try
        nb_tot_try+=F[m].nb_try;        nb_tot_visit+=F[m].nb_visit;
        (*R).nb_try[m]+=F[m].nb_try;    (*R).nb_visit[m]+=F[m].nb_visit;
        
        (*R).Pi[m] += ((float)F[m].nb_try)/total_nb.tries;    // Pi_m
        (*R).Q[m] +=((float)F[m].nb_visit)/total_nb.visits;  //Q_m
        //e_m
        if (F[m].nb_try == 0)   { (*R).e[m]+=1; }    //the site m is always empty
        else { (*R).e[m]+=((float)F[m].nb_visit)/F[m].nb_try; }
        //P_m
        if (F[m].P_nb_elements==0)
        {
            if (F[m].nb_visit!=0) { cout << "error in persistence" << endl;}
            (*R).P[m]+=0.;
        }
        else
        {   (*R).P[m]+=((float)F[m].nb_visit)/F[m].P_nb_elements;  }
        //e_real[m]
        (*R).e_real[m]+=1-(F[m].T_occ/t);
    }
    cout << "For store_data function:" << endl << "total nb of try = " << nb_tot_try << " ; total nb of visit = " << nb_tot_visit << endl << endl;
    
    double eff_visit_Ai=0;
    double eff_try_Ai=0;
    double eff_visit_m=0;
    double eff_try_m=0;
    //store for the histogrammes
    for(int i=0;i<N;i++)
        {
        // for each agent
        effective_nb_Ai(A, i, &eff_visit_Ai, &eff_try_Ai);
        make_distrib(eff_visit_Ai, (*H).histo_visit);
        make_distrib(eff_try_Ai, (*H).histo_try);
        // for each food
        effective_nb_at_m(A, F[i], &eff_visit_m, &eff_try_m);
        (*R).eff_nb_try[i]+=eff_try_m;    (*R).eff_nb_visit[i]+=eff_visit_m;
        make_distrib(eff_visit_m, (*H).histo_visit_m);
        make_distrib(eff_try_m, (*H).histo_try_m);
            
        (*H).N_samples++;
        //stat distrib nb visit per agent i
        (*H).mean_try+=eff_try_Ai;
        (*H).mean_visit+=eff_visit_Ai;
        (*H).m2_visit+=eff_visit_Ai*eff_visit_Ai;
        (*H).m2_try+=eff_try_Ai*eff_try_Ai;
        //stat distrib nb visit per food m
        (*H).mean_try_m+=eff_try_m;
        (*H).mean_visit_m+=eff_visit_m;
        (*H).m2_visit_m+=eff_visit_m*eff_visit_m;
        (*H).m2_try_m+=eff_try_m*eff_try_m;
        }
}
/******************************************************************************/
/******************************** PRINT DATA  *********************************/
/******************************************************************************/
void print_R_and_H(results R, histo H, food* F, int Nit)
{
    //print R
    fstream file(file_name2(Nit, "Fav_m"), ios::out);
    file << "# m:1  \t  u_m:2  \t  nb_try:3  \t  nb_visit:4  \t  Pi_m:5  \t  q_m:6  \t  P_m:7  \t  Q_m:8 \t e_m:9 \t eff_nb_try:10 \t eff_nb_visit:11 \n";
    for(int m = 0; m < N; m++)
    {
        file << m << "\t" << F[m].u << "\t" << ((double)R.nb_try[m])/R.N_av << "\t" << ((double)R.nb_visit[m])/R.N_av << "\t" << R.Pi[m]/R.N_av  << "\t"<< R.e[m]/R.N_av << "\t" << R.P[m]/R.N_av << "\t" << R.Q[m]/R.N_av << "\t" << R.e_real[m]/R.N_av << "\t" << R.eff_nb_try[m]/R.N_av << "\t" << R.eff_nb_visit[m]/R.N_av << "\n";
    }
    file.close();
    
    //print H
    fstream file2(file_name2(Nit, "distrib_eff_nb"), ios::out);
    file2 << "# eff_nb:1  \t  pdf_try_perA:2  \t  pdf_visit_perA:3 \t  pdf_try_perF:4  \t  pdf_visit_perF:5\n";
    double norm=H.N_samples*histo_Dx;
    cout << "nb_samples=" << H.N_samples << " " << N*N_average << endl;
    for(int z=0; z<histo_N_steps; z++)
    {   file2 << histo_x_min+z*histo_Dx << "\t\t" << (double)H.histo_try[z]/norm << "\t\t" << (double)H.histo_visit[z]/norm << "\t\t" << (double)H.histo_try_m[z]/norm << "\t\t" << (double)H.histo_visit_m[z]/norm <<endl;  }
    file2.close();
}

void print_stored_data(results R, histo H, food* F, int Nit)
{
    print_R_and_H(R, H, F, Nit);
    
    //statistique distrib nb visit per Agent Ai
    fstream file3(file_name2(Nit, "eff_nb_stat_per_Ai"), ios::out);
    file3 << "# Nit:1 \t <nb_try>:2  \t  <m2_try>:3  \t  var_try:4 \t <nb_visit>:5  \t  <m2_visit>:6  \t  var_visit:7 \n";
    
    double m1 =H.mean_try/H.N_samples;      double m2=H.m2_try/H.N_samples;         double sig= m2-m1*m1;
    file3 << Nit << "\t" << m1 << "\t" << m2 << "\t" << sig << "\t";
    
    m1=H.mean_visit/H.N_samples;            m2=H.m2_visit/H.N_samples;              sig= m2-m1*m1;
    file3 << m1 << "\t" << m2 << "\t" << sig << endl;
    file3.close();
    
    //statistique distrib nb visit per Food F[m]
    fstream file4(file_name2(Nit, "eff_nb_stat_per_Fm"), ios::out);
    file4 << "# Nit:1 \t <nb_try>:2  \t  <m2_try>:3  \t  var_try:4 \t <nb_visit>:5  \t  <m2_visit>:6  \t  var_visit:7 \n";
    
    m1 =H.mean_try_m/H.N_samples;           m2=H.m2_try_m/H.N_samples;              sig= m2-m1*m1;
    file4 << Nit << "\t" << m1 << "\t" << m2 << "\t" << sig << "\t";
    
    m1=H.mean_visit_m/H.N_samples;          m2=H.m2_visit_m/H.N_samples;            sig= m2-m1*m1;
    file4 << m1 << "\t" << m2 << "\t" << sig << endl;
    file4.close();
}

void print_stored_data_tab(results *R, histo *H, food *F, int Nit)
{
    for (int i_Nit=min_i_Nit; i_Nit <= Nit; i_Nit++)
    {
        print_R_and_H(R[i_Nit-1], H[i_Nit-1], F, i_Nit);
        /*
        fstream file(file_name2(i_Nit, "Fav_m"), ios::out);
    
        file << "# m:1  \t  u_m:2  \t  nb_try:3  \t  nb_visit:4  \t  Pi_m:5  \t  e_m:6  \t  P_m:7  \t  Q_m:8 \n";
    
        for(int m = 0; m < N; m++)
        {
            file << m << "\t" << F[m].u << "\t" << ((double)R[i_Nit-1].nb_try[m])/R[i_Nit-1].N_av << "\t" << ((double)R[i_Nit-1].nb_visit[m])/R[i_Nit-1].N_av << "\t" << R[i_Nit-1].Pi[m]/R[i_Nit-1].N_av  << "\t"<< R[i_Nit-1].e[m]/R[i_Nit-1].N_av << "\t" << R[i_Nit-1].P[m]/R[i_Nit-1].N_av << "\t" << R[i_Nit-1].Q[m]/R[i_Nit-1].N_av << "\n";
        }
        file.close();
        
        // histo
        fstream file2(file_name2(i_Nit, "distrib_eff_nb"), ios::out);
        
        file2 << "# eff_nb:1  \t  pdf_try:2  \t  pdf_visit:3 \n";
        
        double norm=H[i_Nit-1].N_samples*histo_Dx;
        cout << "nb_samples=" << H[i_Nit-1].N_samples << " " << N*N_average << endl;
        for(int z=0; z<histo_N_steps; z++)
        {   file2 << histo_x_min+z*histo_Dx << "\t\t" << (double)H[i_Nit-1].histo_try[z]/norm << "\t\t" << (double)H[i_Nit-1].histo_visit[z]/norm <<endl;  }
        
        file2.close();
         */
    }
    
    fstream file3(file_name2(Nit, "eff_nb_stat"), ios::out);
    file3 << "# Nit:1 \t <nb_try>:2  \t  <m2_try>:3  \t  var_try:4 \t <nb_visit>:5  \t  <m2_visit>:6  \t  var_visit:7 \n";
    double m1, m2, sig;
    for (int i_Nit=min_i_Nit; i_Nit <= Nit; i_Nit++)
    {
    m1 = H[i_Nit-1].mean_try/H[i_Nit-1].N_samples;
    m2 = H[i_Nit-1].m2_try/H[i_Nit-1].N_samples;
    sig = m2-m1*m1;
    file3 << i_Nit << "\t" << m1 << "\t" << m2 << "\t" << sig << "\t";
    m1 =H[i_Nit-1].mean_visit/H[i_Nit-1].N_samples;
    m2=H[i_Nit-1].m2_visit/H[i_Nit-1].N_samples;
    sig= m2-m1*m1;
    file3 << m1 << "\t" << m2 << "\t" << sig << endl;
    }
    file3.close();
}
/******************************************************************************/
/******************************  PRINTS FOOD **********************************/
/******************************************************************************/
void print_data_food(food *F, int Nit, int total_nb_try, int total_nb_visit)
{
    int nb_tot_try=0, nb_tot_visit=0;
    double Pim, Qm, em=0, Pm;
    fstream file(file_name2(Nit, "F_m"), ios::out);
    
    file << "# m:1  \t  u_m:2  \t  nb_try:3  \t  nb_visit:4  \t  Pi_m:5  \t  e_m:6  \t  P_m:7  \t  Q_m:8 \n";
    // Pi_m = #_try[m] / total_nb_try
    // Q_m = #_visit[m] / total_nb_visit    --> gives the most occupied sites
    // e_m = #_visit[m] / #_try[m] or 0 if #_try[m]=0  == pb m is empty   (visit = success try)
    // P_m = #_visit[m] / #_diff_successive_agent visiting m = persistence of the site m
    // nb_A = #_diff_agent = total nb of different agents visiting m
    
    for (int m = 0; m < N; m++)
    {
        nb_tot_try+=F[m].nb_try;        nb_tot_visit+=F[m].nb_visit;
        
        Pim = ((float)F[m].nb_try)/total_nb_try;    // Pi_m
        Qm =((float)F[m].nb_visit)/total_nb_visit;  //Q_m
        //e_m
        if (F[m].nb_try == 0)   { em=1; }    //the site m is always empty
        else { em=((float)F[m].nb_visit)/F[m].nb_try; }
        //P_m
        if (F[m].P_nb_elements==0)
        {
            if (F[m].nb_visit!=0) { cout << "error in persistence" << endl;}
            Pm=0.;
        }
        else
        {   Pm=((float)F[m].nb_visit)/F[m].P_nb_elements;  }
        
        //print
        file << m << "\t" << F[m].u << "\t" << F[m].nb_try << "\t" << F[m].nb_visit << "\t" << Pim  << "\t"<< em << "\t" << Pm << "\t" << Qm << "\n";
    }
    
    cout << "For print_Food function with Nit=10**" << Nit << ":" << endl << "total nb of try = " << nb_tot_try << " ; total nb of visit = " << nb_tot_visit << endl << endl;
    
    file.close();
}

void print_nb_visit(agent *A, int Nit)
{
    fstream file1(file_name2(Nit, "histo_visit"), ios::out);
    fstream file2(file_name2(Nit, "eff_nb_site_visited"), ios::out);
    
    //double total_nb_try_Ai=0;
    double p_visit=0, p_try=0;
    double I_tried_Ai=0;    //sum(p_try**2)
    double I_visited_Ai=0;  //sum(p_visit**2)
    
    file2 << "# Ai \t eff_nb_sites_tried \t eff_nb_sites_visited" << endl;
    for(int i=0; i<N; i++)
    {
        I_tried_Ai=0;    //sum(p_try**2)
        I_visited_Ai=0;  //sum(p_visit**2)
        file1 << "# Ai \t site m \t nb_visit Pba_visit_m \t nb_try Pba_try_m" << endl;
        for(int m=0; m<N; m++)
        {
            p_visit=((double)A[i].nb_visit[m])/A[i].nb_visit[N];        p_try=((double)A[i].nb_try[m])/A[i].nb_try[N];
            
            file1 << i << "\t" << m << "\t" << A[i].nb_visit[m] << "\t" << p_visit << "\t" << A[i].nb_try[m] << "\t" << p_try << endl;
            
            I_tried_Ai+=(p_try*p_try);                  I_visited_Ai+=(p_visit*p_visit);
        }
        file1 << endl;
        file2 << i << "\t" << 1./I_tried_Ai << "\t" << 1./I_visited_Ai << endl;
    }
    
    file1.close();
    file2.close();
}

void print_kn_i_m(agent *A, int Nit, int N_file)
{
    fstream file1(file_name_rk(Nit, "k_i_m", N_file), ios::out);
    file1 << "# k_i_m corresponds to the nb of actual visits of agent i to site m (the site was free -> successful try)" << endl;
    file1 << "# lignes : Ai ; columns : m   --> both from 0 to N" << endl;

    fstream file2(file_name_rk(Nit, "n_i_m", N_file), ios::out);
    file2 << "# n_i_m corresponds to the nb of tries of agent i to site m (successfull and unsuccessfull tries)" << endl;
    file2 << "# lignes : Ai ; columns : m   --> both from 0 to N" << endl;
    
    for(int i=0; i<N; i++) //agents
    {
        for(int m=0; m<N; m++) // food spots
        {
            file1 << A[i].nb_visit[m] << "\t";
	    file2 << A[i].nb_try[m] << "\t";
        }
        file1 << endl;
        file2 << endl;
    }
    file1.close();
    file2.close();
}

void print_e_m(food *F, int Nit, double t)
{   //t is the total simulation time
    fstream file(file_name2(Nit, "e_m"), ios::out);
    double *e=(double*)malloc(N*sizeof(double));
    
    file << "#m \t e_m"<< endl;
    
    for(int m=0; m<N; m++){
        e[m]=1-(F[m].T_occ/t);
        file << m << " \t " << e[m] << endl;
    }
    file.close();
}
/******************************************************************************/
/*********************************  PRINTS 1 **********************************/
/******************************************************************************/
/*
void print_Persistence(food *F, int Nit, int N, double eta)
{
    fstream file(file_name(Nit, "P_m"), ios::out);
    for(int m=0; m<N; m++)
    {   if (F[m].P_nb_elements==0)
    {
        if (F[m].nb_visit!=0) { cout << "error in persistence" << endl;}
        else { file << m << "\t" << 0 << endl; }
    }
    else
    {   file << m << "\t" << ((float)F[m].nb_visit)/F[m].P_nb_elements << endl;  }
    }
    file.close();
}

void Pi_m_fake(food *F, int Nit)   // print_nb_visit
{
    int nb_tot=0;
    fstream file(file_name(Nit, "Pi_m_bis"), ios::out);
    for (int i = 0; i < N; i++)
    {
        nb_tot+=F[i].nb_visit;
        file << i << "\t" << ((float)F[i].nb_visit)/total_nb_visit << "\t" << F[i].u << "\n";
    }
    cout << "total number of visit = " << nb_tot << " " << total_nb_visit << endl << endl;
    file.close();
}

void Pi_m(food *F, int Nit)   // print_nb_visit
{
    int nb_tot=0;
    fstream file(file_name(Nit, "Pi_m"), ios::out);
    for (int i = 0; i < N; i++)
    {
        nb_tot+=F[i].nb_try;
        file << i << "\t" << ((float)F[i].nb_try)/total_nb_try << "\t" << F[i].u << "\n";
    }
    cout << "for Pi_m: total nb of try = " << nb_tot << " " << total_nb_try << endl << endl;
    file.close();
}

void e_m(food *F, int Nit)
{
    int nb_tot=0; double e=0;
    fstream file(file_name(Nit, "e_m"), ios::out);
    file << "# m \t e_m \t U_m \t nb_visit \n";
    for (int i = 0; i < N; i++)
    {
        nb_tot+=F[i].nb_try;
        if (F[i].nb_try == 0)   { e=1; } //the site m is always empty
        else { e=((float)F[i].nb_visit)/F[i].nb_try; }
        file << i << "\t" << e << "\t" << F[i].u  << "\t"<< F[i].nb_visit << "\n";
    }
    cout << "for e_m: total nb of try = " << nb_tot << " " << total_nb_try << endl << endl;
    file.close();
}

*/