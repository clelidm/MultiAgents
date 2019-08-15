//g++ -std=c++11 library.cpp main.cpp print.cpp
// or g++ -std=c++0x library.cpp main.cpp print.cpp
//time ./a.out
/////-lm -lgsl -lgslcblas
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <string>
#include <list>
#include "data.h"

using namespace std;

/******************************************************************************/
/************************** GRAPH NAME ****************************************/
/******************************************************************************/
string file_graph_name(int Nit)
{
    return "Nit" + to_string(Nit) + "_" + graph_name;
    //return dir + to_string(N) + "_e"+to_string(eta) + "/Nit" + to_string(Nit) + "_" + graph_name + ".dat";
}
/******************************************************************************/
/*******************************  CREATE GRAPH ********************************/
/******************************************************************************/
void create_gephi_file(int Nit)
{
  fstream file(file_graph_name(Nit), ios::out);
  file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
  file << "<gexf version=\"1.1\">" << endl;
  file << "<meta lastmodifieddate=\"2010-03-03+23:44\"> </meta>" << endl;
  file << "<graph defaultedgetype=\"undirected\" idtype=\"string\" type=\"static\">" << endl;
  file.close();
}

void print_nodes(int Nit)
{	
  fstream file(file_graph_name(Nit), ios::out|ios::app);
  file << "<nodes count=\"" << N << "\">" << endl;
  for(int i=0; i<N; i++)
  {      file << "<node id=\"" << i << "\" label=\"" << i << "\"/>" << endl;	}
  file << "</nodes>" << endl;
  file.close(); 
}

void create_edges(list< pair<int, double> >& li_A, int nb_eff, list<edge>& li_E) //int nb_edges
{
  edge E;
  list< pair<int, double> >::iterator it1=li_A.begin();
  list< pair<int, double> >::iterator it2=it1;
  
  //list< pair<int, double> >::const_iterator it = li_food.begin();
  
  int n_id=li_E.size()+1;
  
  for(int i1=0; i1<nb_eff; i1++)
  {
    it2=it1;
    //it1=li_A.begin()+i1;
    for(int i2=i1+1; i2<nb_eff; i2++)
    {
      it2++;
      //it2=li_A.begin()+i2;  //it2=it1+1;
      //edge:
      E.id=n_id++;
      E.source=(*it1).first;
      E.target=(*it2).first;
      E.weight=1.;
      li_E.push_back(E);
    }
    it1++;
  }
}

void print_edges(list<edge>& li_E, int Nit) //int nb_edges
{
  fstream file(file_graph_name(Nit), ios::out|ios::app);
    
  list<edge>::const_iterator it;
  int i=0; 
  file << "<edges count=\"" << li_E.size() << "\">" << endl; //li_E.size();
  for(it=li_E.begin(); it!=li_E.end(); ++it)
  {
    file << "<edge id=\"" << i <<"\" source=\"" << (*it).source << "\" target=\"" << (*it).target << "\" weight=\"" << (*it).weight << "\"/>" << endl;
    i++;
  }
  file << "</edges>" << endl;
  file.close();
}

void close_graph(int Nit)
{
  fstream file(file_graph_name(Nit), ios::out|ios::app);
  file << "</graph>" << endl << "</gexf>" << endl;
  file.close();
}	

void create_graph(list<edge> &li_E, int Nit)
{
  create_gephi_file(Nit);
  print_nodes(Nit);
  print_edges(li_E, Nit);
  close_graph(Nit);
}







