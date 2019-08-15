# README

## Requirements:

Uses the C++11 version of C++.

## Usage

### Set the parameters in `data.h`:

Specify the choice of the parameters of the system:
 - `int N` = number of agents = number of food spots;
 - `double eta` = rate at which agents leave their home;
 - `double c` = cost for changing food spot.

Specify the choice of the simulation parameters: 
 - `double epsi` = time resolution;
 - `int64_t max_Nit` = total number of actions taken by all the agents (an action = going out or going home);
 -  `int N_average` = number of systems simulated, over which measured quantities are averaged.

## Compile:

`g++ -std=c++11 -O2 library.cpp main.cpp print.cpp graph.cpp`