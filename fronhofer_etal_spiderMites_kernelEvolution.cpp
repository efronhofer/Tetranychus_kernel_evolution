// code for
// Fronhofer, Stelz, Lutz, Poethke & Bonte
// Spatially correlated extinctions select for less emigration but larger dispersal distances in the spider mite Tetranychus urticae.
// Evolution

/*
	Copyright (C) 2015  Emanuel A. Fronhofer

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
	
*/

#include <iostream>
#include <cstdlib>								//standard C library
#include <ctime>								//access system time library
#include <fstream>								//file streaming library
#include <string>								//string library included
#include <sstream>								//string streaming for reading numbers from
#include <vector>								//vectors
#include <cmath>								//standard math library
#include <gsl/gsl_rng.h>						//gsl random number generator
#include <gsl/gsl_randist.h>					//gsl random distributions
#include <gsl/gsl_statistics.h>					//gsl some statistical methods
#include <gsl/gsl_statistics_double.h> 			//gsl some double statistical methods
#include <gsl/gsl_sort_double.h> 				//gsl sort double arrays
#include <algorithm>

using namespace std;

#include "include/procedures.h"					//procedure simplifications
#include "include/classes.h"					//class definitions

//_____________________________________________________________________________
//------------------------------------------------------------ global variables
unsigned int sim_time;															// actual time in simulation
int max_runs;																	// no. of repeats
float mut_sd;																	// variance/sd used for mutations
float mut_rate;																	// mutation rate
int capacity;																	// habitat capacity
float lambda_null;																// fertility
float mu0;																		// dispersal mortality
float alpha;																	// Allee effect strength
float sigma;																	// environmental stochasticity
int epsilon;																	// random patch extinctions
bool nnd;																		// nearest neighbour dispersal (yes=true, no=false)
bool clumped_ext;																// clumped extinctions (yes/no)
TPatch world[WORLDDIM][WORLDDIM];												// simulated world
float rel_metapopsize;															// relative metapopulation size
float occupancy;																// metapopulation occupancy
float rel_emigrants;															// relative number of emigrants
float rel_addmovers;															// relative number of dispersers that do more than one step

//_____________________________________________________________________________
//------------------------------------------------------------------ procedures

//------------------------------------------------------------- read parameters
void readParameters(){
	ifstream parinfile("input/parameters.in");							//parameter input file
	string buffer;
	istringstream is;

	getline(parinfile,buffer); getline(parinfile,buffer);
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> sim_time;																		//simulation time
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> max_runs;																		//no. of repeats
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> mut_sd;																		//variance used for mutations
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> mut_rate;																		//mutation rate
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> capacity;																		//habitat capacity
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> lambda_null;																	//fertility
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> sigma;																		//uncorrelated evironmental stochasticity
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> alpha;																		//Allee effect strength
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> epsilon;																		//random patch extinction probability
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	is >> mu0;																			//dispersal mortality
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	if (buffer=="yes") {nnd = true;} else {nnd = false;};								//nearest neighbour dispersal
	getline(parinfile,buffer); getline(parinfile,buffer); is.clear(); is.str(buffer);
	if (buffer=="yes") {clumped_ext = true;} else {clumped_ext = false;};				//clumped extinctions
	parinfile.close();
}

//------------------------------------------------------- initialize simulation
void Initialize(){
	// initialize patches and individuals in patches
	for (int x = 0; x < WORLDDIM; ++x) {
		for (int y = 0; y < WORLDDIM; ++y) {
			
			// clear the world
			world[x][y].females.clear();
			world[x][y].newFemales.clear();

			// initialize individuals in this patch
			// females
			for (int f = 0; f < capacity; ++f) {
				TIndiv newfemale;
				for (int a = 0; a < N_CLASSES; ++a) {
					newfemale.dispKernel[a] = ran();
				}
				world[x][y].females.push_back(newfemale);
			}
		}
	}
}

// ------------------------------------------------ analyze population dynamics
void Analyze(unsigned int acttime, int actrun){
	//reset metapopulation size and occupancy
	rel_metapopsize = 0;
	occupancy = 0;

	unsigned int numberoccupied = 0;
	unsigned int metapopsize = 0;

	for (int x = 0; x < WORLDDIM; ++x) {
		for (int y = 0; y < WORLDDIM; ++y) {
			unsigned int localpopsize = world[x][y].females.size();
			metapopsize += localpopsize;
			if (localpopsize > 0) {
				++numberoccupied;
			}
		}
	}
	// calculate occupancy
	occupancy = float(numberoccupied) / float(WORLDDIM*WORLDDIM);
	// calculate relative metapopulation size
	rel_metapopsize = float(metapopsize) / float(WORLDDIM*WORLDDIM*capacity);
}

// ---------------------------------------------------- save individual results
void saveResults(int actrun){
	// output file: individuals
	stringstream outputindiv_path_stream;
	outputindiv_path_stream << "output/output_individuals_run" << actrun << ".out";
	string outputindiv_path = outputindiv_path_stream.str();
	ofstream outputindiv(outputindiv_path.c_str());
	// headers
	outputindiv << "x" << "    " << "y" << "    " << "sex" << "    " << "dispKernel_0" << "    " << "dispKernel_1" << "    " << "dispKernel_2" << "    " << "dispKernel_3" << "    " << "dispKernel_4" << "    " << "dispKernel_5" << "    " << "dispKernel_6" << endl;
	for (int x = 0; x < WORLDDIM; ++x) {
		for (int y = 0; y < WORLDDIM; ++y) {
			for (unsigned int f = 0; f < world[x][y].females.size(); ++f) {
				// write metapop results to file
				outputindiv << x << "    " << y << "    " << "f" << "    " << world[x][y].females.at(f).dispKernel[0]<< "    " << world[x][y].females.at(f).dispKernel[1]<< "    " << world[x][y].females.at(f).dispKernel[2]<< "    " << world[x][y].females.at(f).dispKernel[3] << "    " << world[x][y].females.at(f).dispKernel[4] << "    " << world[x][y].females.at(f).dispKernel[5] << "    " << world[x][y].females.at(f).dispKernel[6] << endl;
			}
		}
	}
	// close indiv output file
	outputindiv.close();
}

// --------------------------------------------- find new patch during dispersal
vector<int> findNewPatch(int x, int y, int no_steps){
	vector<int> res;
	// initialize res with present x and y
	res.push_back(x);
	res.push_back(y);
	// initialize dir with a starting value
	int dir = floor(ran()*4);
	// go through all steps
	// check for boundary; following experimental setting explicitly
	for (int d = 0; d < no_steps; ++d) {
		if (res.at(0) == 0){
			if (res.at(1) == 0){
				if (dir != 1 && dir != 0){
					if (dir == 2){dir = 1;}
					if (dir == 3){dir = 0;}
				}
			}
			if (res.at(1) == WORLDDIM-1){
				if (dir != 1 && dir != 2){
					if (dir==0){dir=1;}
					if (dir==3){dir=2;}
				}
			}
			if (res.at(1) != 0 && res.at(1) != WORLDDIM-1){
				if ( dir == 3){
					dir = floor(ran()*2);
					if (dir==1){dir=2;}
				}
			}
		}
		if (res.at(0) == WORLDDIM-1){
			if (res.at(1) == 0){
				if (dir !=0 && dir !=3){
					if (dir ==2){dir=3;}
					if (dir==1){dir=0;}
				}
			}
			if (res.at(1) == WORLDDIM-1){
				if (dir != 2 && dir != 3){
					if (dir ==0){dir=3;}
					if(dir==1){dir=2;}
				}
			}
			if (res.at(1) != 0 && res.at(1) != WORLDDIM-1){
				if (dir == 1){
					dir=floor(ran()*2);
					if (dir==1){dir=2;}
				}
			}
		}
		if (res.at(0) != 0 && res.at(0) != WORLDDIM-1 && res.at(1)==0){
			if ( dir == 2){
				dir = floor(ran()*2);
				if (dir == 0){dir = 3;}
			}
		}
		if (res.at(0) != 0 && res.at(0) != WORLDDIM-1 && res.at(1)==WORLDDIM-1){
			if (dir == 0){
				dir = floor(ran()*2);
				if(dir==0){dir=3;}
			}
		}
		// use this direction to make a step
		switch (dir) {
			case 0:
				res.at(1) = res.at(1) +1;
				break;
			case 1:
				res.at(0) = res.at(0)+1;
				break;
			case 2:
				res.at(1) = res.at(1)-1;
				break;
			case 3:
				res.at(0) = res.at(0)-1;
				break;
			default:
				cout << "Error in NND: wrong direction value" << endl;
				break;
			}
			if (res.at(0) >= WORLDDIM){
				cout << "Error in NND: indiv. over right boundary" << endl;
			}
			if (res.at(0) < 0){
				cout << "Error in NND: indiv. over left boundary" << endl;
			}
			if (res.at(1) >= WORLDDIM){
				cout << "Error in NND: indiv. over upper bounday" << endl;
			}
			if (res.at(1) < 0){
				cout << "Error in NND: indiv. over lower bounday" << endl;
			}
	}
	return(res);
}

// -------------------------------------------------------- dispersal procedure
void Dispersal(){

	unsigned int no_emigrants = 0;
	unsigned int no_addmovers = 0;
	unsigned int metapopsize = 0;
	rel_emigrants = 0;
	rel_addmovers = 0;
	for (int x = 0; x < WORLDDIM; ++x) {
		for (int y = 0; y < WORLDDIM; ++y) {
			metapopsize += world[x][y].females.size();
			for (unsigned int f = 0; f < world[x][y].females.size(); ++f) {
				// normalize the kernel
				float helpKernel[N_CLASSES];
				float kernelSum=0;
				for (int cl = 0; cl < N_CLASSES; ++cl) {
					kernelSum = kernelSum+world[x][y].females.at(f).dispKernel[cl];
				}
				for (int cl = 0; cl < N_CLASSES; ++cl) {
					helpKernel[cl] = world[x][y].females.at(f).dispKernel[cl]/kernelSum;
				}
				//determine distance
				float limit = ran();
				float threshold = 0;
				int no_steps=0;
				for (int cl = 0; cl < N_CLASSES; ++cl) {
					threshold = threshold + helpKernel[cl];
					if(limit <= threshold){
						no_steps = cl;
						break;
					}
				}
				//emigration
				if(no_steps > 0){
					// increase counter
					++no_emigrants;
					if(no_steps > 1) ++no_addmovers;
					// calculate mortality
					float morality = 1 - exp(-mu0*no_steps);
					// check whether the individual survives
					if (ran() > morality){
						//find new patch
						vector<int> newPatch = findNewPatch(x, y, no_steps);
						// copy disperser into new patch
						TIndiv Disperser = world[x][y].females.at(f);
						world[newPatch.at(0)][newPatch.at(1)].newFemales.push_back(Disperser);
					}
					// delete emigrant from natal patch
					world[x][y].females.at(f) = world[x][y].females.back();
					world[x][y].females.pop_back();
					// decrease female loop counter
					--f;
				}
			}
		}
	}
	// merge philopatrics and residents
	for (int x = 0; x < WORLDDIM; ++x) {
		for (int y = 0; y < WORLDDIM; ++y) {
			for (unsigned int f = 0; f < world[x][y].newFemales.size(); ++f) {
				world[x][y].females.push_back(world[x][y].newFemales.at(f));
			}
			world[x][y].newFemales.clear();
		}
	}
	rel_emigrants = float(no_emigrants) / float(metapopsize);
	rel_addmovers = float(no_addmovers) / float(no_emigrants);
}

// ------------------------------------------------------------------ mutations
float mutate(float allele){
	if(ran()< mut_rate){
		float newallele = allele + gauss(mut_sd);
		if (newallele < 0) newallele = 0;
		return(newallele);
	} else {
		return(allele);
	}
}

// ------------------------------------------------------------ larval survival
float larvalSurvival(unsigned int x, unsigned int y){
	float survival = 0;
	int n_adults = world[x][y].females.size();
	float a = (float(lambda_null) - 1) / float(capacity);
	if (a < 0) {a = 0;}
	//calculation of square-density for allee-effect
	float sqr_dens = float(n_adults)/float(capacity) * float(n_adults)/float(capacity) ;	
	//survival-probability of newborns based on logistic growth
	survival = (sqr_dens/(sqr_dens+(alpha*alpha)))/(1+a*float(n_adults));				
	return(survival);
}

// --------------------------------------------------------------- reproduction
void Reproduction(){
	for (int x = 0; x < WORLDDIM; ++x) {
		for (int y = 0; y < WORLDDIM; ++y) {
			world[x][y].newFemales.clear();
			if (world[x][y].females.size() > 0) {
				float survival = larvalSurvival(x,y);
				float lambda_local = lognorm(lambda_null, sigma);
				for (unsigned int f = 0; f < world[x][y].females.size(); ++f) {
					// calculate number of offspring
					int no_offspring = poisson(lambda_local*survival);
					// loop over offspring
					for (int o = 0; o < no_offspring; ++o) {
						// initialize new individual
						TIndiv newOffspring;
						// inherit emigration rate allele from parent
						for (int dc = 0; dc < N_CLASSES; ++dc) {
							newOffspring.dispKernel[dc] = mutate(world[x][y].females.at(f).dispKernel[dc]);
						}
						// add new individual
						world[x][y].newFemales.push_back(newOffspring);
					}
				}
			}
		}
	}
}

// -------------------------------------------------- death of annual organisms
void Death(){
	for (int x = 0; x < WORLDDIM; ++x) {
		for (int y = 0; y < WORLDDIM; ++y) {
			int local_offspring_no = world[x][y].newFemales.size();
			//cout << local_offspring_no << endl;
			if (local_offspring_no > 0) {
				// now clear adult vector
				world[x][y].females.clear();
				// now copy new females into females
				for (unsigned int nf = 0; nf < world[x][y].newFemales.size(); ++nf) {
					world[x][y].females.push_back(world[x][y].newFemales.at(nf));
				}
				// clear new females vector
				world[x][y].newFemales.clear();
			} else {
				world[x][y].females.clear();
			}
		}
	}
}

// ------------------------------------------------------------ patch extinctions
void Extinction(){

	if (clumped_ext == false){
		//sampling without replacement in order to choose blocks for extinction
		//create array from which sampling will take place
		int all[WORLDDIM*WORLDDIM];
		for (int i=0; i < (WORLDDIM*WORLDDIM); i++){
			all[i] = i;
		}
		int sample[epsilon];
		gsl_ran_choose(gBaseRand,sample,epsilon,all,WORLDDIM*WORLDDIM,sizeof(int));
		int cnt = 0;
		int cnt2 = 0;
		for (int x = 0; x < WORLDDIM; ++x) {
			for (int y = 0; y < WORLDDIM; ++y) {
				cnt = cnt +1;
				if (cnt == sample[cnt2]) {
					// extinction
					world[x][y].females.clear();
					cnt2 = cnt2 + 1;
					if (cnt2 >= epsilon) cnt2 = epsilon-1;
				}
			}
		}
	}
	if (clumped_ext == true){
		//sampling without replacement in order to choose blocks for extinction
		int all[9];
		for (int i=0; i < (9); i++){
			all[i] = i;
		}
		int sample2[1];
		int sample[epsilon];
		gsl_ran_choose(gBaseRand,sample2,1,all,WORLDDIM*WORLDDIM,sizeof(int));
		// define patches to go extinct following experimental settings
		if (sample2[0] == 0) {
			sample[0] = 0;
			sample[1] = 1;
			sample[2] = 4;
			sample[3] = 5;
		}
		if (sample2[0] == 1) {
			sample[0] = 1;
			sample[1] = 2;
			sample[2] = 5;
			sample[3] = 6;
		}
		if (sample2[0] == 2) {
			sample[0] = 2;
			sample[1] = 3;
			sample[2] = 6;
			sample[3] = 7;
		}
		if (sample2[0] == 3) {
			sample[0] = 4;
			sample[1] = 5;
			sample[2] = 8;
			sample[3] = 9;
		}
		if (sample2[0] == 4) {
			sample[0] = 5;
			sample[1] = 6;
			sample[2] = 9;
			sample[3] = 10;
		}
		if (sample2[0] == 5) {
			sample[0] = 6;
			sample[1] = 7;
			sample[2] = 10;
			sample[3] = 11;
		}
		if (sample2[0] == 6) {
			sample[0] = 8;
			sample[1] = 9;
			sample[2] = 12;
			sample[3] = 13;
		}
		if (sample2[0] == 7) {
			sample[0] = 9;
			sample[1] = 10;
			sample[2] = 13;
			sample[3] = 14;
		}
		if (sample2[0] == 8) {
			sample[0] = 10;
			sample[1] = 11;
			sample[2] = 14;
			sample[3] = 15;
		}
		int cnt = 0;
		int cnt2 = 0;
		for (int x = 0; x < WORLDDIM; ++x) {
			for (int y = 0; y < WORLDDIM; ++y) {
				cnt = cnt +1;
				if (cnt == sample[cnt2]) {
					//extinction
					world[x][y].females.clear();
					cnt2 = cnt2 + 1;
					if (cnt2 >= epsilon) cnt2 = epsilon-1;
				}
			}
		}
	}
}

//_____________________________________________________________________________
//------------------------------------------------------------------------ main

int main() {
	// initialize random number generator
	specify_rng(RS);
	//read parameters for all simulation runs
	readParameters();
	// repeat loop
	for (int actrun = 0; actrun < max_runs; ++actrun) {
		// output file: metapopulation
		stringstream outputmetapop_path_stream;
		outputmetapop_path_stream << "output/output_metapop_run" << actrun << ".out";
		string outputmetapop_path = outputmetapop_path_stream.str();
		ofstream outputmetapop(outputmetapop_path.c_str());
		// outputfile header
		outputmetapop << "time" << "    " << "rel_metapopsize" << "    " << "occupancy"<< "    " << "emirate" << endl;
		Initialize();
		// time loop
		for (unsigned int acttime = 0; acttime < sim_time; ++acttime) {
			// dispersal
			Dispersal();
			// reproduction
			Reproduction();
			// death of adults
			Death();
			// random patch extinctions
			Extinction();
			// analyze metapopulation
			Analyze(acttime, actrun);
			// write metapop results to file
			outputmetapop << acttime << "    " << rel_metapopsize << "    " << occupancy << "    " << rel_emigrants << endl;
			// write run and time to console
			cout << actrun << "   " << acttime << "   " << rel_metapopsize << "    " << occupancy << "   " << rel_emigrants << "   " << rel_addmovers <<  endl;
			// save individual results once at simulation end
			if (acttime == sim_time-1){
				saveResults(actrun);
			}
		}
		// close metapop output file
		outputmetapop.close();
	}
	return 0;
}
