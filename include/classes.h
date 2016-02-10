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

//____________________________________________________________________________
//----------------------------------------------------------- define constants

const int WORLDDIM = 4;														// world dimensions, always quadratic
const int RS = 2;															// random seed
const int N_CLASSES = 7;


//____________________________________________________________________________
//------------------------------------------------------------- define classes

// one individual ------------------------------------------------------------
class TIndiv {
public:
	TIndiv();
	float dispKernel[N_CLASSES];

};

TIndiv::TIndiv() { //constructor for TIndiv
    for (int i = 0; i < N_CLASSES; ++i) {
		dispKernel[i] = 0;
	}
}

// one patch -----------------------------------------------------------------
class TPatch {
public:
	TPatch();
	vector<TIndiv> females;
	vector<TIndiv> newFemales;
};

TPatch::TPatch() {
	females.clear();
	newFemales.clear();
}
