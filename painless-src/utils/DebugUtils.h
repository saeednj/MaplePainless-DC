// -----------------------------------------------------------------------------
// Copyright (C) 2017  Ludovic LE FRIOUX
//
// This file is part of PaInleSS.
//
// PaInleSS is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see <http://www.gnu.org/licenses/>.
// -----------------------------------------------------------------------------

#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector>

using namespace std;

//DEBUG MACROS
//#define DEBUG(X) X
#define DEBUG(X)

/// Print array of a given length.
inline void printArray(int * s, int len)
{
	for (int i = 0; i < len; i++) {
		printf("%d ", s[i]);
	}

	printf("\n");
}

/// Print array of length parts*size, each part on a line.
inline void printArray(int * s, int parts, int size)
{
	for (int part = 0; part < parts; part++) {
		printf("Part %d:", part);

		for (int j = 0; j < size; j++) {
			printf("%d ", s[j+(part*size)]);
		}

		printf("\n");
	}
}

/// Print a clause (vector of integers).
inline void printVector(const vector<int> & vec)
{
	printf("[");

	if (vec.size() > 0) {
		printf("%d",vec[0]);
	}

	for (size_t i = 1; i < vec.size(); i++) {
		printf(", %d", vec[i]);
	}

	printf("]\n");
}

/// Print zero-terminated array.
inline void printArrayZT(int * s, int sep = 0)
{
	while(*s != sep) {
		printf("%d ", *s);
		s++;
	}

	printf("\n");
}
