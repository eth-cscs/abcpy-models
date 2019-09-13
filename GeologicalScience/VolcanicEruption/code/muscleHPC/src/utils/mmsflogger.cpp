/**
* @author  Mohamed Ben Belgacem <Mohamed.BenBelgacem@gmail.com>

* MUSCLE-HPC communication module
* Copyright (C) 2016  University of Geneva, Switzerland
*
* MUSCLE-HPC is free software: you can redistribute it and/or
* modify it under the terms of the GNU Affero General Public License as
* published by the Free Software Foundation, either version 3 of the
* License, or (at your option) any later version.
*
* The library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Affero General Public License for more details.
*
* You should have received a copy of the GNU Affero General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "musclehpc/utils/mmsflogger.h"

/* *************** Class plb_ofstream ******************************** */

plb_ofstream::plb_ofstream(int rank)
    : rank (rank), devNullStream(&devNullBuffer),
      original (
           (rank==0) ? new std::ofstream : 0 )
{ }

plb_ofstream::plb_ofstream(int rank, const char* filename, std::ostream::openmode mode)
    : rank(rank), devNullStream(&devNullBuffer),
      original (
          (rank==0) ?
              new std::ofstream(filename,mode) : 0 )
{ }

plb_ofstream::plb_ofstream(plb_ofstream const& rhs)
    : rank(rhs.rank), devNullStream(&devNullBuffer),
      original(0)
{ }

plb_ofstream& plb_ofstream::operator=(plb_ofstream const& rhs) {
    return *this;
}



plb_ofstream::~plb_ofstream() {
    delete original;
}

std::ostream& plb_ofstream::getOriginalStream()
{
    if (rank==0) {
        return *original;
    }
    else {
        return devNullStream;
    }
}

plb_ofstream &plb_ofstream::getWriter()
{
    return *this ;
}

bool plb_ofstream::is_open() {
//#ifdef PLB_MPI_PARALLEL
    int open = false;
    if (rank==0) {
        open = original->is_open();
    }
    //this->mpi->bCast(&open, 1);
    return open;
//#else
//    return original->is_open();

//#endif
}

void plb_ofstream::open(const char* filename, std::ostream::openmode mode)
{
    if (rank==0) {
        original->open(filename, mode);
    }
}

void plb_ofstream::close() {
   if (rank==0) {
        original->close();
    }
}
void plb_ofstream::flush() {
   if (rank==0) {
        original->flush();
    }
}
