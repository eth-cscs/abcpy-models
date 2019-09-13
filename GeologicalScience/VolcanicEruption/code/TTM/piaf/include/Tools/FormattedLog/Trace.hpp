/*
Particles in Advection Field (PIAF)
Copyright (C) 2015  University of Geneva, Switzerland

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef TRACE_H_
#define TRACE_H_

#include <sstream>
#include <boost/format.hpp>
#include <iostream>
#include <fstream>
#include <string>

#include "boost/date_time/posix_time/posix_time.hpp"

//namespace piaf{

//using namespace std;

// Logging class, taken from :
// http://stackoverflow.com/questions/1255576/what-is-good-practice-for-generating-verbose-output/1255904#1255904


namespace TraceImpl {

class Trace {
public:
	Trace( const char* msg, const char* filename ) : fmt(msg), fname( std::string( filename ) ) {}
	~Trace() {
		std::ofstream file;
		file.open (fname,std::ios::app);
		file<<boost::posix_time::microsec_clock::local_time()<<" : "<<fmt<<std::endl<<std::endl;
		file.close();

	}
	template <typename T>
	Trace& operator %(T value) {
		fmt % value;
		return *this;
	}
protected:
	boost::format fmt;
	std::string fname;
};
}

inline TraceImpl::Trace trace( const char *msg, const char *filename ) {
	return TraceImpl::Trace( msg, filename );
}

//} // namespace piaf

#endif
