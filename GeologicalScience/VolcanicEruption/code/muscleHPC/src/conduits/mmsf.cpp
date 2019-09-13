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

#include "mmsf.hpp"

template Binder<char> * MMSF::connect( string fromStr ) ;
template Binder<int> *  MMSF::connect( string fromStr ) ;
template Binder<long> *  MMSF::connect( string fromStr ) ;
template Binder<long long> *  MMSF::connect( string fromStr ) ;
template Binder<unsigned long long> *  MMSF::connect( string fromStr ) ;
template Binder<float> *  MMSF::connect( string fromStr ) ;
template Binder<double> *  MMSF::connect( string fromStr ) ;
template Binder<long double> *  MMSF::connect( string fromStr ) ;

template void  MMSF::connect( string fromStr, string toStr, Conduit<char> * conduit, string conduitID );
template void  MMSF::connect( string fromStr, string toStr, Conduit<int> * conduit, string conduitID );
template void  MMSF::connect( string fromStr, string toStr, Conduit<long> * conduit, string conduitID );
template void  MMSF::connect( string fromStr, string toStr, Conduit<long long> * conduit, string conduitID );
template void  MMSF::connect( string fromStr, string toStr, Conduit<unsigned long long> * conduit, string conduitID );
template void  MMSF::connect( string fromStr, string toStr, Conduit<float> * conduit, string conduitID );
template void  MMSF::connect( string fromStr, string toStr, Conduit<double> * conduit, string conduitID );
template void  MMSF::connect( string fromStr, string toStr, Conduit<long double> * conduit, string conduitID );

template Kernel * MMSF::assignConduit(string portString, Conduit<char> * conduit, bool isIn);
template Kernel * MMSF::assignConduit(string portString, Conduit<int> * conduit, bool isIn);
template Kernel * MMSF::assignConduit(string portString, Conduit<long> * conduit, bool isIn);
template Kernel * MMSF::assignConduit(string portString, Conduit<long long> * conduit, bool isIn);
template Kernel * MMSF::assignConduit(string portString, Conduit<unsigned long long> * conduit, bool isIn);
template Kernel * MMSF::assignConduit(string portString, Conduit<float> * conduit, bool isIn);
template Kernel * MMSF::assignConduit(string portString, Conduit<double> * conduit, bool isIn);
template Kernel * MMSF::assignConduit(string portString, Conduit<long double> * conduit, bool isIn);
