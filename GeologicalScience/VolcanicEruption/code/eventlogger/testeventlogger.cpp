/*
Copyright (C) 2018 Pierre Kunzli, University of Geneva, Switzerland

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

#include "eventlogger.hpp"
#include <iostream>
#include <math.h>


inline double func(double val){
  return sqrt(sqrt(val)*pow(val, 2.0))/(log10(val)*500.0);
}

inline double iterativeFunc(double val){
  double sum = val;
  for(int i=0; i<1000; i++) sum += func(sum);
  return sum;
}

void testEventLogger(unsigned long long events, EventLogger &logger){
    double res = 0.0;
    logger.log("startloop");
    for(unsigned long long i=0; i<events; i++){
      logger.log("startfunc");
      res += iterativeFunc(100.0);
    }
    logger.log("endloop");
}

void test(unsigned long long events){
    double res = 0.0;
    for(unsigned long long i=0; i<events; i++){
      res += iterativeFunc(100.0);
    }
}

void printUsage(){
  std::cout << "Usage : testeventlogger rentention events" << std::endl;
}

void parametersParseFailure(){
  printUsage();
  exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
{

  unsigned long long retention;
  unsigned long long events;

  if(argc != 3) parametersParseFailure();

  try{
    retention = std::stoi(argv[1]);
    events = std::stoi(argv[2]);
  } catch(const std::exception& e){
    parametersParseFailure();
  }

  std::string filenameGlobal = "log_logger.log";
  std::string filenameLocal = "log_func_"+std::to_string(retention)+".log";
  remove(filenameGlobal.c_str());
  remove(filenameLocal.c_str());

  EventLogger globalEventLogger(filenameGlobal, 10);
  EventLogger localEventLogger(filenameLocal, retention);

  globalEventLogger.log("starttesteventlogger");
  testEventLogger(events, localEventLogger);
  globalEventLogger.log("startttest");
  test(events);
  globalEventLogger.log("endtest");

  exit(EXIT_SUCCESS);
}