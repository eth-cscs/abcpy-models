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


#ifndef OPTIONPARSER_H
#define OPTIONPARSER_H

#include <regex>
#include <string>

#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <boost/program_options/config.hpp>
#include <boost/program_options/option.hpp>
#include <boost/program_options/detail/cmdline.hpp>
#include <boost/algorithm/string.hpp>

#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include <algorithm>
namespace po = boost::program_options;

#include <iostream>
#include <fstream>
#include <iterator>
using namespace std;


//**************************************************** CouplingInfo **************************************

struct CouplingInfo{

private:
    string conduitName;
    string dataType;
    string senderName;
    string senderPort;
    string receiverName;
    string receiverPort;

public:
    CouplingInfo(string conduitName, string dataType, string  senderName, string senderPort, string  receiverName, string receiverPort):
        conduitName(conduitName), dataType(dataType),
        senderName(senderName), senderPort(senderPort),
        receiverName(receiverName),receiverPort(receiverPort){

    }

    string getConduitName()const{
        return this->conduitName;
    }
    string getDataType()const{
        return this->dataType;
    }
    string getSenderName()const{
        return this->senderName;
    }
    string getSenderPort()const{
        return this->senderPort;
    }
    string getReceiverName()const{
        return this->receiverName;
    }
    string getReceiverPort()const{
        return this->receiverPort;
    }

    string toString(){
        stringstream ss;
        ss<<conduitName<<"<"<<dataType<<">:"<<senderName<<"."<<senderPort<<"->"<<receiverName<<"."<<receiverPort;
        return ss.str();
    }
};

//**************************************************** CommandLineInfo **************************************

class CommandLineInfo{

private:
    string kernelName;
    string commandLine;
    string programName;
    int argc;
    char ** argv;

public:
    /**
     * @brief CommandLineInfo
     * @param kernelName name of kernel
     * @param commandLine the command line of to run by a kernel
     * @param programName the binary name to be obtained from the  argv[0] in the main method
     */
    CommandLineInfo(string kernelName, string commandLine, string programName);
    virtual~CommandLineInfo();
    // converts commandLine into argc & argv C style.
    void parse();
    // gets the char ** argv
    char ** getargv();
    //gets the int argc
    int getargc();
    string getKernelName() const;

private:
    //
    /**
     * @brief convertCommandLineToArguments splits a command line "strArg" into vector of arguments
     * @param strArg the command line of to run by a kernel
     * @param arguments a empty vector which will stotes the arguments
     */
    void convertCommandLineToArguments (string commandLine, std::vector<std::string> & arguments );
    /**
     * @brief convertCommandLineToArgvStyle  splits a command line "strArg" into argc & argv C style.
     * @param strArg the command line of to run by a kernel
     */
    void convertCommandLineToArgvStyle (string commandLine);

};

//**************************************************** Singleton RegExprCoupling **************************************

class RegExprCoupling{


public:
    static RegExprCoupling& getInstance() {
        static RegExprCoupling    instance; // Guaranteed to be destroyed.
        // Instantiated on first use.
        return instance;
    }
private:
    RegExprCoupling() {}                  // Constructor? (the {} brackets) are needed here.
    // C++ 11
    // =======
    // We can use the better technique of deleting the methods
    // we don't want.
    RegExprCoupling(RegExprCoupling const&) ;
    void operator=(RegExprCoupling const&)  ;

public:
    CouplingInfo *parseLine(string input);
    CommandLineInfo *parseCommandLine(string cmdLine, string programName);
    void parseCoresRequest(string line, string &kname, int &cores);
};




//******************************************************** Option Parser ************************************
class OptionParser
{
private:
    int argc;
    char **argv;
    po::variables_map vm;
    //
    string cxaFile;
    string managerUri;
    string managerOutFile;
    string forwarderFrontEnd;
    vector<string> kernelsTorun;
    bool isAllKerneltorall;    //
    bool runManagerHere;    //
    typedef vector<string> string_vector;
    //
    bool isAlreadyParsed;
    bool useTCP_Manager;
    bool generateProfiling;
    string displayHelp;

    map<string, int> cores_map;
    static OptionParser * instance;

public:
    OptionParser(int argc, char** argv);
    virtual ~OptionParser();

    static OptionParser *getInstance(int argc, char **argv)
    {
        if (! instance){
            instance = new OptionParser(argc, argv);
        }
        return instance;
    }

    static OptionParser *getInstance(){return instance;}

    OptionParser(OptionParser const&)    = delete;
    void operator=(OptionParser const&)  = delete;

    int getArgc() const;
    char ** getArgv ()const;

    string getManagerUrl() const;
    string getForwarderFrontEnd() const;
    string getManagerOutFile() const;
    vector<string> getkernelsToRun() const ;
    bool allowAllKernelToRun() const;
    bool isManager() const;
    bool useManagerTCP() const;
    bool isProfiling() const;
    string getDisplayHelp() const;

    const map<string, int> &getRequestedCores();
    bool getCouplingFromCxaFile(vector<CouplingInfo*> & vectInfo);
    bool getCouplingFromCxaFile(vector<CouplingInfo*> & vectInfo, vector<CommandLineInfo*> & vectCmdInfo);
    void parse();

private:
    void constructCoresMapRequest();

};


#endif // OPTIONPARSER_H
