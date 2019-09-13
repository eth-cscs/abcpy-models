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


#include "musclehpc/utils/optionparser.h"




OptionParser*  OptionParser::instance = nullptr;

OptionParser::OptionParser(int argc, char **argv):
    argc(argc), argv(argv), isAllKerneltorall(false), runManagerHere(false),
    isAlreadyParsed(false),useTCP_Manager(false), generateProfiling(false)
{
    this->parse();
    this->constructCoresMapRequest();
}

OptionParser::~OptionParser(){}


int OptionParser::getArgc() const {return this->argc;}
char ** OptionParser::getArgv() const{return this->argv;}

string OptionParser::getManagerUrl() const{
    return this->managerUri;
}

string OptionParser::getForwarderFrontEnd() const{
    return this->forwarderFrontEnd;
}

string OptionParser::getManagerOutFile() const{
    return this->managerOutFile;
}

vector<string> OptionParser::getkernelsToRun() const{
    return this->kernelsTorun;
}

bool OptionParser::allowAllKernelToRun() const{
    return this->isAllKerneltorall;
}
bool OptionParser::isManager() const{
    return this->runManagerHere;
}

bool OptionParser::useManagerTCP() const{
    return this->useTCP_Manager;
}

bool OptionParser::isProfiling() const{
    return this->generateProfiling;
}

string OptionParser::getDisplayHelp() const{
    return this->displayHelp;
}



void OptionParser::constructCoresMapRequest()
{
    std::ifstream file(this->cxaFile);
    std::string str;
    if(!file){
        //cout<<"Can not open CXA config file: "<<this->cxaFile<<endl;
    }else{
        while (std::getline(file, str))
        {
            RegExprCoupling &coupler = RegExprCoupling::getInstance();
            int cores;
            string kname;
            coupler.parseCoresRequest(str, kname, cores);
            if(!kname.empty()){
                cores_map.insert(std::make_pair(kname, cores));
            }
        }
    }
}

map<string, int> const&  OptionParser::getRequestedCores()
{
    return cores_map;
}

bool OptionParser::getCouplingFromCxaFile(vector<CouplingInfo*> & vectInfo){
    std::ifstream file(this->cxaFile);
    vectInfo.clear();
    std::string str;
    if(!file){
        //cout<<"Can not open CXA config file: "<<this->cxaFile<<endl;
        return false;
    }else{
        while (std::getline(file, str))
        {
            RegExprCoupling &coupler = RegExprCoupling::getInstance();
            CouplingInfo* info=coupler.parseLine(str);
            if(info)
                vectInfo.push_back(info);
        }
    }
    return true;

}
bool OptionParser::getCouplingFromCxaFile(vector<CouplingInfo*> & vectInfo, vector<CommandLineInfo*> & vectCmdInfo){
    std::ifstream file(this->cxaFile);
    vectInfo.clear();
    std::string str;
    if(!file){
        //cout<<"Can not open CXA config file: "<<this->cxaFile<<endl;
        return false;
    }else{
        string programName = this->argv[0];
        RegExprCoupling &coupler = RegExprCoupling::getInstance();
        while (std::getline(file, str))
        {
            CouplingInfo* info=coupler.parseLine(str);
            CommandLineInfo * cmdInfo= coupler.parseCommandLine(str, programName);
            if(info){
                vectInfo.push_back(info);
            }
            if (cmdInfo){
                vectCmdInfo.push_back(cmdInfo);
            }
        }
    }
    return true;
}



void OptionParser::parse(){
    if(isAlreadyParsed)
        return;
    isAlreadyParsed=true;

    //int opt;
    //std::string ConfigFile;

    // Declare a group of options that will be
    // allowed only on command line
    po::options_description generic("Generic options");
    generic.add_options()
            ("version,v", "print version string")
            ("help,h", "produce help message")
            // ("config,c",po::value<std::string>(&ConfigFile)->default_value("multiple_sources.cfg"), "Name of a file of a configuration.")
            ;

    // Declare a group of options that will be
    // allowed both on command line and in
    // config file
    po::options_description config("Configuration");
    config.add_options()
            ("cxa-file,c", po::value<string>(&cxaFile), "MMSF CXA configuration file to execute")
            ("forwarder,f", po::value<string>(&forwarderFrontEnd), "Forwarder front-end server address. Syntax: IP:port")
            ("tcp,t", "run Manager with TCP socket Endpoint (low performance compared to MPI (default))")
            ("allkernels,a", "launch all kernels")
            ("main,m", "run  the Manager here")
            ("profiling,p", "generate profiling in kernelName_timestamp.profiler file")
            ("manager,M", po::value<string>(&managerUri), "url of  the Manager")
            ("outFileUrl,o", po::value<string>(&managerOutFile), "file where to sdore the url of the Manager (should be in a shared partition between all computing nodes)")
            ("kernel-name,k", po::value< vector<string> >(&kernelsTorun), "list of names of kernels to launch here")

            ;

    po::options_description ignored("Ignored options");
    ignored.add_options()
            ("extra-options", po::value< vector<string> >(), "list of ignored options")
            ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    /* po::options_description hidden("Hidden options");
    hidden.add_options()
            ("kernels-name", po::value< vector<string> >(), "list of names of kernels to launch here")
            ;*/


    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config).add(ignored);

    po::options_description config_file_options;
    config_file_options.add(config);//.add(hidden);

    po::options_description visible("Allowed options");
    visible.add(generic).add(config);



    po::positional_options_description p;
    // p.add("kernels-name", -1);
    p.add("extra-options", -1);

    try{
        store(po::command_line_parser(argc, argv).
              options(cmdline_options).positional(p).allow_unregistered().run(), vm);//.allow_unregistered()
        notify(vm);

        if (vm.count("help")) {
                stringstream ss;
                ss<< visible << "\n";
                this->displayHelp = ss.str();
        }
        if (vm.count("version")) {
            cout << "Muscle-HPC, CUI/UNIGE, version 1.0\n";
        }
        if (vm.count("tcp"))
            this->useTCP_Manager=true;
        else
            this->useTCP_Manager=false;

        if (vm.count("allkernels"))
            isAllKerneltorall=true;
        else
            isAllKerneltorall=false;

        if (vm.count("main"))
            runManagerHere=true;
        else
            runManagerHere=false;

        if (vm.count("profiling"))
            this->generateProfiling=true;
        else
            this->generateProfiling=false;

    } catch( const std::exception& e ) { // reference to the base of a polymorphic object
        std::cout<<"Parsing:" << e.what()<<endl; // information from length_error printed
    }
    //vector<string> to_pass_further = collect_unrecognized(parsed.options, include_positional);
    /*ifstream ifs(ConfigFile.c_str());
    if(!ifs)
    {
        cout<<"can not open config file\n";
        return 0;
    }
    else
    {
        store(parse_config_file(ifs, config_file_options), vm);
        notify(vm);
    }*/

    /* if (vm.count("help")) {
        cout << visible << "\n";
        return 0;
    }

    if (vm.count("version")) {
        cout << "Multiple sources example, version 1.0\n";
        return 0;
    }

    if (vm.count("include-path"))
    {
        cout << "Include paths are: "
             << vm["include-path"].as< vector<string> >() << "\n";
    }

    if (vm.count("kernels-name"))
    {
        cout << "Input files are: "
             << vm["kernels-name"].as< vector<string> >() << "\n";
    }
*/
}



//******************************************************* RegExprCoupling **********************************
CouplingInfo *RegExprCoupling::parseLine(string input){
    input.erase(std::remove(input.begin(),input.end(),' '),input.end());
    string s=input;
    static const boost::regex connectionExpr("(\\w+<\\w+>:\\w+\\.\\w+->\\w+\\.\\w+)");

    boost::match_results<std::string::const_iterator> results;
    if (boost::regex_match(s, results, connectionExpr))
    {
        boost::regex expr("(\\w+)");
        boost::regex_token_iterator<std::string::iterator> it{s.begin(), s.end(), expr, 1};
        //boost::regex_token_iterator<std::string::iterator> end;
        CouplingInfo * info = new CouplingInfo(*it++,*it++,*it++,*it++,*it++,*it++);
        return info;

    }else{
        return 0;
    }
}

CommandLineInfo * RegExprCoupling::parseCommandLine(string cmdLine, string programName){
    static const boost::regex commandExpr("(cmdline<\\w+>:.*)");
    boost::match_results<std::string::const_iterator> results;
    if (boost::regex_match(cmdLine, results, commandExpr))
    {
        std::vector<std::string> elems;
        boost::split(elems, cmdLine, boost::is_any_of(":"));  // Note this is boost::split
        assert(elems.size() >1);
        string s1= elems.at(0);
        string kernelCmdLine= elems.at(1);
        boost::regex expr("([^cmdline<]\\w*)");
        boost::regex_token_iterator<std::string::iterator> it{s1.begin(), s1.end(), expr, 1};
        boost::regex_token_iterator<std::string::iterator> end;
        string kernelName;
        if(it != end){
            kernelName=*it;
        }
        assert(!kernelName.empty());
        //std::cout<<"["<<kernelName<<"]>: "<< cmdLine<<endl;
        CommandLineInfo * cmdParser= new CommandLineInfo( kernelName, kernelCmdLine, programName);
        cmdParser->parse();
        return cmdParser;
    }else{
        //std::cerr<<"["<<cmdLine<<"] not a valid (cmdline<\\w+>:.*)"<<endl;
        return 0;
    }
}


    void RegExprCoupling::parseCoresRequest(string line, string & kname, int & cores){

        //map<string, int> coresmap;


        static const boost::regex commandExpr("(cores<\\w+>:\\s*[0-9]+)");
        boost::match_results<std::string::const_iterator> results;
        if (boost::regex_match(line, results, commandExpr))
        {
            std::vector<std::string> elems;
            boost::split(elems, line, boost::is_any_of(":"));  // Note this is boost::split
            assert(elems.size() >1);
            kname= elems.at(0);
            string values= elems.at(1);
            values.erase(std::remove(values.begin(),values.end(),' '),values.end());
            boost::regex expr("([^cores<]\\w*)");
            boost::regex_token_iterator<std::string::iterator> it{kname.begin(), kname.end(), expr, 1};
            boost::regex_token_iterator<std::string::iterator> end;
            if(it != end){
                kname=*it;
            }
            assert(!kname.empty());
            std::string::size_type sz;   // alias of size_t
            cores = std::stoi (values,&sz);
        }else{
            //std::cerr<<"["<<cmdLine<<"] not a valid (cmdline<\\w+>:.*)"<<endl;
            cores= 0;
            kname="";
        }

    }

//******************************************************* CommandLineInfo **********************************

CommandLineInfo::CommandLineInfo(string kernelName, string commandLine, string programName):
    kernelName(kernelName), commandLine(commandLine), programName(programName){
}

CommandLineInfo::~CommandLineInfo(){
    delete[] argv;
}


char **CommandLineInfo::getargv(){
    return this->argv;
}

int CommandLineInfo::getargc(){
    return this->argc;
}


string CommandLineInfo::getKernelName() const{
    return this->kernelName;
}

void CommandLineInfo::parse(){
    convertCommandLineToArgvStyle(this->commandLine);
}

void CommandLineInfo::convertCommandLineToArgvStyle (string commandLine){
    std::vector<std::string>  arguments;
    convertCommandLineToArguments (	commandLine, arguments);
    this->argv = new char*[arguments.size()];
    for(size_t i=0; i<arguments.size(); i++){
        std::string& str = arguments.at(i);
        char * buffer= new char[str.length() + 1];
        std::copy(str.c_str(), str.c_str()+str.length()+1, buffer);
        this->argv[i] =buffer;
    }

    this->argc = arguments.size();
    //std::transform(_arguments.begin(), _arguments.end(), myargv.begin(), [](std::string& str){return str.c_str();});
}

void CommandLineInfo::convertCommandLineToArguments (string commandLine, std::vector<std::string> & arguments ){

    // push the programName
    arguments.push_back(this->kernelName);
    boost::regex word_regex( "(\"[^\"]+\"|[^\\s\"]+)" );
    auto words_begin =  boost::sregex_iterator(commandLine.begin(), commandLine.end(), word_regex);
    auto words_end = boost::sregex_iterator();

    for (boost::sregex_iterator i = words_begin; i != words_end; ++i)
    {
        string s= i->str();
        // Remove all double-quote characters
        s.erase( remove( s.begin(), s.end(), '\"' ), s.end() );
        //std::cout << s  << '\n';
        arguments.push_back(s);
    }
    /*for (std::vector<string>::iterator it = arguments.begin();  it !=  arguments.end(); ++it){
        std::cout <<"->"<< *it  << endl;
    } */

}
