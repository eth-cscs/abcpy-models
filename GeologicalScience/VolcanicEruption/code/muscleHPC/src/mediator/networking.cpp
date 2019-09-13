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

#include "musclehpc/mediator/networking.h"

Networking::Networking():chunk_size(4096){}

void Networking::getHostIpAddresses(vector<string> & ipAdresses, bool isIPV4, bool isIPV6){
	struct ifaddrs *ifaddr, *ifa;
	int family, s, n;
	char host[NI_MAXHOST];

	if (getifaddrs(&ifaddr) == -1) {
		perror("getifaddrs");
		exit(EXIT_FAILURE);
	}

	/* Walk through linked list, maintaining head pointer so we
	   can free list later */

	for (ifa = ifaddr, n = 0; ifa != NULL; ifa = ifa->ifa_next, n++) {
		if (ifa->ifa_addr == NULL)
			continue;

		family = ifa->ifa_addr->sa_family;

		/* Display interface name and family (including symbolic
		   form of the latter for the common families) */

		/* printf("%-8s %s (%d)\n",
		   ifa->ifa_name,
		   (family == AF_PACKET) ? "AF_PACKET" :
		   (family == AF_INET) ? "AF_INET" :
		   (family == AF_INET6) ? "AF_INET6" : "???",
		   family);*/

		/* For an AF_INET* interface address, display the address */

		if (family == AF_INET || family == AF_INET6) {

			s = getnameinfo(ifa->ifa_addr,
					(family == AF_INET) ? sizeof(struct sockaddr_in) : sizeof(struct sockaddr_in6),
					host, NI_MAXHOST,
					NULL, 0, NI_NUMERICHOST);
			if (s != 0) {
				printf("getnameinfo() failed: %s\n", gai_strerror(s));
				//exit(EXIT_FAILURE);
			}

			//printf("\t\taddress: <%s>\n", host);
			string ip(host);
			if(isIPV4 && family == AF_INET)
				ipAdresses.push_back(ip);
			if(isIPV6 && family == AF_INET6)
				ipAdresses.push_back(ip);

		} else if (family == AF_PACKET && ifa->ifa_data != NULL) {
			/*struct rtnl_link_stats *stats = (rtnl_link_stats *)ifa->ifa_data;

			  printf("\t\ttx_packets = %10u; rx_packets = %10u\n"
			  "\t\ttx_bytes   = %10u; rx_bytes   = %10u\n",
			  stats->tx_packets, stats->rx_packets,
			  stats->tx_bytes, stats->rx_bytes);*/
		}
	}

	freeifaddrs(ifaddr);
}

int Networking::connectToFrontEndRelayer( Endpoint & myrelayer){
	for (auto socketurl: myrelayer.urls){
		int sockfd=connectToFrontEndSocket(socketurl);
		if (sockfd>0){
			return sockfd;
		}
	}
	return -1;
}

int Networking::connectToFrontEndForwarder(Endpoint &forwarder, string name, string remoteEndpoint, string correspondingForwarding){

	int sock=this->connectToFrontEndRelayer(forwarder);
	if (sock<0){
        //throw std::invalid_argument( "can not connect to the forwarder: "+forwarder.rawEndpoints );
        return -1;
	}

	int needManager=1;
	int len=name.size()+1;
	int len1=remoteEndpoint.size()+1;
    int len2=correspondingForwarding.size()+1;
    int buffsize=sizeof(int)*4+len+len1+len2;
	vector<char>  buffer;
	buffer.resize(buffsize);
	int idData=0;
	memcpy((char*) &buffer[idData], (const char*)&needManager, sizeof(needManager));
	idData+=sizeof(needManager);
	memcpy((char*)  &buffer[idData], (const char*)&len, sizeof(len));
	idData+=sizeof(len);
	std::copy(name.c_str(), name.c_str()+len,  &buffer[idData]);
	idData+=len;
	memcpy((char*)  &buffer[idData], (const char*)&len1, sizeof(len1));
	idData+=sizeof(len1);
	std::copy(remoteEndpoint.c_str(), remoteEndpoint.c_str()+len1,  &buffer[idData]);
    idData+=len1;
    memcpy((char*)  &buffer[idData], (const char*)&len2, sizeof(len2));
    idData+=sizeof(len2);
    std::copy(correspondingForwarding.c_str(), correspondingForwarding.c_str()+len2,  &buffer[idData]);

	Networking::getInstance().sock_write(sock, buffer.data(), buffer.size());
	vector<char>().swap(buffer);
	return sock;

}

int Networking::connectToFrontEndSocket(SocketURL & frontEndUrl){

    int sockfd;
	struct addrinfo hints, *servinfo, *p;
	int rv;
	char s[INET6_ADDRSTRLEN];

	memset(&hints, 0, sizeof hints);
	hints.ai_family = AF_UNSPEC;
	hints.ai_socktype = SOCK_STREAM;
    stringstream tmp; tmp<<frontEndUrl.port;
	const char * chport = tmp.str().c_str();
	if ((rv = getaddrinfo(frontEndUrl.ipAddress.c_str(), chport, &hints, &servinfo)) != 0) {
		cerr<<"getaddrinfo: "<< gai_strerror(rv) <<endl;
		return -1;
	}

    //this->setnonblocking(sockfd);
    int connectionResult=-1;
	// loop through all the results and connect to the first we can
	for(p = servinfo; p != NULL; p = p->ai_next) {

		if ((sockfd = socket(p->ai_family, p->ai_socktype, p->ai_protocol)) == -1) {
			perror("client: socket");
			continue;
		}

        // strat trying connecting
        int ttr=2;
		do{
            struct sockaddr *sa = (struct sockaddr *)p->ai_addr;
            void * res;
            if (sa->sa_family == AF_INET) {
                res= &(((struct sockaddr_in*)sa)->sin_addr);
            }else{
                res= &(((struct sockaddr_in6*)sa)->sin6_addr);
            }

            inet_ntop(p->ai_family, res ,s, sizeof s);

			stringstream ss;
            ss<<"try connecting to: "<< s<<":"<<frontEndUrl.port;
            connectionResult= connect(sockfd, p->ai_addr, p->ai_addrlen);
			if ( connectionResult < 0){
                ss<<" --> failed to connect"<<endl;
				cout<<ss.str();
			}else{
                ss<<" --> connection ok"<<endl;
                cout<<ss.str();
				break;
			}
            //sleep(3);
        } while ( (connectionResult < 0) && (ttr-- >0) );

        if (connectionResult <0) {
            cout<<"client: connection refused to "<<s<<endl;
			close(sockfd);
			continue;
        }

        break;
	}

	if (p == NULL) {
		cerr <<"client: failed to connect"<<endl;;
		return -1;
	}

	struct sockaddr *sa = (struct sockaddr *)p->ai_addr;
	void * res;
	if (sa->sa_family == AF_INET) {
		res= &(((struct sockaddr_in*)sa)->sin_addr);
	}else{
		res= &(((struct sockaddr_in6*)sa)->sin6_addr);
	}

	inet_ntop(p->ai_family, res ,s, sizeof s);
    cout<<"client: is connected to "<< s<<endl;
	freeaddrinfo(servinfo); // all done with this structure


	// -- set blocking --
	//this->setBlocking(sockfd);


	if(connectionResult < 0)
		sockfd=-1;
	return sockfd;

}

int Networking::openfrontEndSocket(Endpoint & endPoint, int rank){

	int frontEnd_sockfd=-1;

	string IPAddress;
	int portno=-1;

	int server_len;

	if(! endPoint.urls.empty()){
		SocketURL sockurl= endPoint.urls.at(0);
		IPAddress=sockurl.ipAddress;
		portno = sockurl.port;
	}

	//if(this->myMpiManager->isMainProcessor()){
	struct sockaddr_in server_address;
	struct hostent *server;

	//Create and name a socket for the server
	frontEnd_sockfd = socket(AF_INET, SOCK_STREAM, 0);
	if (frontEnd_sockfd < 0) {
		perror("socket");
		exit(EXIT_FAILURE);
	}
	int reuse_addr = 1;  /* Used so we can re-bind to our port
				while a previous connection is still
				in TIME_WAIT state. */
	/* So that we can re-bind to it without TIME_WAIT problems */
	setsockopt(frontEnd_sockfd, SOL_SOCKET, SO_REUSEADDR, &reuse_addr, sizeof(reuse_addr));

	/* Set socket to non-blocking with our setnonblocking routine */
	setnonblocking(frontEnd_sockfd);

	server_address.sin_family = AF_INET;
	// server_address.sin_addr.s_addr = *(IN_ADDR *) server->h_addr;//inet_addr(frontEnd_sockfd) ;
	if(! endPoint.urls.empty()){
		if(IPAddress == "0.0.0.0"){
			server_address.sin_addr.s_addr =INADDR_ANY;
			portno=5001;
		}else{
			server = gethostbyname(IPAddress.c_str());
			if (server == NULL) {
				cerr<<"ERROR, no such host: "<<IPAddress<<endl;
				return false;
			}
			bcopy((char *)server->h_addr,
					(char *)&server_address.sin_addr.s_addr,
					server->h_length);
		}
	}else{
		server_address.sin_addr.s_addr =INADDR_ANY;
		portno=7001;//+(rand()%100);
	}

	bool isPortOk=false;
	this->m_mutex.lock();

	while(!isPortOk){
		server_address.sin_port = htons(portno);
		server_len = sizeof(server_address);
		int res= bind(frontEnd_sockfd, (struct sockaddr *)&server_address, server_len);

		if(res<0){
			if(errno == EADDRINUSE) {
				portno++;
			} else {
				printf("could not bind to process (%d) %s\n", errno, strerror(errno));
				break;
			}
		}else{
			isPortOk=true;
		}
	}
	this->m_mutex.unlock();

	if(isPortOk){
		stringstream ss;
		vector<string> ips;
		this->getHostIpAddresses(ips);
		for(auto s: ips){
			if((s.find("127.0.0.") != 0) && s!="0.0.0.0" && s!="::1")
				ss<<s<<";";
		}
		ss<<portno<<";"<<rank;
		endPoint.rawEndpoints = ss.str();
		endPoint.parse();
		listen(frontEnd_sockfd, 200);
	}
	return frontEnd_sockfd;
}

string Networking::openfrontendMPI(MPI_Comm globalCommunicator){

	string endpointSocketMPI;
	char port_name[MPI_MAX_PORT_NAME];
    MPI_Comm_set_errhandler(globalCommunicator, MPI_ERRORS_RETURN);
	int res=MPI_ERR_UNKNOWN;
	try{
        res=  MPI_Open_port(MPI_INFO_NULL, port_name); //<------
	}catch(const std::exception& e) {
		cout<<e.what()<<endl;
	}
    MPI_Comm_set_errhandler(globalCommunicator, MPI_ERRORS_ARE_FATAL);
    if(res== MPI_SUCCESS){
		endpointSocketMPI=port_name;
	}
	return endpointSocketMPI;
}

bool Networking::canOpenMPIPort(MPI_Comm globalCommunicator){
	bool canOpen=true;
	string endpoint = Networking::getInstance().openfrontendMPI(globalCommunicator);
	// an endpoint must contains ;
	std::size_t found=endpoint.find(";");
	if(found!=std::string::npos){
		char port_name[MPI_MAX_PORT_NAME];
		strncpy (port_name, endpoint.c_str(), endpoint.length()+1);
		MPI_Close_port(port_name);
	}else{
		canOpen=false;
	}
	return canOpen;
}

void Networking::setnonblocking(int sock)
{
	int opts;

	opts = fcntl(sock,F_GETFL);
	if (opts < 0) {
		perror("fcntl(F_GETFL)");
		exit(EXIT_FAILURE);
	}
	opts = (opts | O_NONBLOCK);
	if (fcntl(sock,F_SETFL,opts) < 0) {
		perror("fcntl(F_SETFL)");
		exit(EXIT_FAILURE);
	}
	return;
}

void Networking::setBlocking(int sock)
{
	int opts;

	opts = fcntl(sock,F_GETFL);
	if (opts < 0) {
		perror("fcntl(F_GETFL)");
		exit(EXIT_FAILURE);
	}
	opts &= (~O_NONBLOCK);
	if( fcntl(sock, F_SETFL, opts) < 0) {
		fprintf(stderr, "Error fcntl(..., F_SETFL) (%s)\n", strerror(errno));
		exit(0);
	}
	return;
}

int Networking::sock_read(int sock, void *buffer, int count, int header){

	register char *ptr = (char *) buffer;
	register int bytesleft = count;
	int chunk = this->chunk_size;
	do
	{
		register int rc;
		if(bytesleft < chunk)
			chunk=bytesleft;
		do{
			rc = read(sock, ptr, chunk);
		}while (rc < 0 && errno == EINTR);

		if (rc == 0)
			return count - bytesleft;
		else if (rc < 0)
			return -1;

		bytesleft -= rc;
		ptr += rc;
	}
	while (bytesleft);

	return count;
}

int Networking::sock_write(int sock, const void *buffer, int count, int header){

	register const char *ptr = (const char *) buffer;
	register int bytesleft = count;
	int chunk = this->chunk_size;

	do
	{
		register int rc;
		if(bytesleft < chunk)
			chunk=bytesleft;

		do{
			rc= write(sock, ptr, chunk);
		}while (rc < 0 && errno == EINTR);

		if (rc < 0)
			return -1;

		bytesleft -= rc;
		ptr += rc;
	}
	while (bytesleft);

	return count;
}

/*int Networking::sock_read(int sock, void *buffer, int count, int header){

  register char *ptr = (char *) buffer;
  register int bytesleft = count;
//int chunk = this->chunk_size;
int headersize=0;
if(header>=0){
headersize=sizeof(header);
}
//chunk-=headersize;
int packSize;
do
{
register int rc;
//if(bytesleft < chunk)
//    chunk=bytesleft;
do{
rc=read(sock, &header, headersize);
rc=read(sock, &packSize, sizeof(packSize));
rc = read(sock, ptr, packSize);
}while (rc < 0 && errno == EINTR);

if (rc == 0)
return count - bytesleft;
else if (rc < 0)
return -1;

bytesleft -= rc;
ptr += rc;
}
while (bytesleft);

return count;
}

int Networking::sock_write(int sock, const void *buffer, int count, int header){

register const char *ptr = (const char *) buffer;
register int bytesleft = count;
int chunk = this->chunk_size;
int headersize=0;
if(header>=0){
headersize=sizeof(header);
}
chunk = chunk-( headersize + sizeof(chunk));

do
{
register int rc;
if(bytesleft < chunk)
chunk=bytesleft;

do{
rc= write(sock, &header, headersize);
rc= write(sock, &chunk, sizeof(chunk));
rc= write(sock, ptr, chunk);
}while (rc < 0 && errno == EINTR);

if (rc < 0)
return -1;

bytesleft -= rc;
ptr += rc;
}
while (bytesleft);

return count;
}*/




/*********************************** Messaging *******************************************/
Messaging::~Messaging(){}

Messaging::Messaging(){}


int Messaging::sendMessageHeader(shared_ptr<MpiManager> &mpiManager, string senderName, int senderCoreId, string receiverName,
		int receiverCoreId, int operation, int dataType, int MPI_tag, int endMessage, int dest){

    // int relayerGlobalTag=mpiManager->getRank();
    // mpiManager->send<int>(&relayerGlobalTag, 1, dest, MPI_RELAYER_GLOBAL_TAG);

    int relayerGlobalTag= MPI_RELAYER_GLOBAL_TAG ;
	int len=senderName.size()+1;

	// send senderName
	mpiManager->send<int>(&len, 1, dest, relayerGlobalTag);
	mpiManager->send<char>((char*) senderName.c_str(), len, dest, relayerGlobalTag);

	// send receiverName
	len=receiverName.size()+1;
	mpiManager->send<int>(&len,1, dest, relayerGlobalTag);
	mpiManager->send<char>((char*) receiverName.c_str(), len, dest, relayerGlobalTag);

	// send senderCoreId + receiverCoreId + operation
	mpiManager->send<int>(&senderCoreId, 1, dest, relayerGlobalTag);
	mpiManager->send<int>(&receiverCoreId, 1, dest,  relayerGlobalTag);
	mpiManager->send<int>(&operation, 1, dest, relayerGlobalTag);
	mpiManager->send<int>(&dataType, 1, dest, relayerGlobalTag);
	mpiManager->send<int>(&MPI_tag, 1, dest, relayerGlobalTag);
	mpiManager->send<int>(&endMessage, 1, dest, relayerGlobalTag);
	return relayerGlobalTag;
}


void Messaging::receiveMessageHeader(shared_ptr<MpiManager> &mpiManager, Message &m, int sender, int relayerGlobalTag){
	int len;
	vector<char> buffer;
	// receive
	mpiManager->receive<int>(&len, 1, sender, relayerGlobalTag);
	buffer.resize(len);
	mpiManager->receive<char>(buffer.data(), len, sender, relayerGlobalTag);
	m.senderName=buffer.data();
	// receive
	mpiManager->receive<int>(&len,1, sender, relayerGlobalTag);
	buffer.resize(len);
	mpiManager->receive<char>(buffer.data(), len, sender, relayerGlobalTag);
	m.receiverName=buffer.data();
	vector<char>().swap(buffer);

	mpiManager->receive<int>(&m.senderCoreId, 1, sender, relayerGlobalTag);
	mpiManager->receive<int>(&m.receiverCoreId, 1, sender,  relayerGlobalTag);
	mpiManager->receive<int>(&m.operation, 1, sender, relayerGlobalTag);
	mpiManager->receive<int>(&m.datatype, 1, sender, relayerGlobalTag);
	mpiManager->receive<int>(&m.MPI_tag, 1, sender, relayerGlobalTag);
	mpiManager->receive<int>(&m.endMessage, 1, sender, relayerGlobalTag);
}

void Messaging::send(shared_ptr<MpiManager> &mpiManager, string senderName, int senderCoreId, string receiverName, int receiverCoreId,
		int operation, int dataType, int MPI_tag, int endMessage, char *data, int count, int dest){

	int relayerGlobalTag= sendMessageHeader(mpiManager, senderName, senderCoreId, receiverName, receiverCoreId, operation, dataType,
			MPI_tag, endMessage, dest);
	// iSend Data
	mpiManager->send<int>(&count, 1, dest, relayerGlobalTag);
	MPI_Request request;
	MPI_Status status;
	mpiManager->iSend<char>(data, count, dest, &request, relayerGlobalTag);
	MPI_Wait (&request, &status);
}

void Messaging::iSend(shared_ptr<MpiManager> &mpiManager, string senderName, int senderCoreId, string receiverName, int receiverCoreId,
		int operation, int dataType, int MPI_tag, int endMessage, char *data, int count, int dest, MPI_Request *request){

	int relayerGlobalTag= sendMessageHeader(mpiManager, senderName, senderCoreId, receiverName, receiverCoreId, operation, dataType,
			MPI_tag, endMessage, dest);
	// iSend Data
	mpiManager->send<int>(&count, 1, dest, relayerGlobalTag);
	mpiManager->iSend<char>(data, count, dest, request, relayerGlobalTag);
}

void Messaging::receive(shared_ptr<MpiManager> &mpiManager, Message &m, int sender, int relayerGlobalTag){

	receiveMessageHeader (mpiManager, m, sender, relayerGlobalTag);
	int len;
	mpiManager->receive<int>(&len, 1, sender, relayerGlobalTag);
	m.data.resize(len);
	MPI_Request request;
	mpiManager->iRecv<char>(m.data.data(), len, sender, &request, relayerGlobalTag);
	MPI_Status status;
	MPI_Wait (&request, &status);
}

void Messaging::iRecv(shared_ptr<MpiManager> & mpiManager, Message &m, int sender, int relayerGlobalTag, MPI_Request *request){

	receiveMessageHeader (mpiManager, m, sender, relayerGlobalTag);
	//iRecv data
	int len;
	mpiManager->receive<int>(&len, 1, sender, relayerGlobalTag);
	m.data.resize(len);
	mpiManager->iRecv<char>(m.data.data(), len, sender, request, relayerGlobalTag);
}


void Messaging::send(shared_ptr<MpiManager> &mpiManager, Message &m, int dest){
	send(mpiManager, m.senderName, m.senderCoreId, m.receiverName, m.receiverCoreId, m.operation, m.datatype, m.MPI_tag, m.endMessage, m.data.data(), m.data.size(), dest);
}


SocketURL::SocketURL(string socketurl){
    if(! socketurl.empty()){
        std::vector<std::string> elems;
        std::size_t found=socketurl.find(";");
        char delim=';';
        if(found == std::string::npos){
            delim=':';
        }
        std::stringstream ss(socketurl);
        std::string item;
        while (std::getline(ss, item, delim)) {
            elems.push_back(item);
        }
        assert(elems.size()==2);
        this->ipAddress = elems.at(0);
        std::string::size_type sz;   // alias of size_t
        port = std::stoi (elems.at(elems.size()-1),&sz);
    }
}


Endpoint::Endpoint(): grank(0){
}


Endpoint::Endpoint(string rawEndpoints): rawEndpoints(rawEndpoints){
    parse(rawEndpoints);
}

void Endpoint::parse(string rawEndpoints){

    if(rawEndpoints.empty()) return;
    this->rawEndpoints=rawEndpoints;
    urls.clear();
    std::vector<std::string> elems;
    std::size_t found=this->rawEndpoints.find(";");
    char delim=';';
    if(found == std::string::npos){
        delim=':';
    }
    std::stringstream ss(rawEndpoints);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    if(delim == ':'){
        assert(elems.size()==2);
        std::string::size_type sz;   // alias of size_t
        port = std::stoi (elems.at(1),&sz);
        SocketURL sk(elems.at(0), port);
        urls.push_back(sk);
    }else{
        assert(elems.size()>=3);
        std::string::size_type sz;   // alias of size_t
        port = std::stoi (elems.at(elems.size()-2),&sz);
        grank = std::stoi (elems.at(elems.size()-1),&sz);
        // try to parse all the ipaddresses and connect
        for (size_t i  =0; i< elems.size()-2; i++ ){
            SocketURL sk(elems.at(i), port);
            urls.push_back(sk);
        }
    }

}

void Endpoint::parse(){
    parse( rawEndpoints);
}


void Message::toBytes(vector<char> &buffer){

    int totalLength =  sizeof(int)+ this->senderName.length()+1
            +sizeof(int)+this->receiverName.size()+1
            +sizeof(this->senderCoreId)
            +sizeof(this->receiverCoreId)
            +sizeof(this->operation)
            +sizeof(this->datatype)
            +sizeof(this->MPI_tag)
            +sizeof(this->endMessage)
            +sizeof(int)+this->data.size();

    buffer.resize(totalLength);

    int iData=0;
    int ss= this->senderName.length()+1;
    memcpy((char*) &buffer[iData], (const char*)&ss, sizeof(int));
    iData+=sizeof(int);
    memcpy((char*) &buffer[iData], (const char*)this->senderName.c_str(), ss);
    iData+=ss;

    ss= this->receiverName.length()+1;
    memcpy((char*) &buffer[iData], (const char*)&ss, sizeof(int));
    iData+=sizeof(int);
    memcpy((char*) &buffer[iData], (const char*)this->receiverName.c_str(), ss);
    iData+=ss;

    memcpy((char*) &buffer[iData], (const char*)&senderCoreId, sizeof(int));
    iData+=sizeof(int);
    memcpy((char*) &buffer[iData], (const char*)&receiverCoreId, sizeof(int));
    iData+=sizeof(int);
    memcpy((char*) &buffer[iData], (const char*)&operation, sizeof(int));
    iData+=sizeof(int);
    memcpy((char*) &buffer[iData], (const char*)&datatype, sizeof(int));
    iData+=sizeof(int);
    memcpy((char*) &buffer[iData], (const char*)&MPI_tag, sizeof(int));
    iData+=sizeof(int);
    memcpy((char*) &buffer[iData], (const char*)&endMessage, sizeof(int));
    iData+=sizeof(int);
    ss= this->data.size();
    memcpy((char*) &buffer[iData], (const char*)&ss, sizeof(int));
    iData+=sizeof(int);
    memcpy((char*) &buffer[iData], (const char*) data.data(), data.size());
}

void Message::fromByte(const vector<char> &buffer){

    int iData=0;
    int ss;
    memcpy ((char*) &ss,(const char*) &buffer[iData], sizeof(int));
    iData+=sizeof(ss);
    senderName.assign(&buffer[iData], ss-1);
    iData+=ss;

    memcpy ((char*) &ss,(const char*) &buffer[iData], sizeof(int));
    iData+=sizeof(ss);
    receiverName.assign(&buffer[iData], ss-1);
    iData+=ss;

    memcpy ((char*) &senderCoreId,(const char*) &buffer[iData], sizeof(int));
    iData+=sizeof(int);
    memcpy ((char*) &receiverCoreId,(const char*) &buffer[iData], sizeof(int));
    iData+=sizeof(int);
    memcpy ((char*) &operation,(const char*) &buffer[iData], sizeof(int));
    iData+=sizeof(int);
    memcpy ((char*) &datatype,(const char*) &buffer[iData], sizeof(int));
    iData+=sizeof(int);
    memcpy ((char*) &MPI_tag,(const char*) &buffer[iData], sizeof(int));
    iData+=sizeof(int);
    memcpy ((char*) &endMessage,(const char*) &buffer[iData], sizeof(int));
    iData+=sizeof(int);

    memcpy ((char*) &ss,(const char*) &buffer[iData], sizeof(int));
    iData+=sizeof(int);
    this->data.resize(ss);
    memcpy((char*) data.data(), (const char*) &buffer[iData], ss);


}
