/*
TEphra TRAsport Simulator (TETRAS)
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


#include "ConsoleInterface.hpp"


ConsoleInterface::ConsoleInterface():nFlyingParticle_(0),
nDepositedParticle_(0),
nOutOfDomainParticle_(0),
iteration_(0),
time_(0.0),
dt_(1.0),
ips_(0.0),
isWaitingForKey_(false),
barrier_(new boost::barrier(2)),
isAsserted_(false),
keyListenerThread_(boost::thread(boost::bind(&ConsoleInterface::keyListener,this))){
	mainWindow_ = initscr();
	curs_set(0);
	timeout(-1);
	noecho();
	barrier_->wait();
	refreshDisp();
}

ConsoleInterface::~ConsoleInterface() {
	delwin(mainWindow_);
	endwin();
	refresh();
}


void ConsoleInterface::refreshDisp(){
	boost::mutex::scoped_lock lock(mutex_);
	erase();
	boost::posix_time::time_duration time = boost::posix_time::seconds((int)time_);
	mvprintw(5,5,"press <q> to quit (quit with <ctrl-c> could remain the terminal in unpredictible state)");
	if(isWaitingForKey_)mvprintw(7,5,"press <%c> to continue (<q> still exits the program)",waitedKey_);
	mvprintw(9,5,"number of flying particles : %i",nFlyingParticle_);
	mvprintw(10,5,"number of deposited particles : %i",nDepositedParticle_);
	mvprintw(11,5,"number of out of domain particles : %i",nOutOfDomainParticle_);
	if(isAsserted_)mvprintw(13,5,"%i + %i + % i asserted to be equal to %i",nFlyingParticle_, nDepositedParticle_, nOutOfDomainParticle_, nAssertedPart_);
	mvprintw(15,5,"iteration : %i   -   time : %ih %im %is   -   dt : %f", iteration_, time.hours(), time.minutes(), time.seconds(), dt_);
	mvprintw(17,5,"ips : %f",ips_);
	move(0,0);
	refresh();
}

void ConsoleInterface::waitForEnd(){
	keyListenerThread_.join();
}

void ConsoleInterface::registerSimulation(LocalPlume* sim){
	simulation_ = sim;
	simulation_->spawnSimulationThread();
	ipsCounterThread_ = boost::thread(boost::bind(&ConsoleInterface::ipsCounter, this));
}

void ConsoleInterface::ipsCounter(){
	int currentIteration;
	int lastIteration = 0;
	int diff;
	while(true){
		currentIteration = simulation_->getCurrentIteration();
		diff = currentIteration - lastIteration;
		if(diff==0){
			ips_ = 0.0;
		}
		else{
			ips_ = (double)diff;
		}
		lastIteration = currentIteration;
		refreshDisp();
		sleep(1);
	}
}

// this function is used as a thread
void ConsoleInterface::keyListener(){
	// wait for getch to be blocking (timeout(-1))
	barrier_->wait();
	while(true){
		refreshDisp();
		char c = getch();
		if(c == 'q') {
			// if another thread is waitingon the interface, release it
			if(isWaitingForKey_) barrier_->wait();
			// exit the simulation
			simulation_->stopSim();
			// exit the loop and the thread
			break;
		}
		else if(c == 'p'){
			simulation_->pause();
		}
		// if another thread is waiting on the interface, release it
		if(isWaitingForKey_ && c==waitedKey_) {
			barrier_->wait();
			isWaitingForKey_ = false;
			refreshDisp();
		}
	}
}

void ConsoleInterface::waitForKey(char c){
	waitedKey_ = c;
	isWaitingForKey_ = true;
	refreshDisp();
	barrier_->wait();
}

void ConsoleInterface::setNFlyingParticle(int n){
	nFlyingParticle_ = n;
	refreshDisp();
}

void ConsoleInterface::setNDepositedParticle(int n){
	nDepositedParticle_ = n;
	refreshDisp();
}

void ConsoleInterface::setNOutOfDomainParticle(int n){
	nOutOfDomainParticle_ = n;
	refreshDisp();
}

void ConsoleInterface::setIteration(int n){
	iteration_ = n;
	refreshDisp();
}

void ConsoleInterface::setTime(double t){
	dt_ = t - time_;
	time_ = t;
	refreshDisp();
}

void ConsoleInterface::assertNParticle(int n){
	assert((nFlyingParticle_+nDepositedParticle_+nOutOfDomainParticle_) == n);
	nAssertedPart_ = n;
	isAsserted_ = true;
	refreshDisp();
}
