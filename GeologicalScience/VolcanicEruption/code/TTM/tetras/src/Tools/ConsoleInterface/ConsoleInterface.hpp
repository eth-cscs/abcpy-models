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


#ifndef CONSOLEINTERFACE_H_
#define CONSOLEINTERFACE_H_

#include <curses.h>
#include <boost/thread/thread.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "../Simulation/Simulation.hpp"

/* forward declaration of simulation */
class LocalPlume;

class ConsoleInterface {
public:
	ConsoleInterface();
	virtual ~ConsoleInterface();

	void waitForKey(char c);

	void waitForEnd();

	void registerSimulation(LocalPlume* sim);

	void setNFlyingParticle(int n);
	void setNDepositedParticle(int n);
	void setNOutOfDomainParticle(int n);

	void setIteration(int n);
	void setTime(double t);

	void assertNParticle(int n);
private:
	WINDOW* mainWindow_;
	LocalPlume* simulation_;

	void refreshDisp();
	void keyListener();
	void ipsCounter();

	int nFlyingParticle_;
	int nDepositedParticle_;
	int nOutOfDomainParticle_;

	int iteration_;
	double time_;
	double dt_;
	double ips_;

	bool isWaitingForKey_;
	char waitedKey_;
	boost::barrier* barrier_;

	boost::mutex mutex_;

	int nAssertedPart_;
	bool isAsserted_;

	boost::thread keyListenerThread_;
	boost::thread ipsCounterThread_;
};

#endif /* CONSOLEINTERFACE_H_ */
