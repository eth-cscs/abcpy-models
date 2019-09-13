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

#ifndef MEDIATOR_H
#define MEDIATOR_H

#include "musclehpc/parallelism/mpiManager.h"
#include "musclehpc/conduits/conduit.h"
#include <iostream>     // std::cout
#include <algorithm>    // std::find
#include <vector>       // std::vector


using namespace std;
//------------------------------------------ Mediator Pattern -------------------------------

// ---------------------------------------------------- MpiSubGroup ------------------------------------------
class MpiSubGroup
{
protected:
    int localLeaderRank, remoteLeaderRank;
    MpiManager * mpiManager;
public:
    MpiSubGroup(){}
    int getLocalLeaderRank() const{return this->localLeaderRank;}
    int getRemoteLeaderRank() const{return this->remoteLeaderRank;}

    void init(int localLeaderRank, int remoteLeaderRank,  MpiManager * mpiManager){
        this->localLeaderRank=localLeaderRank;
        this->remoteLeaderRank=remoteLeaderRank;
        this->mpiManager=mpiManager;
    }
};

// ---------------------------------------------------- EventCommon ------------------------------------------

class EventCommon
{
protected:
    MpiSubGroup * eventContext; //Context Event who fires the event

protected:
    EventCommon (MpiSubGroup * source): eventContext(source){}
    virtual ~EventCommon(){}
};

/*
 * Forward declaration of ColleagueEvent class
 */
template<typename E>
class EventColleague;


//------------------------------------------ Mediator Pattern -------------------------------
/*
 * Mediator class is a singleton class and its only one instance of it will -
 * be created for each type of colleague collaboration.
 * An Instance of this class hold the reference to all ColleagueEvent objects -
 * which fire/ receive same Event Argument Type. Since this class is not -
 * meant to used public, it is exposed only to ColleagueEvent class as -
 * a friend.
 */

class Mediator
{

protected:
    typedef vector<EventCommon*> EventList;
    //List of ColleagueEvents in this community.
    EventList colleagues;

public:
    Mediator(){}
    virtual ~Mediator(){}
    static Mediator * instance; //Singleton implementation
    /*
     * Access to the only one possible type specific object of this class.
     */
    static Mediator & GetInstance(){
        if (!instance){
            instance = new Mediator();
        }
        return *instance;
    }


    /*
     * Register a ColleagueEvent in the list of colleagues in this community.
     */
    void RegisterCollegue(EventCommon *colEvent)
    {
        colleagues.push_back(colEvent);
    }
    /*
     * Remove a ColleagueEvent from the list of colleagues in this community.
     */
    void UnregisterCollegue(EventCommon *colEvent)
    {
        typename EventList::iterator itr = find(colleagues.begin(), colleagues.end(), colEvent);
        if (itr != colleagues.end())
        {
            colleagues.erase(itr);
        }
    }

    /*
     * When a ColleagueEvent is fired, notify all other Events in -
     * this community
     */
    template<typename E>
    void FireEvent(EventColleague<E> *source, E & data)
    {
        for (unsigned int i = 0; i < colleagues.size(); i++)
        {
            //Notify all Events other than the one who fired the event.
            EventColleague<E> *dest = dynamic_cast< EventColleague<E>* > (colleagues[i]);
            if (dest     != source){

                source->handlerProc(source->eventContext, data, dest->eventContext);
            }
        }
    }

    friend class EventCommon ;
};





/*
 * The CollegueEvent template class is instantiated by means of -
 * an EventArgumentType.When an object is created, it registers itself with -
 * Mediator<> matching its EventArgAtype. Thus it becomes a member of the -
 * community of ColleagueEvents with same EventArgType.
 * When a ColleagueEvent object is fired, it is informed to the community -
 * Mediator<> and Mediator broadcast the event to all other CollegueEvents -
 * in that community.
 */
template<typename E>
class EventColleague: public EventCommon
{
private:
    typedef void (*EventHandler)(MpiSubGroup *source, E & eventArg, MpiSubGroup* context);
    EventHandler handlerProc; //Event handler function pointer

public:
    /*
     * Constructor receives the event source context and the EventHandler
     * function pointer. Also register this object with the Mediator<> -
     * to join the community.
     */
    EventColleague(MpiSubGroup *source, EventHandler eventProc) : EventCommon(source),
        handlerProc(eventProc)
    {
        //Register with mediator
        Mediator::GetInstance().RegisterCollegue((EventCommon*)this);
    }
    /*
     * Destructor - unregister the object from community.
     */
    virtual ~EventColleague()
    {
        Mediator::GetInstance().UnregisterCollegue((EventCommon*)this);
    }
    /*
     * FireEvent - Inform the Mediator<> that the event object is triggered.
     * The Mediator<> then will broadcast this to other ColleagueEvents -
     * in this community.
     */

    void FireEvent(E & data)
    {
        Mediator::GetInstance().FireEvent<E>(this, data);
    }
    friend class Mediator;
};


//Define the static member of Mediator<> class
Mediator* Mediator::instance = 0;


//---

class MpiKernelWarapper: public MpiSubGroup {
protected:
    int color; //MPI color
    MpiKernel * kernel;

public:
    virtual ~MpiKernelWarapper(){}
    MpiKernelWarapper(int color, MpiKernel * kernel): color(color), kernel(kernel){
    }
    template <typename E>
    void sendData(E & data)
    {
        //create the event
        EventColleague<E> event(this, onEventFunc);
        //Fire sending data
        event.FireEvent(data);
    }

    // Handle an event notification..
    template <typename E>
    void onEventFunc(MpiSubGroup *source, E & data, MpiSubGroup* destination)
    {

    }
};

#endif // MEDIATOR_H
