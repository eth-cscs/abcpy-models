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

#ifndef EVENTLOGGER_H_
#define EVENTLOGGER_H_

#include <chrono>
#include <sys/stat.h>
#include <unordered_map>
#include <vector>
#include <fstream>

#include <iostream>

// from https://stackoverflow.com/questions/25225948/how-to-check-if-a-file-exists-in-c-with-fstreamopen?rq=1
inline bool fileExists (const std::string& filename) {
  struct stat buffer;   
  return (stat (filename.c_str(), &buffer) == 0); 
}

struct Event{
    Event(){
        timeStamp_ = 0;
        name_ = "bidon";
    }
    Event(unsigned long long ts, std::string n): timeStamp_(ts), name_(n){}
    unsigned long long timeStamp_;
    std::string name_;
};

class EventLogger {
public :

    EventLogger(std::string fn, unsigned long long toAccumulate): fileName_(fn), eventsToAccumulate_(toAccumulate) {
        #ifdef LOG
        if(fileExists(fileName_)){
            std::throw_with_nested( std::runtime_error("EventLogger : file " + fileName_ + " already exists") );
        }
        std::ofstream writer;
        writer.open(fileName_);
        writer << "time,event,eventcounter" << std::endl;
        writer.close();
        events_ = std::vector<Event>(eventsToAccumulate_);
        eventIndex_ = 0;
        #endif
    }

    ~EventLogger(){ 
        #ifdef LOG
        flush(); 
        #endif
    }

    inline void flush(){
        #ifdef LOG
        std::ofstream writer;
        writer.open(fileName_, std::ios::app);
        for(unsigned long long i=0; i<eventIndex_; i++){
            auto e = events_[i];
            if(eventsCounter_.count(e.name_) == 0) eventsCounter_[e.name_] = 0;
            else eventsCounter_[e.name_]++; 
            writer << e.timeStamp_ << "," << e.name_ << "," << eventsCounter_[e.name_] << std::endl;
        }
        events_ = std::vector<Event>(eventsToAccumulate_);
        writer.close();
        eventIndex_ = 0;
        #endif
    }

    inline void log(const char* name){
        #ifdef LOG
        // time from https://stackoverflow.com/questions/19555121/how-to-get-current-timestamp-in-milliseconds-since-1970-just-the-way-java-gets
        auto t = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
        events_[eventIndex_].name_ = std::string(name);
        events_[eventIndex_].timeStamp_ = t;

        eventIndex_++;
        if( eventIndex_ == eventsToAccumulate_ ) flush();
        #endif
    }

protected :

    std::string fileName_;
    unsigned long long eventsToAccumulate_;
    std::vector<Event> events_;
    std::unordered_map<std::string, unsigned long long> eventsCounter_;
    unsigned long long eventIndex_;
    struct timespec tspec_;

};  

#endif