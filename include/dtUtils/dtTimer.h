// This file is part of dtCore, a C++ library for robotics software
// development.
//
// This library is commercial and cannot be redistributed, and/or modified
// WITHOUT ANY ALLOWANCE OR PERMISSION OF Hyundai Motor Company.

#ifndef __DTCORE_DTTIMER_H__
#define __DTCORE_DTTIMER_H__

namespace dtCore {

class dtTimer {
public:
    dtTimer() {}
    virtual ~dtTimer() {}

public:
    void Reset();
    void Start();
    void Stop();
    
    double GetElapsedTime_sec();
    double GetElapsedTime_msec();
    double GetElapsedTime_usec();
};

}

#endif // __DTCORE_DTTIMER_H__