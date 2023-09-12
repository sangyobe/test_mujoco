// This file is part of dtCore, a C++ library for robotics software
// development.
//
// This library is commercial and cannot be redistributed, and/or modified
// WITHOUT ANY ALLOWANCE OR PERMISSION OF Hyundai Motor Company.

#ifndef __DTCORE_DTCONTROLLER_H__
#define __DTCORE_DTCONTROLLER_H__

namespace dtCore {

template<typename RealType>
class dtController {
public:
    dtController() {}
    virtual ~dtController() {}

public:
    virtual void ReadDevices() = 0;
    virtual void Estimate() = 0;
    virtual void Reflect() = 0;
    virtual void Compute() = 0;
    virtual void WriteDevices() = 0;

};

}

#endif // __DTCORE_DTCONTROLLER_H__