#ifndef FF_MISC_ELEMENTS_H
#define FF_MISC_ELEMENTS_H

#include "ff_element.h"
#include "ff_drift.h"

// =========== monitor =================
class FF_monitor : public FF_drift
{
public:
    FF_monitor() {}
    virtual ~FF_monitor() {}
};

typedef boost::shared_ptr<FF_monitor > FF_monitor_sptr;

// =========== hmonitor =================
class FF_hmonitor : public FF_drift
{
public:
    FF_hmonitor() { }
    virtual ~FF_hmonitor() {}
};

typedef boost::shared_ptr<FF_hmonitor > FF_hmonitor_sptr;

// =========== vmonitor =================

class FF_vmonitor : public FF_drift
{
public:
    FF_vmonitor() { }
    virtual ~FF_vmonitor() {}
};

typedef boost::shared_ptr<FF_vmonitor > FF_vmonitor_sptr;

// =========== instrument =================

class FF_instrument : public FF_drift
{
public:
    FF_instrument() { }
    virtual ~FF_instrument() {}
};

typedef boost::shared_ptr<FF_instrument > FF_instrument_sptr;

#endif // FF_MISC_ELEMENTS_H
