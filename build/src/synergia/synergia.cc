#include "synergia/synergia.h"

int
synergia_version_major()
{
    const int version_major = 2018;
    return version_major;
}

int
synergia_version_minor()
{
    // version is a two digit number. If the leading digit is 0, it will be
    // incorrectly interpeted as octal. Oi.
    const int version_minor = 102 - 100;
    return version_minor;
}

int
synergia_version_patch()
{
    // version is a two digit number. If the leading digit is 0, it will be
    // incorrectly interpeted as octal. Oi.
    const int version_patch = 120 - 100;
    return version_patch;
}

int
synergia_version_tweak()
{
    // version is a two digit number. If the leading digit is 0, it will be
    // incorrectly interpeted as octal. Oi.
    const int version_tweak = 100 - 100;
    return version_tweak;
}
