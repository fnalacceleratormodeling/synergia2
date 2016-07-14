#ifndef RESTRICT_EXTENSION_H_

#ifdef NO_RESTRICT_EXTENSION
    #define RESTRICT
#else
    #define RESTRICT __restrict__
#endif

#endif // RESTRICT_EXTENSION_H_
