var data = {lines:[
{"lineNum":"    1","line":"/* Optimizing macros and inline functions for stdio functions."},
{"lineNum":"    2","line":"   Copyright (C) 1998-2020 Free Software Foundation, Inc."},
{"lineNum":"    3","line":"   This file is part of the GNU C Library."},
{"lineNum":"    4","line":""},
{"lineNum":"    5","line":"   The GNU C Library is free software; you can redistribute it and/or"},
{"lineNum":"    6","line":"   modify it under the terms of the GNU Lesser General Public"},
{"lineNum":"    7","line":"   License as published by the Free Software Foundation; either"},
{"lineNum":"    8","line":"   version 2.1 of the License, or (at your option) any later version."},
{"lineNum":"    9","line":""},
{"lineNum":"   10","line":"   The GNU C Library is distributed in the hope that it will be useful,"},
{"lineNum":"   11","line":"   but WITHOUT ANY WARRANTY; without even the implied warranty of"},
{"lineNum":"   12","line":"   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU"},
{"lineNum":"   13","line":"   Lesser General Public License for more details."},
{"lineNum":"   14","line":""},
{"lineNum":"   15","line":"   You should have received a copy of the GNU Lesser General Public"},
{"lineNum":"   16","line":"   License along with the GNU C Library; if not, see"},
{"lineNum":"   17","line":"   <https://www.gnu.org/licenses/>.  */"},
{"lineNum":"   18","line":""},
{"lineNum":"   19","line":"#ifndef _BITS_STDIO_H"},
{"lineNum":"   20","line":"#define _BITS_STDIO_H 1"},
{"lineNum":"   21","line":""},
{"lineNum":"   22","line":"#ifndef _STDIO_H"},
{"lineNum":"   23","line":"# error \"Never include <bits/stdio.h> directly; use <stdio.h> instead.\""},
{"lineNum":"   24","line":"#endif"},
{"lineNum":"   25","line":""},
{"lineNum":"   26","line":"#ifndef __extern_inline"},
{"lineNum":"   27","line":"# define __STDIO_INLINE inline"},
{"lineNum":"   28","line":"#else"},
{"lineNum":"   29","line":"# define __STDIO_INLINE __extern_inline"},
{"lineNum":"   30","line":"#endif"},
{"lineNum":"   31","line":""},
{"lineNum":"   32","line":""},
{"lineNum":"   33","line":"#ifdef __USE_EXTERN_INLINES"},
{"lineNum":"   34","line":"/* For -D_FORTIFY_SOURCE{,=2} bits/stdio2.h will define a different"},
{"lineNum":"   35","line":"   inline.  */"},
{"lineNum":"   36","line":"# if !(__USE_FORTIFY_LEVEL > 0 && defined __fortify_function)"},
{"lineNum":"   37","line":"/* Write formatted output to stdout from argument list ARG.  */"},
{"lineNum":"   38","line":"__STDIO_INLINE int"},
{"lineNum":"   39","line":"vprintf (const char *__restrict __fmt, __gnuc_va_list __arg)","class":"lineNoCov","hits":"0","possible_hits":"15",},
{"lineNum":"   40","line":"{"},
{"lineNum":"   41","line":"  return vfprintf (stdout, __fmt, __arg);","class":"lineNoCov","hits":"0","possible_hits":"15",},
{"lineNum":"   42","line":"}"},
{"lineNum":"   43","line":"# endif"},
{"lineNum":"   44","line":""},
{"lineNum":"   45","line":"/* Read a character from stdin.  */"},
{"lineNum":"   46","line":"__STDIO_INLINE int"},
{"lineNum":"   47","line":"getchar (void)"},
{"lineNum":"   48","line":"{"},
{"lineNum":"   49","line":"  return getc (stdin);"},
{"lineNum":"   50","line":"}"},
{"lineNum":"   51","line":""},
{"lineNum":"   52","line":""},
{"lineNum":"   53","line":"# ifdef __USE_MISC"},
{"lineNum":"   54","line":"/* Faster version when locking is not necessary.  */"},
{"lineNum":"   55","line":"__STDIO_INLINE int"},
{"lineNum":"   56","line":"fgetc_unlocked (FILE *__fp)"},
{"lineNum":"   57","line":"{"},
{"lineNum":"   58","line":"  return __getc_unlocked_body (__fp);"},
{"lineNum":"   59","line":"}"},
{"lineNum":"   60","line":"# endif /* misc */"},
{"lineNum":"   61","line":""},
{"lineNum":"   62","line":""},
{"lineNum":"   63","line":"# ifdef __USE_POSIX"},
{"lineNum":"   64","line":"/* This is defined in POSIX.1:1996.  */"},
{"lineNum":"   65","line":"__STDIO_INLINE int"},
{"lineNum":"   66","line":"getc_unlocked (FILE *__fp)"},
{"lineNum":"   67","line":"{"},
{"lineNum":"   68","line":"  return __getc_unlocked_body (__fp);"},
{"lineNum":"   69","line":"}"},
{"lineNum":"   70","line":""},
{"lineNum":"   71","line":"/* This is defined in POSIX.1:1996.  */"},
{"lineNum":"   72","line":"__STDIO_INLINE int"},
{"lineNum":"   73","line":"getchar_unlocked (void)"},
{"lineNum":"   74","line":"{"},
{"lineNum":"   75","line":"  return __getc_unlocked_body (stdin);"},
{"lineNum":"   76","line":"}"},
{"lineNum":"   77","line":"# endif\t/* POSIX */"},
{"lineNum":"   78","line":""},
{"lineNum":"   79","line":""},
{"lineNum":"   80","line":"/* Write a character to stdout.  */"},
{"lineNum":"   81","line":"__STDIO_INLINE int"},
{"lineNum":"   82","line":"putchar (int __c)"},
{"lineNum":"   83","line":"{"},
{"lineNum":"   84","line":"  return putc (__c, stdout);"},
{"lineNum":"   85","line":"}"},
{"lineNum":"   86","line":""},
{"lineNum":"   87","line":""},
{"lineNum":"   88","line":"# ifdef __USE_MISC"},
{"lineNum":"   89","line":"/* Faster version when locking is not necessary.  */"},
{"lineNum":"   90","line":"__STDIO_INLINE int"},
{"lineNum":"   91","line":"fputc_unlocked (int __c, FILE *__stream)"},
{"lineNum":"   92","line":"{"},
{"lineNum":"   93","line":"  return __putc_unlocked_body (__c, __stream);"},
{"lineNum":"   94","line":"}"},
{"lineNum":"   95","line":"# endif /* misc */"},
{"lineNum":"   96","line":""},
{"lineNum":"   97","line":""},
{"lineNum":"   98","line":"# ifdef __USE_POSIX"},
{"lineNum":"   99","line":"/* This is defined in POSIX.1:1996.  */"},
{"lineNum":"  100","line":"__STDIO_INLINE int"},
{"lineNum":"  101","line":"putc_unlocked (int __c, FILE *__stream)"},
{"lineNum":"  102","line":"{"},
{"lineNum":"  103","line":"  return __putc_unlocked_body (__c, __stream);"},
{"lineNum":"  104","line":"}"},
{"lineNum":"  105","line":""},
{"lineNum":"  106","line":"/* This is defined in POSIX.1:1996.  */"},
{"lineNum":"  107","line":"__STDIO_INLINE int"},
{"lineNum":"  108","line":"putchar_unlocked (int __c)"},
{"lineNum":"  109","line":"{"},
{"lineNum":"  110","line":"  return __putc_unlocked_body (__c, stdout);"},
{"lineNum":"  111","line":"}"},
{"lineNum":"  112","line":"# endif\t/* POSIX */"},
{"lineNum":"  113","line":""},
{"lineNum":"  114","line":""},
{"lineNum":"  115","line":"# ifdef\t__USE_GNU"},
{"lineNum":"  116","line":"/* Like `getdelim\', but reads up to a newline.  */"},
{"lineNum":"  117","line":"__STDIO_INLINE __ssize_t"},
{"lineNum":"  118","line":"getline (char **__lineptr, size_t *__n, FILE *__stream)"},
{"lineNum":"  119","line":"{"},
{"lineNum":"  120","line":"  return __getdelim (__lineptr, __n, \'\\n\', __stream);"},
{"lineNum":"  121","line":"}"},
{"lineNum":"  122","line":"# endif /* GNU */"},
{"lineNum":"  123","line":""},
{"lineNum":"  124","line":""},
{"lineNum":"  125","line":"# ifdef __USE_MISC"},
{"lineNum":"  126","line":"/* Faster versions when locking is not required.  */"},
{"lineNum":"  127","line":"__STDIO_INLINE int"},
{"lineNum":"  128","line":"__NTH (feof_unlocked (FILE *__stream))"},
{"lineNum":"  129","line":"{"},
{"lineNum":"  130","line":"  return __feof_unlocked_body (__stream);"},
{"lineNum":"  131","line":"}"},
{"lineNum":"  132","line":""},
{"lineNum":"  133","line":"/* Faster versions when locking is not required.  */"},
{"lineNum":"  134","line":"__STDIO_INLINE int"},
{"lineNum":"  135","line":"__NTH (ferror_unlocked (FILE *__stream))"},
{"lineNum":"  136","line":"{"},
{"lineNum":"  137","line":"  return __ferror_unlocked_body (__stream);"},
{"lineNum":"  138","line":"}"},
{"lineNum":"  139","line":"# endif /* misc */"},
{"lineNum":"  140","line":""},
{"lineNum":"  141","line":"#endif /* Use extern inlines.  */"},
{"lineNum":"  142","line":""},
{"lineNum":"  143","line":""},
{"lineNum":"  144","line":"#if defined __USE_MISC && defined __GNUC__ && defined __OPTIMIZE__ \\"},
{"lineNum":"  145","line":"    && !defined __cplusplus"},
{"lineNum":"  146","line":"/* Perform some simple optimizations.  */"},
{"lineNum":"  147","line":"# define fread_unlocked(ptr, size, n, stream) \\"},
{"lineNum":"  148","line":"  (__extension__ ((__builtin_constant_p (size) && __builtin_constant_p (n)    \\"},
{"lineNum":"  149","line":"\t\t   && (size_t) (size) * (size_t) (n) <= 8\t\t      \\"},
{"lineNum":"  150","line":"\t\t   && (size_t) (size) != 0)\t\t\t\t      \\"},
{"lineNum":"  151","line":"\t\t  ? ({ char *__ptr = (char *) (ptr);\t\t\t      \\"},
{"lineNum":"  152","line":"\t\t       FILE *__stream = (stream);\t\t\t      \\"},
{"lineNum":"  153","line":"\t\t       size_t __cnt;\t\t\t\t\t      \\"},
{"lineNum":"  154","line":"\t\t       for (__cnt = (size_t) (size) * (size_t) (n);\t      \\"},
{"lineNum":"  155","line":"\t\t\t    __cnt > 0; --__cnt)\t\t\t\t      \\"},
{"lineNum":"  156","line":"\t\t\t {\t\t\t\t\t\t      \\"},
{"lineNum":"  157","line":"\t\t\t   int __c = getc_unlocked (__stream);\t\t      \\"},
{"lineNum":"  158","line":"\t\t\t   if (__c == EOF)\t\t\t\t      \\"},
{"lineNum":"  159","line":"\t\t\t     break;\t\t\t\t\t      \\"},
{"lineNum":"  160","line":"\t\t\t   *__ptr++ = __c;\t\t\t\t      \\"},
{"lineNum":"  161","line":"\t\t\t }\t\t\t\t\t\t      \\"},
{"lineNum":"  162","line":"\t\t       ((size_t) (size) * (size_t) (n) - __cnt)\t\t      \\"},
{"lineNum":"  163","line":"\t\t\t/ (size_t) (size); })\t\t\t\t      \\"},
{"lineNum":"  164","line":"\t\t  : (((__builtin_constant_p (size) && (size_t) (size) == 0)   \\"},
{"lineNum":"  165","line":"\t\t      || (__builtin_constant_p (n) && (size_t) (n) == 0))     \\"},
{"lineNum":"  166","line":"\t\t\t/* Evaluate all parameters once.  */\t\t      \\"},
{"lineNum":"  167","line":"\t\t     ? ((void) (ptr), (void) (stream), (void) (size),\t      \\"},
{"lineNum":"  168","line":"\t\t\t(void) (n), (size_t) 0)\t\t\t\t      \\"},
{"lineNum":"  169","line":"\t\t     : fread_unlocked (ptr, size, n, stream))))"},
{"lineNum":"  170","line":""},
{"lineNum":"  171","line":"# define fwrite_unlocked(ptr, size, n, stream) \\"},
{"lineNum":"  172","line":"  (__extension__ ((__builtin_constant_p (size) && __builtin_constant_p (n)    \\"},
{"lineNum":"  173","line":"\t\t   && (size_t) (size) * (size_t) (n) <= 8\t\t      \\"},
{"lineNum":"  174","line":"\t\t   && (size_t) (size) != 0)\t\t\t\t      \\"},
{"lineNum":"  175","line":"\t\t  ? ({ const char *__ptr = (const char *) (ptr);\t      \\"},
{"lineNum":"  176","line":"\t\t       FILE *__stream = (stream);\t\t\t      \\"},
{"lineNum":"  177","line":"\t\t       size_t __cnt;\t\t\t\t\t      \\"},
{"lineNum":"  178","line":"\t\t       for (__cnt = (size_t) (size) * (size_t) (n);\t      \\"},
{"lineNum":"  179","line":"\t\t\t    __cnt > 0; --__cnt)\t\t\t\t      \\"},
{"lineNum":"  180","line":"\t\t\t if (putc_unlocked (*__ptr++, __stream) == EOF)\t      \\"},
{"lineNum":"  181","line":"\t\t\t   break;\t\t\t\t\t      \\"},
{"lineNum":"  182","line":"\t\t       ((size_t) (size) * (size_t) (n) - __cnt)\t\t      \\"},
{"lineNum":"  183","line":"\t\t\t/ (size_t) (size); })\t\t\t\t      \\"},
{"lineNum":"  184","line":"\t\t  : (((__builtin_constant_p (size) && (size_t) (size) == 0)   \\"},
{"lineNum":"  185","line":"\t\t      || (__builtin_constant_p (n) && (size_t) (n) == 0))     \\"},
{"lineNum":"  186","line":"\t\t\t/* Evaluate all parameters once.  */\t\t      \\"},
{"lineNum":"  187","line":"\t\t     ? ((void) (ptr), (void) (stream), (void) (size),\t      \\"},
{"lineNum":"  188","line":"\t\t\t(void) (n), (size_t) 0)\t\t\t\t      \\"},
{"lineNum":"  189","line":"\t\t     : fwrite_unlocked (ptr, size, n, stream))))"},
{"lineNum":"  190","line":"#endif"},
{"lineNum":"  191","line":""},
{"lineNum":"  192","line":"/* Define helper macro.  */"},
{"lineNum":"  193","line":"#undef __STDIO_INLINE"},
{"lineNum":"  194","line":""},
{"lineNum":"  195","line":"#endif /* bits/stdio.h.  */"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-14 19:20:37", "instrumented" : 2, "covered" : 0,};
var merged_data = [];
