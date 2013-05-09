#!/usr/bin/env python


# Sphinx latex output contains an excessive amount of whitespace.
# This script removes empty quote environments in order to 
# reduce the worst whitespace excesses.

import re, string
import sys, os

def line_join(words):
    retval = ''
    line_words = []
    for word in words:
        line_words.append(word)
        if (len(word) > 0) and (word[-1] == '\n'):
            retval += string.join(line_words)
            line_words = []
    if len(line_words) > 0:
        retval += string.join(line_words)
    return retval

f = open(sys.argv[1], 'r')
tmpname = 'remove-empty-quotes-tmp.tex'
o = open(tmpname, 'w')
quote_level = 0
quote = []
bquote = '\\begin{quote}'
equote = '\\end{quote}'
bquoten = '\\begin{quote}\n'
equoten = '\\end{quote}\n'
for line in f.readlines():
    out = []
    for word in line.split():
        length = len(word)
        if word == bquote:
            quote_level += 1
            quote.append(bquoten)
        elif (length > 13) and \
            (word[(length-13):length] == bquote):
            if quote_level == 0:
                out.append(word[:(length-13)])
                o.write(string.join(out))
                out = []
            else:
                quote.append(word[:(length-13)])
            quote_level += 1
            quote.append(bquoten)
        elif word == equote:
            quote_level -= 1
            if quote[-1] == bquoten:
                quote.pop()
            else:
                quote.append(equoten)
            if quote_level == 0:
                if len(quote) > 1:
                    o.write(line_join(quote))
                quote = []
        else:
            if quote_level > 0:
                quote.append(word)
            else:
                out.append(word)
    if (len(out) > 0) or \
        ((len(line.split()) == 0) and (quote_level == 0)):
        o.write(string.join(out)+'\n')
    if (quote_level > 0) and \
        (quote[-1] != bquoten) and (quote[-1] !=equoten):
        quote.append('\n') 
        
f.close()
o.close()
os.rename(tmpname, sys.argv[1])    