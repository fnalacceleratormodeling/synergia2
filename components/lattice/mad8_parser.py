#!/usr/bin/env python

from __future__ import division
from pyparsing import Word, alphas, ParseException, Literal, CaselessLiteral, \
Combine, Optional, nums, Or, Forward, ZeroOrMore, StringEnd, alphanums, \
restOfLine, empty, delimitedList, LineEnd, Group, QuotedString, dblQuotedString, \
removeQuotes, sglQuotedString

import sys
import re
import math
import operator

class Stack_type:
    floatnumber = 1
    ident = 2
    operator = 3
    function = 4
    unary_minus = 5

class Stack_item:
    def __init__(self, type, value):
        self.type = type
        self.value = value

    def __repr__(self):
        retval = '<'
        if self.type == Stack_type.floatnumber :
            retval += 'floatnumber'
        elif self.type == Stack_type.ident:
            retval += 'ident'
        elif self.type == Stack_type.operator:
            retval += 'operator'
        elif self.type == Stack_type.function:
            retval += 'function'
        elif self.type == Stack_type.ident:
            retval += 'unary_minus'
        retval += ':' + str(self.value) + '>'
        return retval

class Expression_parser:
    def __init__(self, uninitialized_warn=True,
                 uninitialized_exception=False):
        self.uninitialized_warn = uninitialized_warn
        self.uninitialized_exception = uninitialized_exception
        self.bnf = self._construct_bnf()
        self.operators = {"+" : operator.add,
                          "-" : operator.sub,
                          "*" : operator.mul,
                          "/" : operator.truediv,
                          "^" : operator.pow}
        # self.functions values are of the form [function, number of arguments]
        self.functions = {"sqrt" : [math.sqrt, 1],
                           "log" : [math.log, 1],
                           "exp" : [math.exp, 1],
                           "sin" : [math.sin, 1],
                           "cos" : [math.cos, 1],
                           "tan" : [math.tan, 1],
                           "asin" : [math.asin, 1],
                           "abs" : [abs, 1],
                           "max" : [max, 2],
                           "min" : [min, 2]}
        self.constants = {"pi" : math.pi,
                          "twopi" : 2.0 * math.pi,
                          "degrad" : 180.0 / math.pi,
                          "raddeg" : math.pi / 180.0,
                          "e" : math.e,
                          "emass" : 0.51099906e-3,
                          "pmass" : 0.93827231,
                          "clight" : 2.99792458e8}

    def _construct_bnf(self):
        point = Literal(".")
        e = CaselessLiteral('e') | CaselessLiteral('d')
        plusorminus = Literal('+') | Literal('-')
        number = Word(nums) 
        integer = Combine(Optional(plusorminus) + number)
        floatnumber = (Combine(integer + \
                               Optional(point + Optional(number)) + \
                               Optional(e + integer)) | \
                       Optional(plusorminus) + Combine(point + number) + \
                        Optional(e + integer))
        ident = Word(alphas, alphanums + '_' + '.' + "'")
     
        plus = Literal("+")
        minus = Literal("-")
        mult = Literal("*")
        div = Literal("/")
        lpar = Literal("(").suppress()
        rpar = Literal(")").suppress()
        addop = plus | minus
        multop = mult | div
        expop = Literal("^")
        
        expr = Forward()
        atom = (Optional(minus) + 
                ((floatnumber).setParseAction(self._push_floatnumber) | \
                 (ident + lpar + delimitedList(expr) + rpar).setParseAction(self._push_function) | \
                 (ident).setParseAction(self._push_ident) | \
                (lpar + expr.suppress() + rpar))).setParseAction(self._push_uminus) 
        
        # by defining exponentiation as "atom [ ^ factor ]..." instead of
        # "atom [ ^ atom ]...", we get right-to-left exponents, instead of
        # left-to-right, that is, 2^3^2 = 2^(3^2), not (2^3)^2.
        factor = Forward()
        factor << atom + ZeroOrMore((expop + factor).setParseAction(self._push_operator))
        
        term = factor + ZeroOrMore((multop + factor).setParseAction(self._push_operator))
        expr << term + ZeroOrMore((addop + term).setParseAction(self._push_operator))
        
        return expr

    def _push_floatnumber(self, strg, loc, toks):
        numstr = toks[0]
        numstr = numstr.replace('d', 'e')
        numstr = numstr.replace('D', 'e')
        self.stack.append(Stack_item(Stack_type.floatnumber, float(numstr)))
    
    def _push_ident(self, strg, loc, toks):
        if not self.functions.has_key(toks[0]) and \
            self.last_ident_loc != loc:
            self.stack.append(Stack_item(Stack_type.ident, toks[0]))
            self.last_ident_loc = loc
                          
    def _push_function(self, strg, loc, toks):
        self.stack.append(Stack_item(Stack_type.function, toks[0]))
        
    def _push_uminus(self, strg, loc, toks):
        if toks and toks[0] == '-': 
            self.stack.append(Stack_item(Stack_type.unary_minus, None))
                          
    def _push_operator(self, strg, loc, toks):
        self.stack.append(Stack_item(Stack_type.operator, toks[0]))
        
    def evaluate_stack(self, s, variables={}, constants=None):
        if constants == None:
            constants = self.constants
        op = s.pop()
        if op.type == Stack_type.unary_minus:
            return - self.evaluate_stack(s, variables, constants)
        if op.type == Stack_type.floatnumber:
            return op.value
        if op.type == Stack_type.operator:
            op2 = self.evaluate_stack(s, variables, constants)
            op1 = self.evaluate_stack(s, variables, constants)
            return self.operators[op.value](op1, op2)
        elif op.type == Stack_type.function:
            ops = []
            for i in range(0, self.functions[op.value][1]):
                ops.append(self.evaluate_stack(s, variables, constants))
            return apply(self.functions[op.value][0], ops)
        elif op.type == Stack_type.ident:
            if op.value in constants:
                return constants[op.value]
            elif op.value in variables:
                return variables[op.value]
            else:
                if self.uninitialized_exception:
                    raise ParseException, 'variable "%s" unitialized' % op.value
                if self.uninitialized_warn:
                    print \
                     'warning: variable "%s" uninitialized, treating as 0.0' \
                     % op.value
                return 0.0
        else:
            print "jfa: huh?", op
    
    def parse(self, text):
        self.stack = []
        self.last_ident_loc = - 1
        result = self.bnf.parseString(text)
        return self.stack

if __name__ == '__main__':
    ep = Expression_parser()
    stack = ep.parse('2*pi+4*5')
    print 'ans =', ep.evaluate_stack(stack)

    stack = ep.parse('sin(2.3)')
    print 'ans =', ep.evaluate_stack(stack)

    stack = ep.parse('-sin(2.3)')
    print 'ans =', ep.evaluate_stack(stack)

    stack = ep.parse('cos(x)')
    print 'ans =', ep.evaluate_stack(stack)


