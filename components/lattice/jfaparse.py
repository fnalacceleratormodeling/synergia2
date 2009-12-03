# /usr/bin/env python

from __future__ import division

import sys
import re
from pyparsing import Word, alphas, ParseException, Literal, CaselessLiteral \
, Combine, Optional, nums, Or, Forward, ZeroOrMore, StringEnd, alphanums, \
restOfLine,empty,delimitedList,LineEnd,Group

import math

# Debugging flag can be set to either "debug_flag=True" or "debug_flag=False"
debug_flag = False

exprStack = []
modifierDict = {}

variables = {}
lines = {}
elements = {}
commands = {}

def pushFirst(str, loc, toks):
#    print "pf:",toks
    if hasattr(toks[0],'asList'):
        exprStack.append(toks[0].asList())
    else:
        exprStack.append(toks[0].lower())

def pushUMinus( strg, loc, toks ):
    if toks and toks[0]=='-': 
        exprStack.append( 'unary -' )

def handleVarAssignment(str,loc,toks):
#    print "hva:",toks[0],exprStack
    var = toks[0].lower()
    global exprStack
#    print "about to evaluate",exprStack
    value = evaluateStack(exprStack)
    variables[var] = value

def handleLineAssignment(str,loc,toks):
#    print "hla:",toks
    line = toks[0]
    elements = toks[1]
    lines[line] = elements

def handleModifier(str,loc,toks):
#    print "hm:",toks
    modifier = toks[0]
    if len(toks) > 1:
        value = evaluateStack(exprStack)
    else:
        value = None
    modifierDict[modifier] = value

def handleElementAssignment(str,loc,toks):
#    print "hea:",toks
    global modifierDict
    element = toks[0]
    type = toks[1]
    elements[element] = (type,modifierDict)
    modifierDict = {}

def handleCommandAssignment(str,loc,toks):
    global modifierDict
#    print "hca:",toks,modifierDict
    command = toks[0]
    commands[command] = modifierDict
    modifierDict = {}
    
# define grammar
point = Literal('.')
e = CaselessLiteral('E')
plusorminus = Literal('+') | Literal('-')
number = Word(nums) 
integer = Combine(Optional(plusorminus) + number)
floatnumber = (Combine(integer + 
                       Optional(point + Optional(number)) + 
                       Optional(e + integer)
                     ) |
            Combine(point + number) + 
                       Optional(e + integer)
                     )

ident = Word(alphas, alphanums + '_') 

plus = Literal("+")
minus = Literal("-")
mult = Literal("*")
div = Literal("/")
lpar = Literal("(").suppress()
rpar = Literal(")").suppress()
lbrack= Literal("[").suppress()
rbrack = Literal("]").suppress()
addop = plus | minus
multop = mult | div
expop = Literal("^")
assign = Literal(":=").suppress()
elementassign = Literal(":").suppress()

subscriptedident = Group(ident + lbrack + ident + rbrack)

#atom = (Optional("-") + ( pi | e | fnumber | ident + lpar + expr + rpar ).setParseAction( pushFirst ) | 
#        ( lpar + expr.suppress() + rpar )).setParseAction(pushUMinus) 

expr = Forward()
atom = (Optional("-") + (floatnumber | integer | ident + lpar + expr + rpar | subscriptedident | ident).setParseAction(pushFirst) | 
         (lpar + expr.suppress() + rpar)
       ).setParseAction(pushUMinus)
        
factor = Forward()
factor << atom + ZeroOrMore((expop + factor).setParseAction(pushFirst))
        
term = factor + ZeroOrMore((multop + factor).setParseAction(pushFirst))
expr << term + ZeroOrMore((addop + term).setParseAction(pushFirst))

modifierdef = ident + Optional(Literal("=") + expr)
modifierdef.setParseAction(handleModifier)
elementdef = ident + Optional(Literal(",") + delimitedList(modifierdef))
multipleident = Optional(Literal("-")) + Optional(integer + Literal('*')) + (ident| lpar + Group(delimitedList(ident)) + rpar)
linedef = CaselessLiteral("line").suppress() + Literal("=").suppress() + lpar + Group(delimitedList(multipleident)) + rpar

bnf = ((ident + assign + expr).setParseAction(handleVarAssignment) |
       ((ident + elementassign) + linedef).setParseAction(handleLineAssignment) |
       ((ident + elementassign) + elementdef).setParseAction(handleElementAssignment) | 
       ((ident + Optional(Literal(",").suppress() + delimitedList(modifierdef))).setParseAction(handleCommandAssignment)) |
       empty)

pattern = bnf + StringEnd()
comment = Literal('!') + restOfLine
pattern.ignore(comment)
continuation = Literal('&') + restOfLine + LineEnd()
pattern.ignore(continuation)


# map operator symbols to corresponding arithmetic operations
opn = { "+" : (lambda a, b: a + b),
        "-" : (lambda a, b: a - b),
        "*" : (lambda a, b: a * b),
        "/" : (lambda a, b: a / b),
        "^" : (lambda a, b: a ** b) }

fn  = { "sin" : math.sin,
        "cos" : math.cos,
        "tan" : math.tan,
        "abs" : abs,
        "sqrt" : math.sqrt}

# Recursive function that evaluates the stack
def evaluateStack(s):
  op = s.pop()
#  print op,type(op),type([]),type(op) == type([])
  if op == 'unary -':
      return -evaluateStack( s )
  elif type(op) == type([]):
      element = op[0]
      modifier = op[1]
      return (elements[element])[1][modifier]
  elif op in "+-*/^":
    op2 = evaluateStack(s)
    op1 = evaluateStack(s)
    return opn[op](op1, op2)
  elif op == "PI":
    return math.pi
  elif op == "E":
    return math.e
  elif op in fn:
#    print "fn[%s]" % op
    return fn[op]( evaluateStack( s ) )
  elif re.search('^[a-zA-Z][a-zA-Z0-9_]*$', op):
    name = op.lower()
    if variables.has_key(name):
      return variables[name]
    else:
#      print "warning: undefined variable:",name
      return name
  elif re.search('^[-+]?[0-9]+$', op):
    return long(op)
  else:
    return float(op)

if __name__ == "__main__":
    filename = sys.argv[1]
    combined_line = ""
    for line in open(filename,"r").readlines():
        combined_line += line
        if not re.search('^[^!]*&',line):
            try:
                L = pattern.parseString(combined_line)
            except ParseException, err:
                print err.line
                print " " * (err.column - 1) + "^"
                print err
                sys.exit(1)
            combined_line = ""

    print variables
    print lines
    print elements
    print commands
