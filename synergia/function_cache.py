#!/usr/bin/env python

import pickle

class Function_cache:
    def __init__(self,file_name):
        self.file_name = file_name
        try:
            cache_file = open(self.file_name,"r")
            self.cache = pickle.Unpickler(cache_file).load()
        except:
            self.cache = {}
    def in_cache(self, *args):
        return self.cache.has_key(args)
    def add(self, value, *args):
        self.cache[args] = value
        file = open(self.file_name,"w")
        pickle.Pickler(file).dump(self.cache)
        file.close()
    def get(self, *args):
        return self.cache[args]
    def list(self):
        for key in self.cache.keys():
            print "key =",key
            
