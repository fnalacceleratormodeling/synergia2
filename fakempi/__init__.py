class World_dummy:
    def __init__(self):
        pass
    def Barrier(self):
        pass
    
class MPI_dummy:
    def __init__(self):
        self.rank = 0
        self.size = 1
        self.WORLD = World_dummy()
    

MPI = MPI_dummy()
