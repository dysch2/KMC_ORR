
class Move(object):
    '''Data structure to store the kinematics of particle movement'''
    def __init__(self,rate,in_state,fi_state,react,position,gcn):
        self.rate = rate
        self.in_state = in_state
        self.fi_state = fi_state
        self.react = react
        self.position = position
        self.gcn = gcn


class reaction(object): 
# Eact and Ereact are activation and reaction energies respectively, ist and fst the initial and final states of the reactants
    def __init__(self,Eact,Ereact,ist,fst,ffreq,bfreq):
        self.freact = Eact
        self.breact = Ereact - Eact
        self.istate = ist
        self.fstate = fst
        self.ffreq = ffreq
        self.bfreq = bfreq

class occ(object):
    def __init__(self,a1,s1,a2,s2,rt,pt,b):
        self.ego = ego # assigned id number to allow us to properly account for the particle
        self.a1 = a1 # gcn type of atom 1
        self.s1 = s1 # position of atom 1
        self.a2 = a2 # gcn type of atom 2
        self.s2 = s2 # position of atom 2
        if self.a1 == 'l' and self.a2 =='l':
            self.rt = 'l' # reaction types for reactant (e.g l,m,h)
            self.pt = 'll' # pair type (e.g ll,lm,mm,mh,hh)
        if self.a1 == 'h' and self.a2 =='h':
            self.rt = 'h' # reaction types for reactant (e.g l,m,h)
            self.pt = 'hh' # pair type (e.g ll,lm,mm,mh,hh)
        if (self.a1 == 'l' and self.a2 =='m') or (self.a1 == 'm' and self.a2 == 'l'):
            self.rt = 'm' # reaction types for reactant (e.g l,m,h)
            self.pt = 'lm' # pair type (e.g ll,lm,mm,mh,hh)
        if (self.a1 == 'h' and self.a2 =='m') or (self.a1 == 'm' and self.a2 == 'h'):
            self.rt = 'm' # reaction types for reactant (e.g l,m,h)
            self.pt = 'mh' # pair type (e.g ll,lm,mm,mh,hh)
        if self.a1 == 'm' and self.a2 =='m':
            self.rt = 'm' # reaction types for reactant (e.g l,m,h)
            self.pt = 'mm' # pair type (e.g ll,lm,mm,mh,hh)
        self.b = b # bonded (1 = yes, 0 = no)

#class erg(object):
#    def __init__(self):
