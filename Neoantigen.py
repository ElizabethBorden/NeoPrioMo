'''
@author: Marta Luksza, mluksza@ias.edu
'''

class Neoantigen(object):
    '''
    classdocs
    '''
    L=1. #default concentration of peptides
    M=1. #default concentration of mutant peptides
    W=1. #default concentration of wildtype peptides
    
    WEPS=0.0003
    WTCAP=float("inf")

    CONSORTIUM_HYDROPHOBIC_RESIDUES="VILFMWC"
    HYDROPHOBIC_RESIDUES="AILMFWYV"
    WEIRD_RESIDUES="CGP"
    
    AGGR="MAX"
    KDnormalize=1
    AGGRnum=float("inf")

    @staticmethod
    def residueChangeClass(res1,res2):
        code=""
        if res1 in Neoantigen.HYDROPHOBIC_RESIDUES:
            code+="H"
        elif res2 in Neoantigen.WEIRD_RESIDUES:
            code+="W"
        else:
            code+="N"
        if res2 in Neoantigen.HYDROPHOBIC_RESIDUES:
            code+="H"
        elif res2 in Neoantigen.WEIRD_RESIDUES:
            code+="W"
        else:
            code+="N"
        return code

    def __init__(self, params):
        '''
        Constructor
        '''
        pparams=params
        if len(params)==9:
            pparams.append("1")
        #[nid,mid,sample,wtPeptide,mtPeptide,allele,wtScore,mtScore,HLA,chopscore]=params
        [nid,mid,sample,wtPeptide,mtPeptide,allele,wtScore,mtScore]=params
        self.id=int(nid)
        self.mid=mid
        self.sample=sample
        self.wtPeptide=wtPeptide
        self.mtPeptide=mtPeptide
        self.allele=allele
        self.potential=-1e10
        try:
            self.wtkD=min(Neoantigen.WTCAP,float(wtScore))            
            self.kD=float(mtScore)
            self.setA()
        except:
            self.kD=Neoantigen.INF
            self.wtkD=Neoantigen.INF
            self.A=1.

    def getSampleName(self):
        return self.sample

    def getHydro(self):
        ''' 
        Returns hydrophobicity fraction of the mutated peptide
        '''
        res = list(self.mtPeptide)
        w=0
        for i in range(0, len(self.mtPeptide)):
            if res[i] in Neoantigen.HYDROPHOBIC_RESIDUES:
                w=w+1
        w=w/len(self.mtPeptide)
        return w

    def getConsortiumHydro(self):
        '''
        Returns hydrophobicity fraction of the mutated peptide using the Consortium's definition of the mutated peptide (Wells 2020 https://doi.org/10.1016/j.cell.2020.09.015)
        '''
        res = list(self.mtPeptide)
        wc=0
        for i in range(0, len(self.mtPeptide)):
            if res[i] in Neoantigen.CONSORTIUM_HYDROPHOBIC_RESIDUES:
                wc=wc+1
        wc=wc/len(self.mtPeptide)
        return wc 

    def getWeight(self):
        '''
        Returns 0 for neoantigens that mutated from a non-hydrophobic residues on position 2 or 9;
        these are excluded from analysis. Returns 1 for all other neoantigens
        '''
        w=1
        if len(self.wtPeptide) == len(self.mtPeptide):
            #[res1,res2]=filter(lambda el: el[0]!=el[1],zip(self.wtPeptide,self.mtPeptide))[0]
            res = list(zip(self.wtPeptide,self.mtPeptide))
            w=1
            for i in range(0, len(self.wtPeptide)):
                if res[i][0] != res[i][1]:
                    self.residueChange=Neoantigen.residueChangeClass(res[0][0], res[0][1])
                    self.position=i
                    self.position=self.position+1
                    if len(self.wtPeptide) == 9:
                        if self.residueChange[0]!="NH" and (self.position==2 or self.position==9):
                            w=0
                    elif len(self.wtPeptide) == 10:
                        if self.residueChange[0]!="NH" and (self.position==3 or self.position==10):
                            w=0
                    elif len(self.wtPeptide) == 11:
                        if self.residueChange[0]!="NH" and (self.position==4 or self.position==11):
                            w=0
        else:
            print("mismatching lengths")
        return w
        #### TO DO: Update for HLA C binding residues. Check if 11mers keep the same patterns
            
    def correctWT(self):
        '''
        Corrects large wildtype kD dissociation constants
        '''
        kd=self.wtkD
        prb=Neoantigen.W/(Neoantigen.W+kd)
        pru=kd/(Neoantigen.W+kd)
        eps=Neoantigen.WEPS
        prb+=eps
        pru+=eps
        z=1+2*eps
        prb/=z
        pru/=z        
        return pru/prb
                                
    def setA(self):
        '''
        Computes MHC amplitude A 
        '''
        self.A=Neoantigen.M/self.kD*self.correctWT()

    def getA(self):
        '''
        Return MHC amplitude A
        '''
        return self.A
        #print(self.A)
