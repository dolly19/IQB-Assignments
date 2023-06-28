
#Sequence
seq="SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDLLIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF"  #Our original sequence

#alpha propensities for each amino acids
pAlpha = {'A':1.45,
          'C':0.77,
          'D':0.98,
          'E':1.53,
          'F':1.12,
          'G':0.53,
          'H':1.24,
          'I':1.00,
          'K':1.07,
          'L':1.34,
          'M':1.20,
          'N':0.73,
          'P':0.59,
          'Q':1.17,
          'R':0.79,
          'S':0.79,
          'T':0.82,
          'V':1.14,
          'W':1.14,
          'Y':0.61}  

#beta propensities for each amino acids
pBeta = {'A':0.97,
          'C':1.30,
          'D':0.80,
          'E':0.26,
          'F':1.28,
          'G':0.81,
          'H':0.71,
          'I':1.60,
          'K':0.74,
          'L':1.22,
          'M':1.67,
          'N':0.65,
          'P':0.62,
          'Q':1.23,
          'R':0.90,
          'S':0.72,
          'T':1.20,
          'V':1.65,
          'W':1.19,
          'Y':1.29}

class prediction:
  def __init__(self):
    self.n = len(seq)
    self.helixList = [0]*self.n
    self.strandList = [0]*self.n
    self.outputList =[0]*self.n
  

  #Function to look in the left and right of the residue to find the complete stretch for predicting helix
  def findHelix(self,idx):
    i, j = idx, idx+5

    #search on left
    while(i>0):
      scoreL = pAlpha.get(seq[i-1])+ pAlpha.get(seq[i])+ pAlpha.get(seq[i+1])+ pAlpha.get(seq[i+2])
      if(scoreL>4):
        i=i-1
      else:
        break

    #search on right
    while(j<self.n-1):
      scoreR = pAlpha.get(seq[j+1])+ pAlpha.get(seq[j])+ pAlpha.get(seq[j-1])+ pAlpha.get(seq[j-2])
      if(scoreR>4):
        j=j+1
      else:
        break

    #assign the 'H' to the fragment
    for assign in range(i,j+1):
      self.helixList[assign] = 'H'


  #Function to scan the stretch of 6 residues that have atleast 4 residues with P(H)>1
  def findHelixResidue(self):
    for i in range(self.n-5):
      count=0
      for j in range(i,i+6):
        propensity = pAlpha.get(seq[j])
        if(propensity>=1):
          count += 1
      if(count>=4):
        self.findHelix(i)


  #Function to look in the left and right of the residue to find the complete stretch for predicting beta strand
  def findStrand(self,idx):
    i, j = idx, idx+4
    #search on left
    while(i>0):
      scoreL = pBeta.get(seq[i-1])+ pBeta.get(seq[i])+ pBeta.get(seq[i+1])+ pBeta.get(seq[i+2])
      if(scoreL>4):
        i=i-1
      else:
        break

    #search on right
    while(j<self.n-1):
      scoreR = pBeta.get(seq[j+1])+ pBeta.get(seq[j])+ pBeta.get(seq[j-1])+ pBeta.get(seq[j-2])
      if(scoreR>4):
        j=j+1
      else:
        break

    #assign the 'S' to the fragment
    for assign in range(i,j+1):
      self.strandList[assign] = 'S'


  #Function to scan the stretch of 5 residues that have atleast 3 residues with P(H)>1
  def findStrandResidue(self):
    for i in range(self.n-4):
      count=0
      for j in range(i,i+5):
        propensity = pBeta.get(seq[j])
        if(propensity>=1):
          count = count+1
      if(count>=3):
        self.findStrand(i)


  #Function to remove the conflicting regions in helix and beta strand Lists
  def removeConflicting(self):
    for i in range(self.n):
      if(self.helixList[i]=='H' and self.strandList[i] == 'S'):
        start =i
        end = i
        helixScore=0
        strandScore=0

        while(end<self.n and (self.helixList[end]=='H' and self.strandList[end] == 'S')):
          helixScore += pAlpha.get(seq[end])
          strandScore += pBeta.get(seq[end])
          end += 1
        
        if(strandScore>helixScore):
          for j in range(start,end):
              self.helixList[j] = '-'
        else:
          for j in range(start,end):
              self.strandList[j] = '-' 

  #Function to display the final output from helix and strand Lists
  def finalOutput(self):
    for i in range(self.n):
      if(self.helixList[i] == "H"):
          self.outputList[i] = "H"
      elif(self.strandList[i] == "S"):
          self.outputList[i] = "S"
    
    for i in range(self.n):
      if(self.outputList[i] == 0):
          self.outputList[i] = "-"

predClass = prediction()
predClass.findHelixResidue()
predClass.findStrandResidue()
predClass.removeConflicting()
predClass.finalOutput()

#Generating the required output
print(seq)
print("|"*len(seq))
for i in range(len(seq)): 
  print(predClass.outputList[i],end = "")

