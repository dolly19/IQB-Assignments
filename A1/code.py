# -*- coding: utf-8 -*-
"""iQB A1.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1XXIRiQpHJhNI_7cgUKVFT-6X7wUConhm
"""

import numpy as np

#Given DNA Sequences
Seq1 = "GATGCGCAG" 
Seq2 = "GGCAGTA"

#Given scoring value
Match = 2 
Mismatch = -3
Gap = -1

#initializing the number of rows and columns for the DP table
m = len(Seq2) + 1
n = len(Seq1) + 1

#function to fill matrix for the global alignments 
def globalAlignment(dp):

#initializing the 0th row and 0th column in DP table
  for i in range(m):
    dp[i][0] = Gap*i
  for j in range(n):
    dp[0][j] = Gap*j

  for i in range(1, m):
    for j in range(1, n):
      if Seq2[i-1] == Seq1[j-1]:
        diagonal = dp[i-1][j-1]+Match  
      else:
        diagonal = dp[i-1][j-1]+Mismatch  
      vertical = dp[i-1][j] + Gap   
      horizontal = dp[i][j-1] + Gap   
      dp[i][j] = max(diagonal, vertical, horizontal)

  return dp

#function to print matrix
def printMatrix(dp):
  S = '_'+Seq1
  T = '_'+Seq2

  for i in range(m+1):
    for j in range(n+1):
      if i == 0 and j == 0:
        print('X', end = '\t')
      elif i == 0:
        print(S[j-1], end = '\t')
      elif j== 0:
        print(T[i-1], end = '\t')
      else:
        print(dp[i-1][j-1], end = '\t')
    print()

#function to traceback for finding the all optimal alignments with their scores in case of global alignment
def globalBacktrack(seq1,seq2,align1,align2,matrix,i,j,score):

  if(i==0 and j==0):
    print(align1)
    print(align2)
    print(f"Score is: {score}")
    print("------------------------------")
    print()
    return

  if (i!=0 and j!=0 and seq2[i-1] == seq1[j-1] and (matrix[i][j] == (matrix[i-1][j-1]+Match))):
    globalBacktrack(seq1, seq2, seq1[j-1]+align1, seq2[i-1]+align2, matrix, i-1,j-1,score+Match)

  if (i!=0 and j!=0 and seq2[i-1] != seq1[j-1] and (matrix[i][j] == (matrix[i-1][j-1]+Mismatch))):
    globalBacktrack(seq1, seq2, seq1[j-1]+align1, seq2[i-1]+align2, matrix, i-1,j-1,score+Mismatch)
  
  if(i!=0 and matrix[i][j] == matrix[i-1][j]+Gap):
    globalBacktrack(seq1, seq2, "-"+align1, seq2[i-1]+align2, matrix, i-1, j,score+Gap)
   
  if(j!=0 and matrix[i][j] == matrix[i][j-1]+Gap):
    globalBacktrack(seq1, seq2, seq1[j-1]+align1, "-"+align2, matrix, i, j-1, score+Gap)

#function to fill matrix for the local alignments 
def localAlignment(dp):
  for i in range(1, m):
    for j in range(1, n):
      if Seq2[i-1] == Seq1[j-1]:
        diagonal = dp[i-1][j-1]+Match  
      else:
        diagonal = dp[i-1][j-1]+Mismatch  
      vertical = dp[i-1][j] + Gap   
      horizontal = dp[i][j-1] + Gap   
      dp[i][j] = max(diagonal, vertical, horizontal,0)
  return dp

#function to traceback for finding the all optimal alignments with their scores in case of local alignment
def localBacktrack(seq1,seq2,align1,align2,matrix,i,j,score):

  if(matrix[i][j]==0):
    print(align1)
    print(align2)
    print(f"Score is: {score}")
    print("------------------------------")
    print()
    return

  if (i!=0 and j!=0 and seq2[i-1] == seq1[j-1] and (matrix[i][j] == (matrix[i-1][j-1]+Match))):
    localBacktrack(seq1, seq2, seq1[j-1]+align1, seq2[i-1]+align2, matrix, i-1,j-1,score+Match)

  if (i!=0 and j!=0 and seq2[i-1] != seq1[j-1] and (matrix[i][j] == (matrix[i-1][j-1]+Mismatch))):
    localBacktrack(seq1, seq2, seq1[j-1]+align1, seq2[i-1]+align2, matrix, i-1,j-1,score+Mismatch)
  
  if(i!=0 and matrix[i][j] == matrix[i-1][j]+Gap):
    localBacktrack(seq1, seq2, "-"+align1, seq2[i-1]+align2, matrix, i-1, j,score+Gap)
   
  if(j!=0 and matrix[i][j] == matrix[i][j-1]+Gap):
    localBacktrack(seq1, seq2, seq1[j-1]+align1, "-"+align2, matrix, i, j-1, score+Gap)

#function to find index of maximum value in the matrix for local alignment
def findMaxIndex(matrix):
  row=0
  col=0
  max=0
  for i in range(1, m):
      for j in range(1, n):
          if(max<matrix[i][j]):
              max,row,col= matrix[i][j],i,j
  return row,col

"""**Question 1**"""

#initializing the DP table with 0's
dp = np.zeros((m, n), dtype = np.int)

scoringMatrix = globalAlignment(dp)
print(f"The matrix for global alignment: ")
print()
printMatrix(scoringMatrix)

print(f"All the optimal alignments with their scores: ")
print()
globalBacktrack(Seq1, Seq2, "", "", scoringMatrix, m-1, n-1,0)

"""**Question 2**"""

#updated scores for this question
Match = 2
Mismatch = -1
Gap = -3

dp = np.zeros((m, n), dtype = np.int)

newScoringMatrix = globalAlignment(dp)
print(f"The matrix for global alignment: ")
print()
printMatrix(newScoringMatrix)

print(f"All the optimal alignments with their scores: ")
print()
globalBacktrack(Seq1, Seq2, "", "", newScoringMatrix, m-1, n-1,0)

"""**Question 3**"""

Seq1 = "GATGCGCAG" 
Seq2 = "GGCAGTA"
Match = 2 
Mismatch = -1
Gap = -3
m = len(Seq2) + 1
n = len(Seq1) + 1
dp = np.zeros((m, n), dtype = np.int)

localScoringMatrix = localAlignment(dp)
print(f"The matrix for local alignment: ")
print()
printMatrix(localScoringMatrix)

print(f"All the optimal alignments with their scores: ")
print()
row,col = findMaxIndex(localScoringMatrix)
localBacktrack(Seq1, Seq2, "", "", localScoringMatrix, row, col,0)