
import numpy as np
from random import *
import sys


def read_file(filename):
	""" the sequence file to be raed """

	fichier = open(filename,"r")
	return fichier.readline().rstrip('\n')




def create_matrix(length):
	"""Create a matrix that contains all the amino acids"""

	mat = np.char.array([['*']*length*2]*length*2)
	return mat

def energy_calculator(matrix,seq,dic):
	"""A function that calculates the total energy of a conformation"""
	e = 0
	for a in range(0,len(seq)):
		if seq[a] != 'h':
			e = e
		else:
			(i,j) = dic[a]
			for k in range(i-1,i+2):
				for l in range(j-1,j+2):
					if matrix[k][l] == 'h':
						e = e - 1
			e = e + 1
	return e

def create_file(matrix,filename):
	"""Write data in file"""

	np.savetxt(filename, matrix, fmt='%-10.1c', delimiter='')
	return


if __name__ == '__main__':
	#is nstep > 1?
	if (int(sys.argv[2]) < 1):
		print('must be > 1 ')
		sys.exit()
	nstep = int(sys.argv[2]) #number of stemps which is the 2nd arg

	#Import the sequence
	sequence = sys.argv[1]
	seq = read_file(sequence)


	#Create a matrix
	seqlen = len(seq)
	mat = create_matrix(seqlen)

	j=int(seqlen/2)

	#locate the sequence in the middle
	for k in range(seqlen):
		mat[seqlen][j]=seq[k]
		j=j+1

	#Create a dict containing the coordinates of aa
	list_aa = range(seqlen)
	dic = {}
	j=int(seqlen/2)
	for i in range(seqlen):
		dic[list_aa[i]] = (seqlen,j)
		j=j+1


	#compute the energy
	Einit = energy_calculator(mat,seq,dic)
	print("initial energy is: " + str(Einit))

	#Monte Carlo

	for k in range(nstep):
		newDic = dic
		a = randint(0,len(list_aa)-1) #random pick of an aa
		(i,j) = newDic[a] #coordinates of the aa
		d = randint(1,8) #C there are 8 different possible moves, one will chosen randomly
		if d==1:
			if (mat[i][j-1] != '*'): #make sur the move can be done
				break
			else:
				for l in range(a,len(list_aa)):
					(i,j) = newDic[l]
					newDic[l] = (i,j-1)
		elif d==2:
			if (mat[i][j+1] != '*'):
				break
			else:
				for l in range(a,len(list_aa)):
					(i,j) = newDic[l]
					newDic[l] = (i,j+1)
		elif d==3:
			if (mat[i-1][j] != '*'):
				break
			else:
				for l in range(a,len(list_aa)):
					(i,j) = newDic[l]
					newDic[l] = (i-1,j)
		elif d==4:
			if (mat[i+1][j] != '*'):
				break
			else:
				for l in range(a,len(list_aa)):
					(i,j) = newDic[l]
					newDic[l] = (i+1,j)
		elif d==5:
			if (mat[i][j-1] != '*'):
				break
			else:
				for l in range(a,len(list_aa)):
					(i,j) = newDic[l]
					newDic[l] = (i-1,j-1)
		elif d==6:
			if (mat[i][j+1] != '*'):
				break
			else:
				for l in range(a,len(list_aa)):
					(i,j) = newDic[l]
					newDic[l] = (i-1,j+1)
		elif d==7:
			if (mat[i-1][j] != '*'):
				break
			else:
				for l in range(a,len(list_aa)):
					(i,j) = newDic[l]
					newDic[l] = (i+1,j-1)
		elif d==8:
			if (mat[i+1][j] != '*'):
				break
			else:
				for l in range(a,len(list_aa)):
					(i,j) = newDic[l]
					newDic[l] = (i+1,j+1)


		#locate the new sequence in the middle
		newMat = create_matrix(seqlen)
		for k in range(len(seq)):
			(i,j) = newDic[k]
			newMat[i][j]=seq[k]

		#compute the new energy
		Enew = energy_calculator(newMat,seq,newDic)
		if Enew <= Einit:
			dic = newDic
			mat2 = newMat
			Einit = Enew

	#results
	create_file(mat2,'results.txt')
	print(mat)
	print("final energy_calculator is: " + str(Einit))
	print(mat2)
