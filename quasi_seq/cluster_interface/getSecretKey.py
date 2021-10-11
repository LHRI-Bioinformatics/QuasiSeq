from sys import argv
from django.utils.crypto import get_random_string
import os.path

def creatSecretKey(path):
	keyFile=path+'/secretkey.txt'
	if not os.path.isfile(keyFile):
		chars = 'abcdefghijklmnopqrstuvwxyz0123456789!@#$%^&*(-_=+)'
		SECRET_KEY = get_random_string(50, chars)
		outfile = open(keyFile, "w")
		outfile.writelines(SECRET_KEY)
		outfile.close()

creatSecretKey(argv[1])
