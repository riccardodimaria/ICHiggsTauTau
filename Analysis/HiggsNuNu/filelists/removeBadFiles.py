#!/usr/bin/python

import os


listOfLists = os.listdir("./161031/medium")
#listOfLists = os.listdir("./161031/short")

# import bad files
badFilesFileName = "failed_files_DATA.dat"
#badFilesFileName = "failed_files_MC.dat"
badFilesFile = open( badFilesFileName )

badFiles = badFilesFile.read().split("\n")

# search in lists
for l in listOfLists :
    f=open("./161031/medium/"+l)
    fOut=open("./161031/noZombie/"+l, 'w+')
    
    content = f.read()
    
    allFiles = content.split("\n")
    
    for f in allFiles :
        if f != '' :
            if f in badFiles :
                print f+" is a bad file"
            else :
                fOut.write(f+'\n')
            
