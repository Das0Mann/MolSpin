import pandas as pd
import multiprocessing
from multiprocessing.pool import ThreadPool
import subprocess
import math
import os
import sys

s_TotalSteps = []
s_TotalStepsPassed = 0
s_SeachAll = True
s_FileIndex = 0

class Variable:
    def __init__(self, name : str, min : float, max : float, hireachy, StepSize, NumSteps, SyncVariable = ""):
        self.name = name
        self.min = min
        self.max = max
        if(StepSize != float(StepSize)):
            StepSize = (max - min)/(NumSteps-1)
        self.stepsize = StepSize
        if(math.isnan(NumSteps)):
            NumSteps = ((max - min)/self.stepsize) + 1
        self.numsteps = int(NumSteps) 
        self.step = 0
        self.hireachy = hireachy
        self.SyncedVariable = SyncVariable

class MSDVariable:
    def __init__(self, path, name, value, line):
        self.path = path
        self.name = name
        self.value = value
        self.line = line

def ReadMSDFile(BaseMSDFile : str):
    depth = 0
    buffer = ""
    InObject = False
    LineNum = 0
    VarString = []
    VariableList = []
    VariableLines = []
    InString = False
    InVariable = False
    f = open(BaseMSDFile, "r")
    file = []
    read = True
    for line in f:
        read = True
        LineNum += 1
        line = line.strip("\n")
        line = line.strip("\t")
        line = line.lstrip()

        #print(line)
        file.append("")
        if(depth == 2):
            InObject = True
        else:
            InObject = False
        for c in line:
            if buffer == "//":
                read = False
                buffer = ""
            if read == False:
                continue
            if c == '"':
                InString = not InString
            if(InString):
                buffer += c
                continue
            if c == ' ' and InObject == False:
                if buffer == "":
                    continue
                file[LineNum-1] = file[LineNum-1] + buffer + ' '
                buffer = ""
            elif c == ' ':
                continue
            elif c == '{':
                depth = depth + 1
                file[LineNum-1] += '{'
                VarString.append(buffer)
                buffer = ""
            elif c == '}':
                depth = depth - 1
                file[LineNum-1] += '}'
                VarString.pop()
                buffer = ""
            elif c == '=':
                VarString.append(buffer)
                file[LineNum-1] += buffer + '='
                VariableList.append(MSDVariable(VarString[0] + '.' + VarString[1] + '.' + VarString[2], VarString[2], "", LineNum-1))
                VariableLines.append(LineNum-1)
                buffer = ""
                VarString.pop()
                InVariable = True
            elif c == ';':
                if InVariable:
                    VariableList[-1].value = buffer
                    InVariable = False
                file[LineNum-1] += buffer + ';'
                buffer = ""
            else:
                buffer += c
        if(not InString):
            file[LineNum-1] += buffer
        else:
            LineNum = LineNum - 1

    removed = 0
    size = len(file)
    lineupdate = []
    for i in range(0,size):
        lineupdate.append([i])
    for i in range(0,size):
        if file[i-removed] == '':
            file.pop(i-removed)
            removed += 1
            continue
        lineupdate[i].append(i-removed)
        #print(file[i-removed])
    
    for i in range(0,len(VariableLines)):
        VariableLines[i] = lineupdate[VariableLines[i]][1]
    for i in range(0, len(VariableList)):
        VariableList[i].line = lineupdate[VariableList[i].line][1]
    
    f.close()
    
    return file, VariableList, VariableLines
            

def WriteMSDFile(file, Varlist, VariableLines, filename):
    f = open("msd-files/" + filename + ".msd", "w")
    
    line = 0
    depth = 0
    var = 0
    for l in file:
        if(l == '}'):
            depth = depth -1

        for i in range(0,depth):
            f.write("\t")
        
        if l == '{':
            depth = depth + 1

        if line in VariableLines:
            line += 1
            f.write(Varlist[var].name + "=" + Varlist[var].value + ";")
            var += 1
            f.write("\n")
            continue
        else:
            f.write(l + "\n")
        line += 1
    f.close()

def UpdateVar(vCSV : Variable):
    global s_TotalStepsPassed
    if s_SeachAll:
        MinSteps = 1
        for i in range(vCSV.hireachy+1, len(s_TotalSteps)):
            MinSteps *= s_TotalSteps[i]
        if(s_TotalStepsPassed == 0):
            val = vCSV.min
            return str(val)
        #print(s_TotalStepsPassed % MinSteps)
        if(s_TotalStepsPassed % MinSteps == 0):
            val = (((s_TotalStepsPassed / MinSteps) - math.floor(s_TotalStepsPassed/(vCSV.numsteps*MinSteps)) * vCSV.numsteps) * vCSV.stepsize) + vCSV.min
        else:
            val = ((math.floor(s_TotalStepsPassed/MinSteps) - math.floor(s_TotalStepsPassed / (vCSV.numsteps*MinSteps)) * vCSV.numsteps) * vCSV.stepsize) + vCSV.min
        return str(val)

def MSDFileModifer(BaseMSDFile : str, MSDFile : str, Varlist : list, file, var , varline):

    def FindVar(name : str):
        for v in var:
            if(v.path == name):
                return v.value

    global s_FileIndex
    for v in Varlist:
        for v2 in var:
            if(v.name == v2.path):
                if(v.SyncedVariable == ""):
                    v2.value = UpdateVar(v)
                else:
                    v2.value = FindVar(v.SyncedVariable)
                continue
                #print(v2.value)
    for v in Varlist:
        if v.name != "logfile" and v.name != "datafile":
            continue
    for v2 in var:
        if(v2.path.find("logfile") == -1 and v2.path.find("datafile") == -1):
            continue
        if(v2.path.find("logfile") != -1):
            v2.value = '"' + MSDFile + "-" + str(s_FileIndex) + '.log"'
        else:
            v2.value = '"' + MSDFile + "-" + str(s_FileIndex) + '.dat"'
    WriteMSDFile(file, var, varline, MSDFile + "/" + BaseMSDFile.lstrip("msd-files/").rstrip(".msd") + "-" + str(s_FileIndex))
    s_FileIndex += 1
    return 0

def GetVariablesToModify(FilePathCSV : str):
    global s_TotalSteps
    csv = pd.read_csv(FilePathCSV)
    df = pd.DataFrame(csv)
    #print(df)
    coldata = []
    columns = df.columns.values
    for i in range(0,len(columns)):
        coldata.append([])
        for e in range(0,len(df[columns[i]])):
            coldata[i].append(df[columns[i]][e])
    #print(coldata)
    Variables = []
    for i in range(0, len(coldata[0])):
        if(not isinstance(coldata[5][i], str)):
            Variables.append(Variable(coldata[0][i], coldata[1][i], coldata[2][i], i, coldata[3][i], coldata[4][i]))
            s_TotalSteps.append(Variables[i].numsteps)
        else:
            Variables.append(Variable(coldata[0][i], coldata[1][i], coldata[2][i], i, coldata[3][i], coldata[4][i], coldata[5][i]))
    return Variables

def Work(i,file):
    r = subprocess.Popen(["./molspin", "-p", "2", file + "-" + str(i) + ".msd"])
    r.wait()

def main():
    global s_TotalStepsPassed
    args = sys.argv
    file = ""
    if(len(args) == 1):
        file = input("Name of file to modify: ")
    else:
        file = args[1]
    var = GetVariablesToModify("VariableControl.csv")
    TotalStepsProd = 1
    filestring = "msd-files/" + file
    if not os.path.exists("msd-files/" + file):
        os.makedirs("msd-files/" + file)
    for i in range(0,len(s_TotalSteps)):
        TotalStepsProd = TotalStepsProd * s_TotalSteps[i]
    filedat, vardat, varlinedat = ReadMSDFile(filestring + ".msd")
    for i in range(0,TotalStepsProd):
        MSDFileModifer(filestring + ".msd", file, var, filedat, vardat, varlinedat)
        s_TotalStepsPassed += 1
    num = int(math.floor(multiprocessing.cpu_count()/2.0))
    print(num)
    tp = ThreadPool(num)
    print("Done")
    work = []
    for i in range(0,s_FileIndex):
        r = tp.apply_async(Work,(i,"msd-files/" + file + "/" + file,))
        work.append(r)
    done = 0
    doneThreads = []
    while(done != len(work)):
        index = 0
        for t in work:
            if(index in doneThreads):
                index = index + 1
                continue
            t.wait(1)
            if(t.ready()):
                done = done + 1
                doneThreads.append(index)
            index = index + 1
    tp.close()
    tp.join()

main()
