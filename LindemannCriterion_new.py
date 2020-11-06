from sys import argv
from numpy import zeros,sqrt,mean
file=open(argv[1],"r")
natoms=int(argv[2])
nlines=0
n_permut=natoms*(natoms-1)/2
r=zeros((n_permut),float)
r2=zeros((n_permut),float)
liner=zeros((500000),float)
liner2=zeros((500000),float)

for line in range(500000):
    for i in range(n_permut):
        dummy=file.readline()
        if (dummy=="\n"):
            break
        if (dummy == ""):
            break
        r[i]+=float(dummy)
        r2[i]+=float(dummy)*float(dummy)
        liner[line]+=float(dummy)
        liner2[line]+=float(dummy)*float(dummy)
    if (dummy=="\n"):
        break
    if (dummy==""):
        break
    nlines+=1
file.close()
print "Number of Lines: ",nlines
d_L=0.0
d_LR = 0.0
d_LQ=0.0
for i in range(n_permut):
    d_L+=sqrt(((r2[i]/float(nlines))-(r[i]/float(nlines))*(r[i]/float(nlines))))/(r[i]/float(nlines))

d_L=d_L/n_permut

r2_tot=0.0
r_tot=0.0
liner2_tot=0.0
liner_tot=0.0
for i in range(n_permut):
    r2_tot+=r2[i]
    r_tot+=r[i]

for line in range(nlines):
    liner_tot+=(liner[line]/n_permut)*(liner[line]/n_permut)
    liner2_tot+=liner2[line]*liner2[line]

liner2_tot=liner2_tot/(nlines)
r2_tot=r2_tot/(n_permut*nlines)
r_tot=r_tot/(n_permut*nlines)

d_LR=sqrt(r2_tot-r_tot*r_tot)/r_tot
d_LQ=sqrt(natoms-1)*sqrt(liner2_tot-(r2_tot*n_permut)*(r2_tot*n_permut))/(r2_tot*n_permut)
d_LF=sqrt((liner_tot/nlines)-(r_tot)*(r_tot))/(r_tot)
print "Distinguishable :    d_L=",d_L
print "Indistinguishable : d_LR=",d_LR
print "Quantum           : d_LQ=",d_LQ
print str(liner_tot/nlines)+" "+str(r_tot)
print "Fluctuation       : d_LF=",d_LF

