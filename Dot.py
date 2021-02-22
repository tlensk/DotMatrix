# # Copyright 2019 Dan Boley and Tatiana Lenskaia
import time as tm
import numpy as np
import matplotlib.pyplot as plt
import sys

#======================================
#Created: July 16 2019 by boley
# Dot.py:   Display dot plot showing common motifs between pairs of genome sequences
# Sample command lines:
#    python3 -i Dot.py 15 MT118835.fasta MT300186.fasta  --labels 100 --markersize [1,10] 
#    python3 -i Dot.py 40 MultiGenome.fasta ''  --labels -10 

# Updates
# Tue Sep 15 09:24:52 CDT 2020 changed default to --labels -11
# Mon Oct  5 15:08:06 CDT 2020 changed default to --fullaxes.  added option '--tightaxes' to undo.
#======================================

# Process command line paramters

args=sys.argv
print('args=',args)
fOutLog = False # sys.stdout # open("DotTiming.txt", "x")
winlen = 40
InFile1 = False #  'NC_003444.fasta'
InFile2 = False #  "NC_020163.fasta"
OutFile = '--noout'
if len(args)>1:
   if args[1] == 'help':
      winlen = 40
   else:
      winlen = int(args[1]) 
   if len(args)>2:
       if args[2]: InFile1=args[2]
       if len(args)>3 and args[3]: InFile2=args[3]
       if len(args)>4 and args[4]: OutFile=args[4]
       if not InFile2: InFile2 = InFile1
else:
   args.append('help')
if OutFile=='-': OutFile = 'stdout'
if OutFile[0]=='-': OutFile = False
labels = ('--labels' in args) and int(args[1+args.index('--labels')])
if labels is False: labels=-11
loweronly = ('--loweronly' in args)
tightaxes = ('--tightaxes' in args)
nofill = ('--nofill' in args)
nogrid = ('--nogrid' in args)
noreverse = ('--noreverse' in args)
countonly = ('--countonly' in args)
rna2dna = ('--rna' in args)
nopropagate = ('--nopropagate' in args)
gap = ('--gap' in args)
if countonly: gap = 0
elif gap:     gap = int(args[args.index('--gap')+1])
else:         gap = 1
noplot = ('--noplot' in args) or countonly
marksizes = ('--markersize' in args) and eval(args[1+args.index('--markersize')])
if not marksizes: marksizes = [2,6]
 
fOutOverlap=False
if OutFile=='stdout':
   fOutOverlap = sys.stdout
elif OutFile:
   fOutOverlap = open(OutFile, "x")
elif countonly:
   Outfile = '-'
   fOutOverlap = sys.stdout

sep = "\t"

print('Command line options ("False" means "omitted"):')
print('    winlen: ....',winlen,'           sliding window length [a number, negative: gtype = circular]')
print('    InFile1: ...',InFile1,'        [input fasta file 1, shorter genomes such as viruses]')
print('    InFile2: ...',InFile2,'        [input fasta file 2, longer genomes such as bacteria]  if empty string, use InFile1.')
print('    OutFile: ...',OutFile,'        [optional output file to store the matches or counts]')
print('switches:')
print('  --labels: ....',labels,"         use labels of length up to this number on plot")
print('                          [a negative number: just show the ID, not the text.  Zero means no labels]')
print('                          [default: -11 meaning just print the accession number part of the FASTA header]')
print('  --loweronly: .',loweronly,  "   show only lower triangle of plots when comparing seqs with themselves")
print('  --tightaxes: .',tightaxes,  "   plot tight axes when matching strings are limited in range")
print('  --nofill: ....',nofill,     "   don't fill in dots to mark matching positions; instead show just starting points of matches")
print('  --nogrid: ....',nogrid,     "   omit the grid on the plots")
print('  --noplot: ....',noplot,     "   omit the plot altogether")
print('  --noreverse: .',noreverse,  "   omit the reverse complement search")
print('  --rna: .......',rna2dna  ,  "   translate RNA strings to DNA form ")
print('  --nopropagate ',nopropagate,"   don't propagate: just mark the ends of the matched substrings")
print('  --countonly: .',countonly,  "   just count the number of matches instead of showing anything")
print('  --gap: .......',gap,    "       mark single letter changes (when gap == 1, the default)")
print('  --markersize .',marksizes,   "  markersizes: list: [match,break] ")

gtype = "linear"
if winlen < 0:
   gtype = 'circular'
   winlen = -winlen
print('    gtype =',gtype,"       [are genomes 'linear' or 'circular'?]")
print('')

#if noplot and not fOutOverlap and not fOutLog: print ('\n=== This will generate no output: are you sure?\n')
if 'help' in args:
   print("\nusage: python3 -i Dot.py  WINDOW_LENGTH  SHORT_GENOMES.fasta  LONG_GENOMES.fasta  [OUTFILE]  switches")
   print("  For LONG-GENOMES, use '' to select the SHORT.GENOME (a self-intersection).") 
   print("  For OUTFILE, use '' for none, '-' for standard output.  Switches cannot be abbreviated.  Examples:")
   print("    python3 -i Dot.py 15 MT118835.fasta MT300186.fasta  --labels 100 --markersize [1,10] ")
   print("    python3 -i Dot.py 40 MultiGenome.fasta ''  --labels -10 ")
   print("    python3 -i Dot.py 20 LYTIC/lactphs.fasta '' '' --nogrid --labels -10 --loweronly\n")
   print("  Purpose: display dot plot: place dots showing positions of matching strings of length at least WINDOW_LENGTH:")
   print("    blue=forward match, ending w/ black; red=reverse match, ending w/ magenta. green=SNP:",gap,"letter difference)")
   sys.exit(1)

# Functions
def AddToDict(dict,item,value,default=0):
   "add item to dict if value is bigger that what is already in the dict"
   if item in dict:
      dict[item]=max(dict[item],value)
   else:
      dict[item]=max(default,value)
   return

def CreateDict(text, m, gtype = "l"):
    """Input: 
                texttext        (genomic sequence) 
                mm                      (screening window size)
                gtype           (type of genome: c = circular, l=linear)
                                second letter = 'r' if the the text is reverse complement.
    Output:
                dictionary  (keys: distinct strings of size mm; values: starting locations within text)"""
    
    if (m <= 0) or (m > len(text)):
        print("\n=== (n = "+str(m)+") is not a valid window size for this genome!!!");
        return {}
    
    d_g = dict()
    nn = len(text)
    print('len of text: ',nn,m,gtype,end="... ",flush=True)
    
    gtype = gtype.lower()
    if gtype[:1] == "c":
        text = text + text[:(m-1)]
        lastpos = nn
    elif gtype[:1] == "l":
        lastpos = nn-m+1
    else:
        print("\n=== Is this genome linear or circular?");
        return d_g
    
    for ii in range (lastpos):
        pos = ii
        if gtype[1:2] == 'r':
           pos = nn-ii-1
           if pos<0: pos=pos+nn
        bl = text[ii:(ii+m)]
        bl = bl.lower()
        #print(bl,pos,ii,text)
        if bl in d_g:
            d_g[bl].append(pos)
        else:
            d_g[bl] = [pos]
    if gtype[1:2] == 'r': print('reversed: ',end="")
    print('len of dict: ',len(d_g))
    return d_g


#======================================
#June 16, 2018
def FindIntersection(d_g11, d_g22):
        """Input: 
                        two dictionaries
        Output: 
                        list of keys appearing in both dictionaries """
        
        #return list(d_g11.keys() &  d_g22.keys());


        nn1 = len(d_g11);
        nn2 = len(d_g22);
        
        if nn1 <= nn2:
           t_ints = [s for s in d_g11 if s in d_g22]
        else:
           t_ints = [s for s in d_g22 if s in d_g11]

        #the following is the old code ... not sure which is more efficient.

        #t_ints = list()
        #if nn1 <= nn2:
        #   for s in d_g11:
        #       if s in d_g22:
        #           t_ints.append(s);
        #else:
        #   for s in d_g22:
        #       if s in d_g11:
        #           t_ints.append(s);

        return t_ints;

#========================================

def preprocess(line):
    "preprocess strings being read in: strip white space, make lower case, RNA->DNA if needed"
    line = line.strip()
    line = line.lower()
    if rna2dna: line = line.translate(str.maketrans('u','t'))
    return line;

    
def Rev_cmp(st):
    "compute the reverse complement of a DNA sequence"
    st = st.lower()
    cmp_st = st.translate(str.maketrans("acgt","tgca"))
    rev_cmp_st = cmp_st[::-1]
    return rev_cmp_st;
#===================================================


def AddToArrayOfDicts_FR(t_dicts, text, r_id,r_name,m=40, gtype = "l",fOut=False,t_vir=[]):
    """Input: 
        t_dicts = list: ([forward_dict_1,rev_dict_1],...,[for_dict_k,rev_dict_k]) or empty []
        text = new sequence 
        m = sliding window length
        gtype = 'l' (linear) or 'c' (circular)
    Output: 
        t_dicts = list: ([forward_dict_1,rev_dict_1],...,[for_dict_k+1,rev_dict_k+1])"""
    #print(r_id,r_name, len(text))
    if not r_id: r_id=InFile1+'_'+str(len(t_vir))
    r_size = len(text)
    d = CreateDict(text,m,gtype)
    if noreverse:
       r_text=""
       r_d = {}
    else:
       r_text = Rev_cmp(text)
       r_d = CreateDict(r_text,m,gtype[0]+'r')
    t_dicts.append((d,r_d))
    t_vir.append((r_id,r_name,r_size,len(d)))
    if fOut: print('ShortSeq:',r_id,r_name, r_size, len(d),len(r_d), sep = "\t", file = fOut)
    text=""
    r_id=""
    return t_vir

def FindIntersectionWithList(d_g11, t_list, l = 1,bact_serial=-1,CountOnly=False):
    """ find intersection of d_g11 with every entry in t_list
    Input:
        d_g11    :  a dictionary from CreateDict
        t_list   :  a list of dictionaries
        l = 1 means t_list is a list of individual dictionaries
        l = 2 means t_list is a list of pairs: ([forward_dict_1,rev_dict_1],...,[for_dict_k,rev_dict_k])
        bact_serial   : label or serial number for d_g11
    Output: list of intersections of the form [bact_serial, k, position in dictionary d_g11, position in dictionary t_list[k], +1 or -1]
          where k is the index of the dictionary being scanned, for k in range(len(t_list))
    """
    if l == 1:
        tt = []
        for i in range(len(t_list)):
            t=({},{})
            it = t_list[i]
            if loweronly and bact_serial>i: t_ints = []
            else:                           t_ints = FindIntersection(it, d_g11)
            if CountOnly:
               tt.append(len(t_ints))
            elif bact_serial>=0:
                for j in t_ints:
                    for k1 in d_g11[j]:
                        for k2 in it[j]:
                            AddToDict(t[0],(k1,k2),4)
                            AddToDict(t[0],(k1+1,k2+1),2)
                tt.append(t)
            else: tt.append(t_ints)
        return tt
    elif l == 2:
        tt = []
        for i in range(len(t_list)):
            t=({},{})
            it = t_list[i]
            if loweronly and bact_serial>i:
               t_ints = []
               r_t_ints = []
            else:
               t_ints = FindIntersection(it[0], d_g11)
               r_t_ints = FindIntersection(it[1], d_g11)
            if CountOnly:
               tt.append((len(t_ints),len(r_t_ints)))
            elif bact_serial>=0:
                for j in t_ints:
                    for k1 in d_g11[j]:
                        for k2 in it[0][j]:
                            AddToDict(t[0],(k1,k2),4)
                            AddToDict(t[0],(k1+1,k2+1),2)
                for j in r_t_ints:
                    for k1 in d_g11[j]:
                        for k2 in it[1][j]:
                            AddToDict(t[1],(k1,k2),4)
                            AddToDict(t[1],(k1+1,k2-1),2)
                tt.append(t)
            else: tt.append((t_ints,r_t_ints))
        return tt
    else:
        print("\n=== Are we computing intersection for only forward (l=1) or for both forward and reverse (l=2)?")
        return []
    
#Create an array of dictionaries for viruses

t_rep = [] # holds the labels for the bacteria (2nd parameter)
t_vir = [] # holds the labels for the viruses  (1st parameter)


t_Dv = []  # holds the list of dictionaries for the viruses
m = winlen # 40 #15  # the sliding window length (with typical values

fIn = open(InFile1,"r")

text = ""  # temporary to hold the DNA sequence itself
r_id = ""   # temporary to hold the ID of the sequence (accession number?)
r_name = "" # temporary to hold other info on the sequence


#fOut = open("2vir.txt","x")

# Now fill in the dictionaries for the viruses and put them into a list
tm_st = tm.time()
for line in fIn:
    line = preprocess(line)
    #print('len of line: ',len(line),': ',line[0:10])
    if line == "":
        if text != "": 
           t_vir = AddToArrayOfDicts_FR(t_Dv, text, r_id,r_name,m, gtype,fOut=fOutLog,t_vir=t_vir)
           text=""
    elif line[0] == ">":
        if text != "": 
           t_vir = AddToArrayOfDicts_FR(t_Dv, text,r_id,r_name, m, gtype,fOut=fOutLog,t_vir=t_vir)
           text=""
        line = line[1:]
        while line[:1]==' ': line=line[1:]
        t_info = line.split(" ",1)
        r_id = t_info[0].split(".")[0]
        print('AddVirus:',line)
        if not r_id: r_id=str(r_id)+InFile1+'_'+str(len(t_vir))
        r_name = len(t_info)>1 and t_info[1]
    else:
        text = text+line
if text != "":
        t_vir = AddToArrayOfDicts_FR(t_Dv, text,r_id,r_name, m, gtype,fOut=fOutLog,t_vir=t_vir)
        text=""
fIn.close()

#print(len(t_Dv))
#print(t_Dv[0])
tm_fn = tm.time()
#fOut.close()

print('To make list of dictionaries, Time:',tm_fn-tm_st)
if fOutLog: print("Array of viral dictionaries:", len(t_Dv),'Time: ',tm_fn-tm_st, sep = sep,file=fOutLog)
if fOutLog: print("Time RepliconID", "RepliconSize,bp", "Time,sec",sep = sep,file=fOutLog)


#Find intersection between a prokaryote replicon and each viral dictionary from the array 


text = ""
r_id = ""
r_name = ""
row_ind=-1

Overlaps = [] # 2d array of intersections #bacteria by #viruses
# Overlaps[b][v] == (intersection_for_forward_strand,intersection_for_reverse_strand) == a tuple of intersections.
# For [bacteria,virus] pair, an intersection is a dictionary with entries of the form (k1,k2): flag
# where k1 is the position of a string (length==winlen) in the bacteria
#       k2 is the position of the same string in the virus
#       flag = 4
# later flag has value 2 to denote the end of a sequence of consecutive positions of matches
#       flag has value 3 if there is a gap of 1 position to the next match (isolated 1 letter difference)
#          the number 1 is the default.
# if --countonly is set, then Overlaps[b][v] == (#common_strings_forward,#common_strings_reverse) == a tuple of integers

dummydict={1:2}
#Dcts=[]
tm_fn1=tm_st
bact_serialno=-1

# define some functions to compute the dictionaries and do the intersection on the fly.
# For each bacteria, assmemble the dictionary, find the intersection, then throw away the dictionary, keeping only the intersection.
# The intersection is
def TextIntersectWithList(text,m,gtype,t_Dv,l=2,bact_counter=bact_serialno,t_repl=t_rep,Count_Only=False):
        """ compute dictionary for genome "text" and find intersection with all genomes in "t_Dv"
        the arguments are like those of FindIntersection.
           Bact_counter  is an identifying index for the genome
           t_repl is a list of statistics for the genome
           Count_Only: true if you want only a count of how many intersections, O.W. false. 
           This fcn appends the result to Overlaps, and returns nothing.
        """
        #print(r_id,r_name, len(text),t_rep)
        r_size = len(text)
        tm_st1 = tm.time()
        if fOutLog: print("Replicon:",r_id,r_name, r_size, sep = "\t", file = fOutLog)
        if type(text) == type('string'):  # if the dict has already been computed....
            Db = CreateDict(text,m, gtype)
        elif type(text) == type(dummydict):
            Db = text
        elif l==2 and (type(text)==type(()) or type(text)==type([])):
            Db = text[0]
        else: print('\n=== need a string or a dictionary',type(text))
        db_size = len(Db)
        #Dcts.append(Db)
        tm_db = tm.time()
        ttt = FindIntersectionWithList(Db,t_Dv,l=2,bact_serial=bact_counter,CountOnly=Count_Only)
        t_repl.append((r_id,r_name,r_size,db_size))
        tm_fn1 = tm.time()
        if fOutLog: print('To make dict #',bact_counter,r_size,db_size,' & find intersections to all in list, Time:', tm_db-tm_st1, tm_fn1-tm_db,sep = sep,file=fOutLog)
        print('To make dict #',bact_counter,r_size,db_size,' & find intersections to all in list, Time:', tm_db-tm_st1, tm_fn1-tm_db,sep = sep)


        #row_ind = [x[0] for x in t_repl].index(r_id)
        #print(row_ind)
        Overlaps.append(ttt)
        return


#fOut = open("2gen_meta.txt","x")

text=""
if InFile2 == InFile1:  # the two genomes are the same, so just use the existing dictionaries

   print('The two genomes are the same, hence using the same dictionaries')
   tm_st1 = tm.time()
   for dct in t_Dv:
       bact_serialno+=1
       #print('Main1: bact_serial',bact_serialno)
       TextIntersectWithList(dct,m,gtype,t_Dv,2,bact_serialno,t_rep,Count_Only=countonly)
   t_rep=t_vir
   tm_fn1 = tm.time()
   if fOutLog: print('To intersect:',r_id, 'Time:' , tm_fn1-tm_st1,sep = sep,file=fOutLog)

   #print(row_ind)

else:

   ## read in each individual genome and compute intersections with all viruses w/out storing the dictionary.
   fIn = open(InFile2,"r")
   text = ""
   for line in fIn:
       line = preprocess(line)
       if line == "":
           if text != "":
               bact_serialno+=1
               if not r_id: r_id = InFile2+'_'+str(bact_serialno)
               #print('Main1: bact_serial',bact_serialno)
               TextIntersectWithList(text,m,gtype,t_Dv,2,bact_serialno,t_rep,Count_Only=countonly)
               text = ""
               r_id = ""
       elif line[0] == ">":
           if text != "":
               bact_serialno+=1
               if not r_id: r_id = InFile2+'_'+str(bact_serialno)
               #print('Main2: bact_serial',bact_serialno)
               TextIntersectWithList(text,m,gtype,t_Dv,2,bact_serialno,t_rep,Count_Only=countonly)
               text = ""
               r_id = ""
           line = line[1:]
           while line[:1]==' ': line=line[1:]
           t_info = line.split(" ",1)
           r_id = t_info[0].split(".")[0]
           if not r_id: r_id = str(r_id)+InFile2+'_'+str(bact_serialno)
           r_name = len(t_info)>1 and t_info[1]
       else:
           text = text+line
   if text != "":
           bact_serialno+=1
           if not r_id: r_id = str(r_id)+InFile2+'_'+str(bact_serialno)
           #print('Main3: bact_serial',bact_serialno)
           TextIntersectWithList(text,m,gtype,t_Dv,2,bact_serialno,t_rep,Count_Only=countonly)
           text=""
           r_id = ""

   fIn.close()

# Now find all gaps of length 'gap' [default 1], marking single letter changes.

def FindGaps(Pair,gap,m=winlen):
    for (k1,k2) in Pair[0]:
        if Pair[0][(k1,k2)]==2:
           if (k1+m+gap,k2+m+gap) in Pair[0]:
               Pair[0][(k1,k2)]=3
               #print('gapf:',b,v,k1,k2)
    for (k1,k2) in Pair[1]:
        if Pair[1][(k1,k2)]==2:
           if (k1+m+gap,k2-m-gap) in Pair[1]:
               Pair[1][(k1,k2)]=3
               #print('gapr:',b,v,k1,k2)
    return Pair
if gap:
   tm_stg=tm.time()
   nrows=len(Overlaps)
   ncols=len(Overlaps[0])
   for b in range(nrows):
      for v in range(ncols):
          FindGaps(Overlaps[b][v],gap)
   if fOutLog: print('find gaps of len ',gap,'in Time:',tm.time()-tm_stg,file=fOutLog)

if fOutLog: print('Len Dicts:',len(t_rep), len(t_vir), 'Total Time:', tm.time()-tm_st, sep = sep,file=fOutLog)
else: print('Len Dicts:',len(t_rep), len(t_vir), 'Total Time:', tm.time()-tm_st, sep = sep)
print('CommandLine = """ ',end="")
for x in args:  # first print a header line identifying the data.
   if not x or ' ' in x: print('"',x,'"',sep='',end="  ")
   else: print(x,end="  ")
print('"""')

# If requested, print out the collection of all intersections, to a file or stdout as requested.
if fOutOverlap:
    print('Title="""',file = fOutOverlap,end="")
    for x in args:  # first print a header line identifying the data.
       if not x or ' ' in x: print('"',x,'"',sep='',end="  ",file=fOutOverlap)
       else: print(x,end="  ",file=fOutOverlap)
    print('"""\nOverlaps=\\',file=fOutOverlap)       # next print out the matrix of counts or common substring locations.
    print('[',Overlaps[0],sep='',file = fOutOverlap,end="")
    for i in range(1,len(Overlaps)):
        print(',\\\n ',Overlaps[i],sep='',file = fOutOverlap,end="")
    print(']\n',file = fOutOverlap) 
    if fOutOverlap != sys.stdout: fOutOverlap.close()
#if fOutOverlap:
#    for i in range(len(t_rep)):
#        for j in range(len(t_vir)):
#            print(Overlaps[i][j], sep = sep, file = fOutOverlap, end = "\t")
#        print("",end="\n",file=fOutOverlap)
#    fOutOverlap.close()

if fOutLog: print(t_vir, t_rep,file=fOutLog)
if fOutLog and fOutLog != sys.stdout: fOutLog.close()
if fOutOverlap and fOutOverlap != sys.stdout: fOutOverlap.close()

## Draw the dot plots.

DotHandles=[]  # used to assemble the handles for all the dots on the plot
EndHandles=[]  # used to assemble the handles for all the dots on the plot
def show(): plt.show(block=False)  # shorthand to force a show on a plot
def markersize(size=0,EndSize=0,Handles=DotHandles,EHandles=EndHandles,goshow=True):
   " adjust the markersize of the dots on the plot (default size is 2,6)"
   if not EndSize: EndSize=size+3
   if size==0:
      return ([z.get_markersize() for z in Handles],[z.get_markersize() for z in EHandles])
   else:
      for z in Handles: z.set_markersize(size)
      for z in EHandles: z.set_markersize(EndSize)
      if goshow: show()

def countdots(bact,virus,nrows=1,ncols=1,M=Overlaps,LabelRep=t_rep,LabelVir=t_vir,item=False,goshow=False):
    "count the number of matches, gaps, breaks for forward and reverse"
    Pair=M[bact][virus]
    lens = ((len([t for t in Pair[0] if Pair[0][t]==4]),len([t for t in Pair[0] if Pair[0][t]==3]),len([t for t in Pair[0] if Pair[0][t]==2])),\
            (len([t for t in Pair[1] if Pair[1][t]==4]),len([t for t in Pair[1] if Pair[1][t]==3]),len([t for t in Pair[1] if Pair[1][t]==2])))
    if goshow: print(lens)
    return lens

def PropagateDots(ptlst,direction=1,m=winlen):
    """ ptlst is a list of beginning positions of first unmatched string after a matching string.
    This fcn fills in the rest of the matching string positions from that last matched string on the plot.
    Output: (x,y,xlast,ylast): x,y holds the positions filled in, except that the last position is placed in xlast,ylast.
    """
    x=[]
    y=[]
    xlast=[]
    ylast=[]
    if nopropagate:
       for pt in ptlst:
          #x.extend([pt[0]])
          #y.extend([pt[1]])
          xlast.append(pt[0])
          ylast.append(pt[1])
    else:
       for pt in ptlst:
          x.extend([pt[0]+j for j in range(m-1)])
          y.extend([pt[1]+j*direction for j in range(m-1)])
          xlast.append(pt[0]+m-1)
          ylast.append(pt[1]+(m-1)*direction)
    return (x,y,xlast,ylast)

def drawdots(bact,virus,nrows=1,ncols=1,M=Overlaps,LabelRep=t_rep,LabelVir=t_vir,item=False,goshow=True):
    """
    draw the dot plot
       bact,virus      = index of bacteria & virus in list, respectively
       nrows,ncols     = how many columns/rows to draw (note the role of column row swapped)
       M               = 2d array of intersections
       LabelRep,LabelVir = identifying labels for the bacteria & viruses, resp.
       item            = if needed to override how the plot is placed
       goshow          = set to False if multiple plots must be assembled before showing the result.
    """
    if nrows==0: 
       nrows=len(M)
    #else:
    #   nrows=1
    if ncols==0:
       ncols=len(M[0])
    #else:
    #   ncols=1
    Pair=M[bact][virus]
    if (not loweronly) or (bact<=virus):
       #print('drawdots: ',[bact,virus],t_rep[bact],t_vir[virus]) # tmp
       print('(col,row)',(bact,virus),'(matches,SNPs,breaks): ',end="")
       countdots(bact,virus,goshow=True)
       if noplot: return
       if nrows*ncols > 1:
          if not item: item=1+bact+nrows*virus
          #print('plot: ',nrows,ncols,item,bact,virus,virus==ncols-1,bact==0)
          plt.subplot(ncols,nrows,item)
       if nrows*ncols == 1:
          plt.subplot(1,1,1)
          DotHandles.clear()
          EndHandles.clear()
       if not tightaxes: z=plt.plot([t_rep[bact][2],0],[0,t_vir[virus][2]],'k.',markersize=0)
       #handles.append(z)
       DotHandles.extend(plt.plot([t[0] for t in Pair[0] if Pair[0][t]==4],[t[1] for t in Pair[0] if Pair[0][t]==4],'b.'))
       DotHandles.extend(plt.plot([t[0] for t in Pair[1] if Pair[1][t]==4],[t[1] for t in Pair[1] if Pair[1][t]==4],'r.'))
       if nofill:
          EndHandles.extend(plt.plot([t[0] for t in Pair[0] if Pair[0][t]==3],[t[1] for t in Pair[0] if Pair[0][t]==3],'g.'))
          EndHandles.extend(plt.plot([t[0] for t in Pair[1] if Pair[1][t]==3],[t[1] for t in Pair[1] if Pair[1][t]==3],'g.'))
          EndHandles.extend(plt.plot([t[0] for t in Pair[0] if Pair[0][t]==2],[t[1] for t in Pair[0] if Pair[0][t]==2],'k.'))
          EndHandles.extend(plt.plot([t[0] for t in Pair[1] if Pair[1][t]==2],[t[1] for t in Pair[1] if Pair[1][t]==2],'m.'))
       else:
          (x,y,xlast,ylast) = PropagateDots([t for t in Pair[0] if Pair[0][t]==3])
          DotHandles.extend(plt.plot(x,y,'b.'))
          EndHandles.extend(plt.plot(xlast,ylast,'g.'))
          (x,y,xlast,ylast) = PropagateDots([t for t in Pair[1] if Pair[1][t]==3],direction=-1)
          DotHandles.extend(plt.plot(x,y,'r.'))
          EndHandles.extend(plt.plot(xlast,ylast,'g.'))
          (x,y,xlast,ylast) = PropagateDots([t for t in Pair[0] if Pair[0][t]==2],m=winlen-1)
          DotHandles.extend(plt.plot(x,y,'b.'))
          EndHandles.extend(plt.plot(xlast,ylast,'k.'))
          (x,y,xlast,ylast) = PropagateDots([t for t in Pair[1] if Pair[1][t]==2],direction=-1,m=winlen-1)
          DotHandles.extend(plt.plot(x,y,'r.'))
          EndHandles.extend(plt.plot(xlast,ylast,'m.'))
    plt.grid(not nogrid)
    markersize(marksizes[0],marksizes[1],goshow=False)
    #plt.title(str((bact,virus)))
    if labels>0:
       if virus==ncols-1 or ncols==1: plt.xlabel(str(LabelRep[bact][1])[:labels])
       if bact==0 or nrows==1: plt.ylabel(str(LabelVir[virus][1])[:labels])
    elif labels<0:
       if virus==ncols-1 or ncols==1: plt.xlabel(str(LabelRep[bact][0])[:abs(labels)])
       if bact==0 or nrows==1: plt.ylabel(str(LabelVir[virus][0])[:abs(labels)])
    if goshow: show()
def draw_all(Overlap=Overlaps,goshow=True):
   """
   draw_all() redraws the entire plot.
   Use drawdot(b,v) to draw a single plot.
   You might want to use plt.figure() to open a new figure.
   """
   nrows=len(Overlap)
   ncols=len(Overlap[0])
   DotHandles.clear()
   EndHandles.clear()
   for i in range(nrows):
      for j in range(ncols):
          drawdots(i,j,nrows,ncols,M=Overlap,goshow=False) 
   if goshow and not noplot: show()
if not countonly: draw_all(Overlaps,goshow=False)
if not noplot:
   plt.suptitle('common motifs of length '+str(m)+((noreverse and ' ') or '. red: reverse compl.; ')+'green: SNP')
   if sys.flags.interactive: plt.show(block=False)
   else: plt.show()
