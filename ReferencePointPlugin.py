import sys
import PyPluMA

class ReferencePointPlugin:
    def input(self, outfile):
       paramfile = open(outfile, 'r')
       self.parameters = dict()
       for line in paramfile:
         contents = line.split('\t')
         self.parameters[contents[0]] = contents[1].strip()

       self.csvfile = open(PyPluMA.prefix()+"/"+self.parameters["csvfile"], 'r')
       self.startpos = int(self.parameters["startpos"])
       self.refgene = self.parameters["refgene"]
       self.refstrand = self.parameters["refstrand"]
       self.genomestart = int(self.parameters["genomestart"])
       self.genomeend = int(self.parameters["genomeend"])


    def run(self):
       #self.header = ["locus_tag", "type", "+/-", "start", "end", "wstart", "wend", "inference", "note", "codon_start", "transl_table", "pseudo", "product", "protein_id", "translation", "gene", "EC_number", "anticodon", "ncRNA_class", "db_xref", "transl_except"]

       self.header = self.csvfile.readline().strip().split(',')

       self.entries = []
       for line in self.csvfile:
           self.entries.append(line.strip().split(','))
 
       gene_idx = self.header.index("gene")
       start_idx = self.header.index("start")
       end_idx = self.header.index("end")
       wstart_idx = self.header.index("wstart")
       wend_idx = self.header.index("wend")
       plus_minus_idx = self.header.index("+/-")

       refgene_idx = -1
       for i in range(len(self.entries)):
           if (self.entries[i][gene_idx] == self.refgene):
              refgene_idx = i
       refgene = self.entries[refgene_idx]
       print(refgene)
       refplusminus = refgene[plus_minus_idx]
       refstart = int(refgene[start_idx])
       refend = int(refgene[end_idx])
       if (refplusminus == self.refstrand):
           reverseStrands = False
       else:
           reverseStrands = True
           print("WARNING REVERSING")

       self.outputentries = []
       if (not reverseStrands): # Easier, add a delta to start and end (watch wraparound)
           delta = int(self.startpos - refstart)
           for i in list(range(refgene_idx, len(self.entries)))+list(range(0,refgene_idx)):
               entry = self.entries[i]
               # Initially, add delta
               # Then check boundaries
               entry[start_idx] = str(int(entry[start_idx])+delta)
               if (len(entry[wstart_idx]) != 0):
                   entry[end_idx] = str(int(entry[wend_idx])+delta)
                   entry[wstart_idx] = ""
                   entry[wend_idx] = ""
               else:
                   entry[end_idx] = str(int(entry[end_idx])+delta)
               if (int(entry[start_idx]) <= 0):  # Out of bounds, negative
                   entry[start_idx] = str(int(entry[start_idx])+self.genomeend)
                   entry[end_idx] = str(int(entry[start_idx])+self.genomeend)
               elif (int(entry[start_idx]) > self.genomeend):
                   entry[start_idx] = str(int(entry[start_idx])-self.genomeend)
                   entry[end_idx] = str(int(entry[end_idx])-self.genomeend)
               self.outputentries.append(entry)
       else: # Tougher, need to reverse direction and index
           donealready = []
           #print(self.entries)
           #print(len(self.entries))
           for i in list(range(refgene_idx, -1, -1))+list(range(len(self.entries)-1,refgene_idx,-1)):
               entry = self.entries[i]
               # Initially, add delta
               # Then check boundaries
               if (entry[plus_minus_idx] == "+"):
                  entry[plus_minus_idx] = "-"
               else:
                  entry[plus_minus_idx] = "+"
               myEnd = entry[end_idx]
               entry[end_idx] = str(self.startpos+(refend-int(entry[start_idx])))
               if (len(entry[wstart_idx]) != 0):
                   print("FOUND A WRAP")
                   print(entry)
                   myStart = entry[start_idx]
                   entry[start_idx] = str(self.startpos+(refend-int(entry[wend_idx])))
                   entry[end_idx] = str(int(entry[start_idx])+int(entry[wend_idx])+int(myEnd)-int(myStart))
                   entry[wstart_idx] = ""
                   entry[wend_idx] = ""
               else:
                   entry[start_idx] = str(self.startpos+(refend-int(myEnd)))
               if (int(entry[start_idx]) <= 0):  # Out of bounds, negative
                   entry[start_idx] = str(int(entry[start_idx])+self.genomeend)
                   entry[end_idx] = str(int(entry[end_idx])+self.genomeend)
               elif (int(entry[start_idx]) > self.genomeend):
                   entry[start_idx] = str(int(entry[start_idx])-self.genomeend)
                   entry[end_idx] = str(int(entry[end_idx])-self.genomeend)
               self.outputentries.append(entry)


    def output(self, filename):
       outfile = open(filename, 'w')
       for i in range(len(self.header)):
           outfile.write(self.header[i])
           if (i == len(self.header)-1):
               outfile.write("\n")
           else:
               outfile.write(',')

       for j in range(len(self.outputentries)):
           for i in range(len(self.outputentries[j])): # Don't do locus_tag
               outfile.write(self.outputentries[j][i])
               if (i == len(self.outputentries[j])-1):
                   outfile.write("\n")
               else:
                   outfile.write(',')


