
import re, sys, string

###########################
## Regular expression
#Old regex
#G4_bulges = re.compile(r"(?=(G(?:[ATC]{0,6}?G){2}(?:[ATCG]{1,7}?G(?:[ATC]{0,6}?G){2}){3}))", re.I)  # find each motif which overlap
#G4_bulges2 = re.compile(r"G([ATC]{0,6}?)G([ATC]{0,6}?)G[ATCG]{1,7}?G([ATC]{0,6}?)G([ATC]{0,6}?)G[ATCG]{1,7}?G([ATC]{0,6}?)G([ATC]{0,6}?)G[ATCG]{1,7}?G([ATC]{0,6}?)G([ATC]{0,6}?)G", re.I)  # find each motif which overlap
#C4_bulges = re.compile(r"(?=(C(?:[ATG]{0,6}?C){2}(?:[ATCG]{1,7}?C(?:[ATG]{0,6}?C){2}){3}))", re.I)  # find each motif which overlap
#C4_bulges2 = re.compile(r"C([ATG]{0,6}?)C([ATG]{0,6}?)C[ATCG]{1,7}?C([ATG]{0,6}?)C([ATG]{0,6}?)C[ATCG]{1,7}?C([ATG]{0,6}?)C([ATG]{0,6}?)C[ATCG]{1,7}?C([ATG]{0,6}?)C([ATG]{0,6}?)C", re.I)  # find each motif which overlap

#New, truncated G4BS regex: max 0-3 nucleotide insertions and linker length between 1-3 bases
G4_bulges = re.compile(r"(?=(G(?:[ATC]{0,3}?G){2}(?:[ATCG]{1,3}?G(?:[ATC]{0,3}?G){2}){3}))", re.I)  # find each motif which overlap
G4_bulges2 = re.compile(r"G([ATC]{0,3}?)G([ATC]{0,3}?)G[ATCG]{1,3}?G([ATC]{0,3}?)G([ATC]{0,3}?)G[ATCG]{1,3}?G([ATC]{0,3}?)G([ATC]{0,3}?)G[ATCG]{1,3}?G([ATC]{0,3}?)G([ATC]{0,3}?)G", re.I)  # find each motif which overlap
C4_bulges = re.compile(r"(?=(C(?:[ATG]{0,3}?C){2}(?:[ATCG]{1,3}?C(?:[ATG]{0,3}?C){2}){3}))", re.I)  # find each motif which overlap
C4_bulges2 = re.compile(r"C([ATG]{0,3}?)C([ATG]{0,3}?)C[ATCG]{1,3}?C([ATG]{0,3}?)C([ATG]{0,3}?)C[ATCG]{1,3}?C([ATG]{0,3}?)C([ATG]{0,3}?)C[ATCG]{1,3}?C([ATG]{0,3}?)C([ATG]{0,3}?)C", re.I)  # find each motif which overlap

#Needed to filter out G4BS with more than one single "G" linkers
G4_linker = re.compile(r"G[ATC]{0,3}?G[ATC]{0,3}?G([ATCG]{1,3}?)G[ATC]{0,3}?G[ATC]{0,3}?G([ATCG]{1,3}?)G[ATC]{0,3}?G[ATC]{0,3}?G([ATCG]{1,3}?)G[ATC]{0,3}?G[ATC]{0,3}?G", re.I)
C4_linker = re.compile(r"C[ATG]{0,3}?C[ATG]{0,3}?C([ATCG]{1,3}?)C[ATG]{0,3}?C[ATG]{0,3}?C([ATCG]{1,3}?)C[ATG]{0,3}?C[ATG]{0,3}?C([ATCG]{1,3}?)C[ATG]{0,3}?C[ATG]{0,3}?C", re.I)  # find each motif which overlap

#Unchanged for either analysis
Gclusters = re.compile(r"(G{3}).+(G{3})", re.I)
Gclusters2 = re.compile(r"(G{3})", re.I)
Cclusters = re.compile(r"(C{3}).+(C{3})", re.I)
Cclusters2 = re.compile(r"(C{3})", re.I)

#The motif lengths are also changed here. For "I" (inserts) [0, 6] -> [0, 3], for "N" [1-7] -> [1-3]
G4_bulges_each = re.compile(r"""(?P<G11>G)
                                (?P<I11>[ATC]{0,3}?)
                                (?P<G12>G)
                                (?P<I12>[ATC]{0,3}?)
                                (?P<G13>G)
                                (?P<N1>[ATCG]{1,3}?)
                                (?P<G21>G)
                                (?P<I21>[ATC]{0,3}?)
                                (?P<G22>G)
                                (?P<I22>[ATC]{0,3}?)
                                (?P<G23>G)
                                (?P<N2>[ATCG]{1,3}?)
                                (?P<G31>G)
                                (?P<I31>[ATC]{0,3}?)
                                (?P<G32>G)
                                (?P<I32>[ATC]{0,3}?)
                                (?P<G33>G)
                                (?P<N3>[ATCG]{1,3}?)
                                (?P<G41>G)
                                (?P<I41>[ATC]{0,3}?)
                                (?P<G42>G)
                                (?P<I42>[ATC]{0,3}?)
                                (?P<G43>G)""",re.I|re.X)  # find each motif which overlap
C4_bulges_each = re.compile(r"""(?P<G11>C)
                                (?P<I11>[ATG]{0,3}?)
                                (?P<G12>C)
                                (?P<I12>[ATG]{0,3}?)
                                (?P<G13>C)
                                (?P<N1>[ATCG]{1,3}?)
                                (?P<G21>C)
                                (?P<I21>[ATG]{0,3}?)
                                (?P<G22>C)
                                (?P<I22>[ATG]{0,3}?)
                                (?P<G23>C)
                                (?P<N2>[ATCG]{1,3}?)
                                (?P<G31>C)
                                (?P<I31>[ATG]{0,3}?)
                                (?P<G32>C)
                                (?P<I32>[ATG]{0,3}?)
                                (?P<G33>C)
                                (?P<N3>[ATCG]{1,3}?)
                                (?P<G41>C)
                                (?P<I41>[ATG]{0,3}?)
                                (?P<G42>C)
                                (?P<I42>[ATG]{0,3}?)
                                (?P<G43>C)""", re.I|re.X)  # find each motif which overlap


##G4_bulges_for_findall = re.compile(r"G(?:[ATCG]{0,6}?G){2}(?:[ATCG]{1,7}?G(?:[ATCG]{0,6}?G){2}){3,}", re.I) # find all join motifs with not overlap
##G4_2G = re.compile(r"G{2}(?:[ATCG]{1,7}?G{2}){3,}", re.I)

#Also changed the regex for generic G4 detection using the same numbers as above (([ATCG]{1,7} -> ([ATCG]{1,3})
Generic_G4_g = re.compile(r'G{3,}([ATCG]{1,3}G{3,}){3,}', re.I)
Generic_G4_c = re.compile(r'C{3,}([ATCG]{1,3}C{3,}){3,}', re.I)

##
############################

chromDict = {'chr1':'01','chr2':'02','chr3':'03','chr4':'04','chr5':'05','chr6':'06','chr7':'07','chr8':'08',
             'chr9':'09','chr10':'10','chr11':'11','chr12':'12','chr13':'13','chr14':'14','chr15':'15',
             'chr16':'16','chr17':'17','chr18':'18','chr19':'19','chr20':'20','chr21':'21','chr22':'22',
             'chrX':'23','chrY':'24'}


def main():
    #####################
    ## Parameter

    import argparse
    parser = argparse.ArgumentParser(description='The program use to search G4-bulge in sequence (fasta file)')
    parser.add_argument('fasta_file', type=argparse.FileType('r'), help='Fasta file')
    parser.add_argument('-g', '--number_of_g_cluster', type=int, action='store',
                        dest='number_cluster', default=int(1),
                        help='number of G-cluster without bulge (GGG) (1-4) (defualt: 1)')
    parser.add_argument('-b', '--max_number_bulge', type=int, action='store',
                        dest='max_number_bulge', default=int(4),
                        help='maximum number of bulge (1-8) (defualt: 4)')
    args = parser.parse_args()


    number_cluster = args.number_cluster
    if number_cluster < 0:
        number_cluster = 1
    max_number_bulge = args.max_number_bulge
    
##    extract_Gclusters = re.compile(r"G{3}(?:.+?G{3}){%s}" %number_cluster, re.I)
##    extract_Cclusters = re.compile(r"C{3}(?:.+?C{3}){%s}" %number_cluster, re.I)

    ##
    #####################

    from Bio import SeqIO
    from Bio.Seq import Seq
    print "\t".join(['#SeqName','Start','End','Sequence','ID','Strand','num_Gcluster_no_bulge','num_bulges','label','I1,I2,N1,I3,I4,N2,I5,I6,N3,I7,I8','I1,I2,N1,I3,I4,N2,I5,I6,N3,I7,I8'])
    for seq_record in SeqIO.parse(args.fasta_file, "fasta"):
        name = seq_record.id
        matches = G4_bulges.finditer(str(seq_record.seq))

        #matches = G4_bulges.finditer(s)
        for match in matches:
            test = "Go"
            G_linkers = G4_linker.finditer(match.group(1))
            for link in G_linkers:
                if link.groups().count("G") > 1 or link.groups().count("g") > 1:
                    test = "Fail"
                    #print("next")
                #else:
                    #print("cont")
            #print(test)
            if test == "Fail":
                continue
            for bulges in G4_bulges2.finditer(match.group(1)):
                num_Gcluster_no_bulge = len(filter(lambda x: x[0] == '' and x[1] == '', [[bulges.groups()[i],bulges.groups()[i+1]] for i in range(0,len(bulges.groups()),2)]))
                num_bulges = len(filter(lambda x: x != "", bulges.groups()))
                if num_bulges <= max_number_bulge and num_Gcluster_no_bulge >= number_cluster and not Generic_G4_g.search(match.group(1)) and not C_filter(match.group(1)):
                    if chromDict.has_key(name):
                        g4id = "G4B%s%s" %(chromDict[name],str(match.start(1)).zfill(9))
                    else:
                        g4id = "G4B%s%s" %(name,str(match.start(1)).zfill(9))
                    G4_bulges_detail = G4_bulges_each.search(match.group(1))
                    sequence = match.group(1)
                    Linkers = ['I11','I12','N1','I21','I22','N2','I31','I32','N3','I41','I42']
                    for link in Linkers:
                        startLink,endLink = getPositionG4BS(G4_bulges_detail,link)
                        sequence = maskStringLowercase(sequence,startLink,endLink)
                    for G in ['N3','N2','N1']:
                        startLink,endLink = getPositionG4BS(G4_bulges_detail,G)
                        sequence = sequence[0:startLink]+"-"+sequence[startLink:endLink].lower()+"-"+sequence[endLink:]
                    print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(name, match.start(1), match.end(1), \
                                                             match.group(1).upper(), g4id, "+",num_Gcluster_no_bulge, num_bulges, \
                                                             sequence, ",".join(map(str,[ G4_bulges_detail.group(x) for x in Linkers])),
                                                                         ",".join(map(str,[ len(G4_bulges_detail.group(x)) for x in Linkers])) )
        matches = C4_bulges.finditer(str(seq_record.seq))
        for match in matches:
            test = "Go"
            C_linkers = C4_linker.finditer(match.group(1))
            for link in C_linkers:
                if link.groups().count("C") > 1 or link.groups().count("c") > 1:
                    test = "Fail"
                    # print("next")
                # else:
                #   print("cont")
            # print(test)
            if test == "Fail":
                continue
            for bulges in C4_bulges2.finditer(match.group(1)):
                num_Gcluster_no_bulge = len(filter(lambda x: x[0] == '' and x[1] == '', [[bulges.groups()[i],bulges.groups()[i+1]] for i in range(0,len(bulges.groups()),2)]))
                num_bulges = len(filter(lambda x: x != "", bulges.groups()))
                if num_bulges <= max_number_bulge and num_Gcluster_no_bulge >= number_cluster and not Generic_G4_c.search(match.group(1)) and not G_filter(match.group(1)):
                    if chromDict.has_key(name):
                        g4id = "G4B%s%s" %(chromDict[name],str(match.start(1)).zfill(9))
                    else:
                        g4id = "G4B%s%s" %(name,str(match.start(1)).zfill(9))
                    G4_bulges_detail = C4_bulges_each.search(match.group(1))
                    sequence = match.group(1)
                    Linkers = ['I11','I12','N1','I21','I22','N2','I31','I32','N3','I41','I42']
                    for link in Linkers:
                        startLink,endLink = getPositionG4BS(G4_bulges_detail,link)
                        sequence = maskStringLowercase(sequence,startLink,endLink)
                    for G in ['N3','N2','N1']:
                        startLink,endLink = getPositionG4BS(G4_bulges_detail,G)
                        sequence = sequence[0:startLink]+"-"+sequence[startLink:endLink].lower()+"-"+sequence[endLink:]
                    linker1 = [ reverse_complement_dna(G4_bulges_detail.group(x).lower()) for x in Linkers ]
                    linker1.reverse()
                    linker2 = [ len(G4_bulges_detail.group(x)) for x in Linkers]
                    linker2.reverse()
                    print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(name, match.start(1), match.end(1), \
                                                             reverse_complement_dna(match.group(1).upper()), g4id, "-",num_Gcluster_no_bulge, num_bulges, \
                                                             reverse_complement_dna(sequence), ",".join(map(str,linker1)),
                                                                         ",".join(map(str,linker2)) )
    args.fasta_file.close()       


def getPositionG4BS(match, name):
    g=['G11','I11','G12','I12','G13','N1','G21','I21','G22','I22','G23','N2','G31','I31','G32','I32','G33','N3','G41','I41','G42','I42','G43']
    start = 0
    for name2 in g:
        if name == name2:
            end = start+len(match.group(name2))
            return start, end
        else:
            start += len(match.group(name2))

def maskStringLowercase(string,start,end):
    return string[0:start]+string[start:end].lower()+string[end:]

def match(seq):
    m = G4_bulges.search(seq)
    if m:
        print "yes", m.group(0)
    else:
        print "no"
def reverse_complement_dna(seq):
    return seq.translate(string.maketrans("ATCGatcg", "TAGCtagc"))[::-1]


def C_filter(str):
    # Compile the ReGex
    p = re.compile(r"([cC])\1+", re.I)

    # If the string is empty
    # return false
    if (str == None):
        return False

    # Return if the string
    # matched the ReGex
    if (re.search(p, str)):
        return True
    else:
        return False

def G_filter(str):
    # Compile the ReGex
    p = re.compile(r"([gG])\1+", re.I)

    # If the string is empty
    # return false
    if (str == None):
        return False

    # Return if the string
    # matched the ReGex
    if (re.search(p, str)):
        return True
    else:
        return False

##inf = open("example_input.txt")
##for l in inf.xreadlines():
##    L = l.strip().split("\t")
##    print L[0], L[1], 
##    m = G4_bulges.search(L[1])
##    if m and Gclusters.search(L[1]):
##        print " -> found -> ", m.group(1)
##    else:
##        print " -> not found"
##inf.close()
##

def example1():
    number_cluster = 1
    max_number_bulge = 4
##    extract_Gclusters = re.compile(r"G{3}(?:.+?G{3}){%s}" %number_cluster, re.I)
##    extract_Cclusters = re.compile(r"C{3}(?:.+?C{3}){%s}" %number_cluster, re.I)
    
    name = "test"
    s = "TTTTTTTTGCAAAAAGAAAAAAGTTTTTTTGACAAAAGAAAAAAGTTTTTTTGAACAAAGAAAAAAGTTTTTTTGAAACAAGAAAAAAGTTTTTTT\
    GAAAACAGAAAAAAGTTTTTTTGAAAAACGAAAAAAGTTTTTTTGCCCCCCGCCCCCCGtttttttGCAAAAAGAAAAAAGTTTTTTTGACAAAAGAAAAAAG\
    TTTTTTTGAACAAAGAAAAAAGTTTTTTTGAAACAAGAAAAAAGTTTTTTTGAAAACAGAAAAAAGTTTTTTTGAAAAACGAAAAAAGTTTTTTTTTT"
    s = "GGGGGGGGGGGGGGGGGAAAAAAAGGGGGGGGGGGGGGGGAGGGGGGGGGGGGGGGGAAGGGGGGGGGGGAGGGGGGGGGGGGG"
    s = "GTTTTTGTGGGTGTGGTGTGGTGGGGGG"
    s = "GTGGAAGGGTTGAGGCGGG"
    s = "GTGGGGAGGTGAGGGGGGGTGGAAGGGTTGAGGCGGG"
    matches = G4_bulges.finditer(s)
    for match in matches:
    ##    print name, match.start(1), match.end(1), reverse_complement_dna(match.group(1)).upper()
    ##    print name, match.start(1), match.end(1), match.group(1).upper()
        print(name, match.start(1), match.end(1), match.group(1).upper())
        asd = G4_linker.finditer(match.group(1))
        for a in asd:
            print(a.groups())
            if a.groups().count("G") > 1:
                print("Skipped")
                continue
            else:
                print("Passed")
    ##    new_s = match.group(1).upper()
        for match2 in G4_bulges2.finditer(match.group(1)):
            #if match2.groups().count("g") > 1 or match2.groups().count("G") > 1:
            #    print("KEK")
            #    continue
            #else: print("WUT")
            print "  ", match2.groups(), len(filter(lambda x: x != "", match2.groups()))
            num_Gcluster_no_bulge =  len(filter(lambda x: x[0] == '' and x[1] == '', [[match2.groups()[i],match2.groups()[i+1]] for i in range(0,len(match2.groups()),2)]))
            print num_Gcluster_no_bulge

        G4_bulges_detail = G4_bulges_each.search(match.group(1))
        sequence = match.group(1)
        for link in ['I11','I12','N1','I21','I22','N2','I31','I32','N3','I41','I42']:
            startLink,endLink = getPositionG4BS(G4_bulges_detail,link)
            sequence = maskStringLowercase(sequence,startLink,endLink)
            #print(sequence)
        for G in ['N3','N2','N1']:
            startLink,endLink = getPositionG4BS(G4_bulges_detail,G)
            sequence = sequence[0:startLink]+"-"+sequence[startLink:endLink].lower()+"-"+sequence[endLink:]
            #print(sequence)
        print(sequence)
        
        print("G%sG%sG-%s-G%sG%sG-%s-G%sG%sG-%s-G%sG%sG"%("I%s"%(len(G4_bulges_detail.group('I11'))),"I%s"%(len(G4_bulges_detail.group('I12'))),
                                                    "N%s"%(len(G4_bulges_detail.group('N1'))),
                                                    "I%s"%(len(G4_bulges_detail.group('I21'))),"I%s"%(len(G4_bulges_detail.group('I22'))),
                                                    "N%s"%(len(G4_bulges_detail.group('N2'))),
                                                    "I%s"%(len(G4_bulges_detail.group('I31'))),"I%s"%(len(G4_bulges_detail.group('I32'))),
                                                    "N%s"%(len(G4_bulges_detail.group('N3'))),
                                                    "I%s"%(len(G4_bulges_detail.group('I41'))),"I%s"%(len(G4_bulges_detail.group('I42')))))

##        for match3 in Gclusters.finditer(s):
##            print "  ", match3.groups(), match3.group(1)

def example2():

    #### EXAMPLE ######
    s = "TTTTTTTTGCAAAAAGAAAAAAGTTTTTTTGACAAAAGAAAAAAGTTTTTTTGAACAAAGAAAAAAGTTTTTTTGAAACAAGAAAAAAGTTTTTTT\
    GAAAACAGAAAAAAGTTTTTTTGAAAAACGAAAAAAGTTTTTTTGCCCCCCGCCCCCCGTTTTTTTTGCAAAAAGAAAAAAGTTTTTTTGACAAAAGAAAAAAG\
    TTTTTTTGAACAAAGAAAAAAGTTTTTTTGAAACAAGAAAAAAGTTTTTTTGAAAACAGAAAAAAGTTTTTTTGAAAAACGAAAAAAGTTTTTTTTTT"

    s = "AAAAAAGGGAAGGGAGGGCGGGAAAAAGGCCGAAGGGAAGGGTGGGTTGTGGTGGAGTGGCGTTAAAAGGCCGTGGGTGGGTGGG"
    matches = re.finditer(r"(?=(G(?:[ATCG]{0,6}?G){2}(?:[ATCG]{1,7}?G(?:[ATCG]{0,6}?G){2}){3}))",s)
    matches = G4_bulges.finditer(s)
    for match in matches:
        if match not in c_filter.finditer(match.group(1)):
            print(match.group(1))
    print('--------------')
    matches = re.findall(r"G(?:[ATCG]{0,6}?G){2}(?:[ATCG]{1,7}?G(?:[ATCG]{0,6}?G){2}){3,}",s)
    for match in matches:
        print(match)

if __name__ == "__main__":
    main()
##    example1()
#G4_bulges = re.compile(r"(?=(G(?:[ATC]{0,3}?G){2}(?:[ATCG]{1,3}?G(?:[ATC]{0,3}?G){2}){3}))", re.I)
#s = "AAAAAAGGGAAGGGAGGGCGGGAAAAAGGCCGAAGGGAAGGGTGGGTTGTGGTGGAGTGGCGTTAAAAGGCCGTGGGTGGGTGGG"
#matches = G4_bulges.finditer(s)

#for match in matches:
#    if not isValidGUID(match.group(1)):
#        print(match.group(1))