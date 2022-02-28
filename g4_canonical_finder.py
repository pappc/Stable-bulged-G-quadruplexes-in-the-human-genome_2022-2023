
#import re, sys, string
import regex as re, sys, string
###########################

#Generic_G4_g = re.compile(r'G{3,}([ATCG]{1,3}G{3,}){3,}', re.I)
#Generic_G4_c = re.compile(r'C{3,}([ATCG]{1,3}C{3,}){3,}', re.I)
Generic_G4_g = re.compile(r'(?=(G{3,}[ATCG]{1,7}?G{3,}[ATCG]{1,7}?G{3,}[ATCG]{1,7}?G{3,}))', re.I)
Generic_G4_c = re.compile(r'(?=(C{3,}[ATCG]{1,7}?C{3,}[ATCG]{1,7}?C{3,}[ATCG]{1,7}?C{3,}))', re.I)
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
    args = parser.parse_args()

    ##
    ##################
    ###

    from Bio import SeqIO
    from Bio.Seq import Seq
    print "\t".join(['#SeqName','Start','End','Sequence','ID','Strand','num_Gcluster_no_bulge','num_bulges','label','I1,I2,N1,I3,I4,N2,I5,I6,N3,I7,I8','I1,I2,N1,I3,I4,N2,I5,I6,N3,I7,I8'])
    for seq_record in SeqIO.parse(args.fasta_file, "fasta"):
        name = seq_record.id
        matches = Generic_G4_g.finditer(str(seq_record.seq), overlapped = True)
        
        for match in matches:
            if chromDict.has_key(name):
                g4id = "G4B%s%s" %(chromDict[name],str(match.start(1)).zfill(9))
            else:
                g4id = "G4B%s%s" %(name,str(match.start(1)).zfill(9))

            print "%s\t%s\t%s\t%s\t%s\t%s" %(name, match.start(1), match.end(1), match.group(1).upper(), g4id, "+")

        matches = Generic_G4_c.finditer(str(seq_record.seq), overlapped = True)
        for match in matches:
            if chromDict.has_key(name):
                g4id = "G4B%s%s" %(chromDict[name],str(match.start(1)).zfill(9))
            else:
                g4id = "G4B%s%s" %(name,str(match.start(1)).zfill(9))
                sequence = match.group(1)
            print "%s\t%s\t%s\t%s\t%s\t%s" %(name, match.start(1), match.end(1), reverse_complement_dna(match.group(1).upper()), g4id, "-")
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
