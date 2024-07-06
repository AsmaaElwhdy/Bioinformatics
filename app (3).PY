import streamlit as st
import numpy as np
import pandas as pd
st.set_page_config(page_title="Bioinformatics")



file=open("dna1.fasta")
l=[i for i in file]
seq=l[1][0:-1]

def GC_Content(seq):
    l=len(seq)
    num_G=seq.count("G")
    num_C=seq.count("C")
    total=num_C+num_G
    return total/l

def AT_Content(seq):
    l=len(seq)
    num_A=seq.count("A")
    num_T=seq.count("T")
    total=num_A+num_T
    return total/l

def Complement(seq):
    dic={"G":"C","C":"G","A":"T","T":"A"}
    s=list(seq)
    for i in range(len(seq)):
        s[i]=str(dic[s[i]])
    s="".join(s)
    return s

def Reverse(seq):
    s=list(seq)
    s=reversed(s)
    s="".join(s)
    return s
    
def Reverse_Complement(seq):
    seq=Reverse(seq)
    seq=Complement(seq)
    return seq 

#########################################################################################

def Translation_Table(seq):
    dic = {"TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V",
           "TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V",
           "TTA" : "L", "CTA" : "L", "ATA" : "I", "GTA" : "V",
           "TTG" : "L", "CTG" : "L", "ATG" : "M", "GTG" : "V",
           "TCT" : "S", "CCT" : "P", "ACT" : "T", "GCT" : "A",
           "TCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
           "TCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
           "TCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
           "TAT" : "Y", "CAT" : "H", "AAT" : "N", "GAT" : "D",
           "TAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
           "TAA" : "*", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
           "TAG" : "*", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
           "TGT" : "C", "CGT" : "R", "AGT" : "S", "GGT" : "G",
           "TGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
           "TGA" : "*", "CGA" : "R", "AGA" : "R", "GGA" : "G",
           "TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" 
           }
    s=""
    sf=""
    flag=0
    for i in range(0,len(seq)-2,3):
        if dic[seq[i:i+3]]=="M":
            flag=1
        elif (dic[seq[i:i+3]]=="*"):
            flag=0
        if flag==1:
            s+=dic[seq[i:i+3]]
        sf+=dic[seq[i:i+3]]
    return sf,s        
#########################################################################################

def match(seq,sub_seq):
    x=-1
    for i in range(len(seq)-len(sub_seq)+1):
        if sub_seq==seq[i:i+len(sub_seq)]:
            x=i
    return x

# print(match('ACGT','dsTgasd'))
######################################################################################
def Badchars(seq,sub_seq):
    table=np.zeros([4,len(sub_seq)])     
    row=["C","G","A","T"]
    for i in range (4):
        num=-1
        for j in range (len(sub_seq)):
            if row[i]==sub_seq[j]:
                table[i,j]=-1
                num=-1
            else:
                num+=1
                table[i,j]=num
    x=-1
    i=0
    while(i<len(seq)-len(sub_seq)+1):
        if sub_seq==seq[i:i+len(sub_seq)]:
            x=i
        else:
            for j in range(i+len(sub_seq)-1,i-1,-1):
                if seq[j] != sub_seq[int(j-i)]:
                    k=row.index(seq[j])
                    i+=table[k,j-i]
                    break
        i=int(i+1)
    return x

# print(Badchars('ACGT','C'))


file=open("dna1.fasta")
l=[i for i in file]
p=l[1][0:-1]

#########################################################################################
def preprocess_good_suffix(sub_seq):
    m = len(sub_seq)
    suffixes = [0] * m
    border = m

    for i in range(m - 1, -1, -1):
        if is_suffix(sub_seq, i):
            border = i
        suffixes[m - 1 - i] = border - i

    return suffixes


def is_suffix(sub_seq, i):
    m = len(sub_seq)
    j = m - 1
    while i >= 0 and sub_seq[i] == sub_seq[j]:
        i -= 1
        j -= 1
    return i == -1


def good_suffix_algorithm(seq, sub_seq):
    n = len(seq)
    m = len(sub_seq)
    suffixes = preprocess_good_suffix(sub_seq)
    skip = [m] * (m + 1)

    for i in range(m):
        skip[m - suffixes[i]] = m - i - 1

    i = 0
    while i <= n - m:
        j = m - 1

        while j >= 0 and sub_seq[j] == seq[i + j]:
            j -= 1

        if j < 0:
            # Pattern found at index i
            print("Pattern found at index", i)

        i += max(skip[j + 1], m - j)
##########################################################################

import bisect

def IndexSorted(seq,ln):
    index = []
    for i in range(len(seq)-ln+1):
        index.append((seq[i:i+ln], i))
    index.sort() 
    return index
# print(IndexSorted('ACGT',3))

def query(t,p,index):
    keys = [r[0] for r in index]
    st = bisect.bisect_left(keys,p[:len(keys[0])])
    
    en = bisect.bisect(keys,p[:len(keys[0])])
    hits = index[st:en] 
    print(hits)
    l=[h[1] for h in hits ]
    offsets=[]
    for i in l:
        if t[i:i+len(p)]==p:
            offsets.append(i)
    return offsets

file=open("dna1.fasta")
l=[i for i in file]
seq=l[1][0:-1]
# p="AAG"
ln=3

index=IndexSorted(seq,ln)
# print(index)
# print(query(t,p,index))

################################################################################

from itertools import permutations

def overlap(a,b,min_length=3):
    start=0
    while True:
        start=a.find(b[:min_length],start)
        if start==-1:
            return 0
        if b.startswith(a[start:]):
            return len(a[start:])
        start+=1
        
        
print(overlap("ACGGTAGGT", "GGTAGGTCC",3))
a="ACGGTA" 
b="GTAGTA"
min_len=3
j=len(a)-1
count=0
for i in range(len(b)-1,-1,-1):
    if a[j]==b[i]:
        count+=1
        j-=1
    else:
        count=0
        j=len(a)-1
        if a[j]==b[i]:
            count+=1
            j-=1
if count>=min_len:
    print('overlap:',count)
else:
    print("no overlap")

        
def native_overlap(reads,k):
    olap={}
    for a,b in permutations(reads,2):
        olen=overlap(a, b,k)
        if olen>0:
            olap[(a,b)]=olen
    return olap

print(native_overlap(["ACGGTA", "GGTACC","GCATACG"],3))

#######################################################################
def compute_lps_array(pattern):
    m = len(pattern)
    lps = [0] * m
    length = 0
    i = 1

    while i < m:
        if pattern[i] == pattern[length]:
            length += 1
            lps[i] = length
            i += 1
        else:
            if length != 0:
                length = lps[length - 1]
            else:
                lps[i] = 0
                i += 1

    return lps

def kmp_search(text, pattern):
    n = len(text)
    m = len(pattern)
    lps = compute_lps_array(pattern)

    i = j = 0
    indices = []

    while i < n:
        if pattern[j] == text[i]:
            i += 1
            j += 1

        if j == m:
            indices.append(i - j)
            j = lps[j - 1]

        elif i < n and pattern[j] != text[i]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1

    return indices
#################################################################
#approximate
def levenshtein_distance(sub_seq, seq):
    m, n = len(sub_seq), len(seq)
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0:
                dp[i][j] = j
            elif j == 0:
                dp[i][j] = i
            elif sub_seq[i-1] == seq[j-1]:
                dp[i][j] = dp[i-1][j-1]
            else:
                dp[i][j] = 1 + min(dp[i-1][j],      # Deletion
                                  dp[i][j-1],      # Insertion
                                  dp[i-1][j-1])    # Substitution

    return dp[m][n]
#############################################################################
#kmer
def create_kmer_index(seq, k):
    kmer_index = {}
    n = len(seq)

    for i in range(n - k + 1):
        kmer = seq[i:i+k]
        if kmer not in kmer_index:
            kmer_index[kmer] = [i]
        else:
            kmer_index[kmer].append(i)

    return kmer_index
###################################################################################

def main():
    menu=['Home', 'Show data', 'Reverse Complement', 'protein','Good Suffix', 'Matching', 'Indexed Sorted', 'Suffix Array', 'KMP','Approximate Matching','K-mer']
    st.sidebar.write('# Select algorithm')
    choice=st.sidebar.selectbox('',menu)
    if choice=='Home':
        st.title('Bioinformatics Project')

        st.markdown("""\n 
        #### Bioinformatics is an interdisciplinary field that combines biology,
        #### computer science, and statistics to analyze and interpret biological 
        #### data .
          """)
        
        expander=st.expander('Team Members ')
        expander.markdown("""
            **Nancy Wael**\n
            **Asmaa Mohamed**\n
            **Esraa Mohamed**\n
            **Ahmed Hussein**\n
            **Ahmed Hassan**\n
            **Abdalla Ashraf**\n
            
        

        """)

    global seq, p, t
    if choice=='Show data':
        st.subheader('DNA Fasta\n')
        infile = open('HAPPENN_dataset.fasta')
        tb=[]
        for line in infile:
            if line[0]==">":
                s=line.split("|lcl|")
            else:
                if s[3]=='non-hemolytic' or s[3]=='non-hemolytic\n':
                    tb.append([line[:-1],0])
                else:
                    tb.append([line[:-1],1])
        # print(tb)
        if st.button('Show sequence'):
            st.write(tb)
        head=['Sequence','y']
        df=pd.DataFrame(tb,columns=head)
        df.to_csv("HAPPENN.csv", index=False)
        if st.button('Show Data Frame'):
            st.dataframe(df)

    if choice=='Reverse Complement':
        st.subheader('Reverse Complement\n')
        if st.button('Show Data And Predict'):
        
            st.write(seq)
            st.success(f"GC Content Is: {GC_Content(seq)}")
            st.success(f"AT Content Is: {AT_Content(seq)}")
            st.success(f"Complement Is: {Complement(seq)}")
            st.success(f"Reverse Is: {Reverse(seq)}")
            st.success(f"Reverse Complement Is: {Reverse_Complement(seq)}")
    
    if choice=='protein':
        st.subheader('DNA TO PROTEIN\n')
        if st.button('Show Data And Predict'):
        
            st.write(seq)
            st.success(f"Amino Asid series Is: {Translation_Table(seq)}")

    if choice=='Matching':
        st.subheader('Matching\n')
        t= st.text_input('Enter String..')
        if st.button('Show Data And Predict'):
        
            st.write(seq)
            st.success(f" Exact match Is: {match(p,t )}")  
            st.success(f" Badcharcter Is: {Badchars(p,t )}")       
   

    if choice== 'Good Suffix':
        st.subheader('Good Suffix\n')
        p = st.text_input("Enter Pattern:")
        t = st.text_input("Enter Text:")

        if st.button("Run Algorithm"):
            good_suffix_algorithm(t, p)





    if choice=='Indexed Sorted':
            st.subheader('Indexed Sorted\n')
            p= st.text_input('Enter String..')
            if st.button('Show Data And Predict'):
                st.write(seq)
                st.success(f" Index Sorted Is: {IndexSorted(seq,ln)}")  
                st.success(f" query Is: {query(seq,p,index)}")

    
                    


    if choice=='Suffix Array':
                st.subheader('Suffix Array\n')
                # T = 'ACGACTACGATAAC$'
                T= st.text_input('Enter String..')
                l=[]
                for i in range(len(T)):
                    l.append(T[i:])
                l2=l.copy()
                l.sort()
                dec={}
                for i in range(len(l)):
                    dec[l[i]]=i
                table=[]
                for i in range(len(l)):
                    table.append([l2[i],i,dec[l2[i]]])
                # print(table,'\n')
                # print(l)
                if st.button('Show Table'):
                    
                    st.success(f" Table Is: {table}")

   
    if choice=='KMP':
            st.subheader('KMP\n')
            t = st.text_input("Enter Text:")
            p = st.text_input("Enter Pattern:")
            if st.button("Run KMP Algorithm"):
               if t and p:
                  result = kmp_search(t, p)
                  if result:
                    st.success(f"Pattern found at indices: {result}")
                  else:
                    st.info("Pattern not found in the given text.")
   


    if choice=='Approximate Matching':
             st.subheader('Approximate Matching\n')
             t = st.text_input("Enter Text:")
             p = st.text_input("Enter Pattrn:")
             
             if st.button("Calculate Levenshtein Distance"):
                  distance = levenshtein_distance(p, t)
                  st.success(f"Levenshtein Distance: {distance}")
    
    if choice=='K-mer':
             st.subheader('K-mer\n')
             seq = st.text_input("Enter Text:")
             k = st.number_input("Enter k:", min_value=1, step=1)
             if st.button("Create k-mer Index"):
                 kmer_index = create_kmer_index(seq, int(k))
                 st.success("k-mer Index created successfully.")
                 st.write("k-mer Index:")
                 st.write(kmer_index)

if __name__ == "__main__":
    main()