import os
import sys
import lzma
import time
import contextlib
import logging
from tqdm import tqdm
from numpy import log2
from Bio.Seq import Seq
from numpy.random import randint
from collections import Counter
from random import choices

logging.basicConfig(format='SIMGENES :: %(levelname)s - %(message)s')

def entropy(seq):
    if '*' in seq:
        seq = seq.replace('*', '')
    if 'X' in seq:
        seq = seq.replace('X', '')
    df = dict(Counter(seq))
    x = 1 / len(df)
    maxH = len(df)*(x*log2(x))
    s = sum(df.values())
    H = 0
    for _, v in df.items():
        H += (v/s)*log2(v/s)
    if maxH == 0: N=0
    else: N=H/maxH
    return -H, -maxH, N


def mtime(L:int) -> float:
    t = (L-30)/(900-30)
    return (t*60)+0.5


def format_orf(seq: str):
    '''
    Alternate starts: PMID 32209305
    Alternate stops: PMID 25217634
    '''
    starts = ['ATG', 'GTG', 'TTG']
    init_weights = [69.47, 10.40, 5.62]
    stops = ['TAA', 'TGA', 'TAG']
    st_weights = [60, 35, 25]
    seq = str(seq).upper()
    a = choices(starts,
                weights=init_weights,
                k=1)[0]            
    b = choices(stops,
                weights=st_weights,
                k=1)[0]
    return a + seq + b        


def check_internal_stops(seq: str):
    stops = ['TAA', 'TGA', 'TAG']
    idx_list = range(0, len(seq), 3)    
    codons = [seq[idx:idx+3] for idx in idx_list]
    return [idx for idx, x in enumerate(codons) if x in stops]
    
    
def get_prot(seq: str):
    x = 3 - len(seq)%3
    if x != 3:
        seq = seq + 'N'*x
    return str(Seq(seq).translate())
    

def det_props(GC=None):
    x = [isinstance(GC, float),
         isinstance(GC, int),
         GC == None]
    x = sum(x)     
    if x == 0:
        logging.critical('GC content is not set properly')
        sys.exit()
    if (GC != None) and (GC > 100):
        logging.critical('GC content is above 100%')
        sys.exit() 
    if (GC != None) and (abs(GC) < 1) and (GC != 0):
        logging.critical('GC content is under 1%')
        sys.exit()
    def prop_gen(GC):
        if (GC == None) or (GC == 0):
            logging.info('Using randomly selected GC')
            GC = randint(214, 749) / 10  # PMID: 23028785
            GC = float(f'{GC:.1f}')
        if GC == -1:
            a, t, c, g = [25, 25, 25, 25]
        else:
            c = randint(0, int(GC*10)) / 10
            g = GC - c
            a = int((100 - GC)*10)
            a = randint(0, a) / 10
            t = 100 - a - GC
            t = float(f'{t:.1f}')
        if sum([a, t, c, g]) == 100:
            return [a, t, c, g]    
    p = None
    while not p:
        p = prop_gen(GC)
    logging.info(f'[a = {p[0]:.2f}%, t = {p[1]:.2f}%, c = {p[2]:.2f}%, g = {p[3]:.2f}%]')    
    return p
    

def get_len(minL: int, maxL: int):
    L = choices(list(range(minL, maxL+3, 3)),
                k=1)[0]            
    return L


def internal_seq(L, gc):
    maxtime = mtime(L)
    tic = time.perf_counter()
    props = det_props(gc)
    props_edit = props.copy()
    props_edit[0] = sum(props_edit[0:2])
    del props_edit[1:2]
    rand_seq = choices(list('ATCG'),
                       weights=props,
                       k=L-6)
    rand_seq = ''.join(rand_seq)
    p = check_internal_stops(rand_seq)
    for x in p:
        loc = x*3
        newbp = choices('ACG',
                        weights=props_edit,
                        k=1)[0]
        rand_seq = rand_seq[0:loc] + newbp + rand_seq[loc+1:]
    p = check_internal_stops(rand_seq)
    toc = time.perf_counter()
    if len(p) == 0:
        rand_seq = format_orf(rand_seq)    
        return rand_seq
    else:
        logging.error('Did not generate a sequence')
        logging.debug('Could not set a seq without internal stops')
        return None    
    

def random_gene(L=None, gc: int = None, trans: bool = False):
    '''
    You can provide a range to the length of random genes,
    by attributing to L a tuple with min and max 
    gene lengths in base pairs, respectively. For example:
    (30, 300)
    
    You also can provide a fixed gene length by using an integer,
    or yet conform with the randomly selected length from the
    range between 300 to 900 base pairs - considering the 
    average gene length of 819 bp -- Rajic et al 2005
    (default).
    
    gc - is the gc content (0-100) of the desired genes, in case
    it is not given, the program will select randomly the
    proportion of nucleotides between 21.4<=gc<=74.9 (typical
    gc content of ORFs). There is also the option -1, which
    means that all nucleotides are expected at same frequency
    (default)
    '''
    if isinstance(L, int):
        if (L%3) != 0:
            logging.warning(f'Lenght = {L} is not divisible by 3')
            logging.info('Assigning random length')
            L = None
    elif isinstance(L, tuple) and len(L) == 2:
        if L[0] < L[1]:
            L = get_len(L[0], L[1])
        else:
            logging.warning('Reordering your tuple')
            logging.info(f'L = ({L[1]}, {L[0]})')
            L = get_len(L[1], L[0])
    if L == None:
        L = get_len(300, 900)
    if not L:
        logging.critical('Length is not set properly')
        sys.exit()    
    logging.info(f'Selected L={L}')
    rand_seq = internal_seq(L, gc)
    if rand_seq and trans:
        return rand_seq, get_prot(rand_seq)
    elif rand_seq:
        return rand_seq, None
    else:
        return None, None
    

def batch_random(ofile:str='output_testseq',
                 L:(int, tuple)=(30, 300),
                 GC:float=50.8,
                 transl:bool=False,
                 n:int=1_000_000):
    devnull = open(os.devnull, 'w')
    ntout = lzma.open(f'{ofile}.fna.xz', 'wt')
    if transl:
        ptnout = lzma.open(f'{ofile}.faa.xz', 'wt')
    for idx in tqdm(range(n)):
        with contextlib.redirect_stdout(devnull):
            s = None
            while s == None:   
                s, p = random_gene(L, GC, transl)                   
            ntout.write(f'>random_gene_{idx}\n{s}\n')
            if p: ptnout.write(f'>random_gene_{idx}\n{p}\n')
    devnull.close()
    ntout.close()
    if transl: ptnout.close()
    
