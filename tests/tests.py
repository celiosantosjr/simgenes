def test_entropy():
    from simgenes.simgenes import entropy
    res = {'ATGCTGAGATGAGAGATTT': (1.7990156706551217, 2.0, 0.8995078353275608),
           'AAAAAAAAAAAAAAAAAAA': (-0.0, -0.0, 0),
           'MAKLKHIPQTTTTAHH': (2.952819531114783, 3.169925001442312, 0.9315108495536184),
           'MAKLKHIPQ*TTTTAHH*': (2.952819531114783, 3.169925001442312, 0.9315108495536184)}
    for k, (v1, v2, v3) in res.items():
        r1, r2, r3 = entropy(k)
        assert (r1, r2, r3) == (v1, v2, v3)


def test_format_orf():
    from simgenes.simgenes import format_orf
    seq = 'XXXXXX'
    seq = format_orf(seq)
    a = seq[0:3]
    b = seq[-3:]
    starts = ['ATG', 'GTG', 'TTG']
    stops = ['TAA', 'TGA', 'TAG']
    assert a in starts
    assert b in stops
    assert len(seq) % 3 == 0


def test_internal_stops():
    from simgenes.simgenes import check_internal_stops
    seq = 'ATGAAAAAATAAAAATGATGTTTTCCCCGCTAG'
    stops = check_internal_stops(seq)
    assert [3, 5, 10] == stops


def test_prot():
    from simgenes.simgenes import get_prot
    seq = ['', 'ATGTTAA', 'ATGAAATTTAAATAG', 'KKK']
    res = [get_prot(x) for x in seq]
    expect = ['', 'MLX', 'MKFK*', 'X']
    assert res == expect


def test_len():
    from simgenes.simgenes import get_len
    x = get_len(10,1000)
    assert isinstance(x, int)
    assert 10 <= x <= 1000


def checkseq(L, gc, trans):
    from simgenes.simgenes import random_gene
    s, p = random_gene(L=L, gc=gc, trans=trans)
    if trans:
        assert p != None, 'Did not translate sequence'
    else:
        assert p == None, 'Translated sequence'
    x = len(s)
    assert (x%3) == 0, 'Sequence is not codon-structured'
    if gc:
        if gc == -1: gc = 50
        g, c = s.count('G'), s.count('C')
        newgc = 100*(g+c)/x
        dgc = abs(gc - newgc)
        assert dgc <= 10, f'GC content deviates from that asked {gc:.2f}:{newgc:.2f}:{dgc:.2f}'
    if isinstance(L, tuple):
        assert L[0] <= x <= L[1], 'Length is not in the specified range'
    elif isinstance(L, int):
        if L%3 == 0:
            assert x == L, 'Length is not the one specified'


def test_rg():
    L = [(600, 1500), 1002, 1000, None]
    gc = [-1, 75.2, None]
    trans = [True, False]
    for m in L:
        for g in gc:
            for t in trans:
                checkseq(L=m, gc=g, trans=t)


def test_batch():
    import os
    import lzma
    from Bio import SeqIO
    from simgenes.simgenes import batch_random
    batch_random(ofile='output_testseq',
                 L=(30, 300),
                 GC=50.8,
                 transl=True,
                 n=1_00)
    hs = []
    for record in SeqIO.parse(lzma.open('output_testseq.fna.xz', 'rt'),
                              'fasta'):
        hs.append((record.id, str(record.seq)))
    l = [1 for _, x in hs if (x == None) or (len(x) == 0)]
    l = sum(l)
    assert len(hs) == 100, 'Produced less nt sequences than expected'
    assert l == 0, 'Produced empty nt sequences'
    hs = []
    for record in SeqIO.parse(lzma.open('output_testseq.faa.xz', 'rt'),
                              'fasta'):
        hs.append((record.id, str(record.seq)))
    l = [1 for _, x in hs if (x == None) or (len(x) == 0)]
    l = sum(l)
    assert len(hs) == 100, 'Produced less aa sequences than expected'
    assert l == 0, 'Produced empty aa sequences'
    os.remove('output_testseq.faa.xz')
    os.remove('output_testseq.fna.xz')
