import random
import numpy

from mkalgo.utilities import utils
from mkalgo.utilities import funcs


class base(object):
    def __init__(self, l, metric, r=None):
        self.l = l
        self.r = r
        if metric == 'euclidean':
            self.metric = funcs.euclidean
        if metric == 'dtw':
            self.metric = funcs.dtw

        self.motif_a_index, self.motif_a, self.motif_a_begin, self.motif_a_end = None, None, None, None
        self.motif_b_index, self.motif_b, self.motif_b_begin, self.motif_b_end = None, None, None, None


    def lorder_order(self, lorder, m):
        lorder = sorted(lorder, key=lambda n: n['dist'])
        for j in xrange(m):
            if j == 0:
                lorder[j]['D'] = 0
            else:
                lorder[j]['D'] = lorder[j]['dist'] - lorder[j-1]['dist']
        return lorder


    def get_splits(self, lts):
        splits = []
        m = 0
        for i in range(0, lts, self.l):
            splits.append([i, i+self.l])
            m += 1
        return splits, m



class mk_eab(base):
    def __init__(self, *args, **kwargs):
        super(mk_eab, self).__init__(*args, **kwargs)


    def search(self, ts):
        lts = len(ts)
        lorder = []
        assert type(ts) == type([]) or type(ts) == type(numpy.array([]))
        if not lts >= self.l*4:
            raise 'very large motif size for this time series'
        splits, m = self.get_splits(lts)
        ref_index = utils.random_sample(m)
        ref = ts[splits[ref_index][0] : splits[ref_index][1]]
        self.motif_a_begin = splits[ref_index][0]
        self.motif_a_end = splits[ref_index][1]
        self.motif_a_index = ref_index
        self.motif_a = ref

        bsf = float('inf')

        #--find the distance of subsequences from the reference point
        for j in xrange(m):
            if j == self.motif_a_index:
                lorder.append({
                    'dist': 0, 
                    'index': j, 
                    'begin': splits[j][0], 
                    'end': splits[j][1]
                })
            else:
                dist = self.metric(ref, ts[splits[j][0] : splits[j][1]])
                if dist < bsf:
                    bsf = dist
                    self.motif_b_index = j
                    self.motif_b = ts[splits[j][0] : splits[j][1]]
                    self.motif_b_begin = splits[j][0]
                    self.motif_b_end = splits[j][1]
                lorder.append({
                    'dist': dist, 
                    'index': j, 
                    'begin': splits[j][0], 
                    'end': splits[j][1]
                })

        #--sorting the linear order of distances and introducing the 
        #--delta between consecutive subsequences
        lorder = self.lorder_order(lorder, m)

        #--now find the best motif
        offset = 0
        abandon = False

        while abandon == False:
            offset += 1
            abandon = True
            for j in xrange(m-offset):
                npos = j+offset
                if (lorder[j]['D'] - lorder[npos]['D']) < bsf:
                    abandon = False
                    dist = self.metric(
                        ts[lorder[j]['begin']: lorder[j]['end']], 
                        ts[lorder[npos]['begin'] : lorder[npos]['end']]
                    )
                    if dist < bsf:
                        bsf = dist
                        self.motif_a_begin = lorder[j]['begin']
                        self.motif_a_end = lorder[j]['end']
                        self.motif_a_index = lorder[j]['index']
                        self.motif_a = ts[lorder[j]['begin']: lorder[j]['end']]

                        self.motif_b_begin = lorder[npos]['begin']
                        self.motif_b_end = lorder[npos]['end']
                        self.motif_b_index = lorder[npos]['index']
                        self.motif_b = ts[lorder[npos]['begin'] : lorder[npos]['end']]
        first = {
            'begin': self.motif_a_begin, 
            'end': self.motif_a_end, 
            'motif': self.motif_a
        }

        second = {
            'begin': self.motif_b_begin, 
            'end': self.motif_b_end, 
            'motif': self.motif_b
        }

        return first, second



class mk(base):
    def __init__(self, *args, **kwargs):
        super(mk, self).__init__(*args, **kwargs)


    def search(self, ts):
        lts = len(ts)
        if not lts >= self.l*4:
            raise 'very large motif size for this time series'

        bsf = float('inf')
        self.r = (self.r or 1)
        splits, m = self.get_splits(lts)
        lorder = []
        maxz = -float('inf')
        ref = None

        #--in theory, we can pick the ref with maximum standard deviation
        #--as per the paper, this just includes an additional iteration here
        for i in xrange(self.r):
            ref_index = utils.random_sample(m)
            ref_candidate = ts[splits[ref_index][0] : splits[ref_index][1]]
            #--standard deviations
            std = utils.std(ref_candidate)
            if std > maxz:
                maxz = std
                ref = ref_candidate
                self.motif_a_begin = splits[ref_index][0]
                self.motif_a_end = splits[ref_index][1]
                self.motif_a = ref
                self.motif_a_index = ref_index

        #--now you can continue with your regular algorithm
        bsf = float('inf')

        #--find the distance of subsequences from the reference point
        for j in xrange(m):
            if j == self.motif_a_index:
                lorder.append({
                    'dist': 0, 
                    'index': j, 
                    'begin': splits[j][0], 
                    'end': splits[j][1]
                })
            else:
                dist = self.metric(ref, ts[splits[j][0] : splits[j][1]])
                if dist < bsf:
                    bsf = dist
                    self.motif_b_begin = splits[j][0]
                    self.motif_b_end = splits[j][1]
                    self.motif_b_index = j
                    self.motif_b = ts[splits[j][0] : splits[j][1]]
                lorder.append({
                    'dist': dist, 
                    'index': j, 
                    'begin': splits[j][0], 
                    'end': splits[j][1]
                })

        #--sorting the linear order of distances and introducing the 
        #--delta between consecutive subsequences
        lorder = self.lorder_order(lorder, m)

        #--now find the best motif
        offset = 0
        abandon = False

        while abandon == False:
            offset += 1
            abandon = True
            for j in xrange(m-offset):
                npos = j+offset
                if (lorder[j]['D'] - lorder[npos]['D']) < bsf:
                    abandon = False
                    dist = self.metric(
                        ts[lorder[j]['begin']: lorder[j]['end']], 
                        ts[lorder[npos]['begin'] : lorder[npos]['end']]
                    )
                    if dist < bsf:
                        bsf = dist
                        self.motif_a_begin = lorder[j]['begin']
                        self.motif_a_end = lorder[j]['end']
                        self.motif_a_index = lorder[j]['index']
                        self.motif_a = ts[lorder[j]['begin']: lorder[j]['end']]

                        self.motif_b_begin = lorder[npos]['begin']
                        self.motif_b_end = lorder[npos]['end']
                        self.motif_b_index = lorder[npos]['index']
                        self.motif_b = ts[lorder[npos]['begin'] : lorder[npos]['end']]
        first = {
            'begin': self.motif_a_begin, 
            'end': self.motif_a_end, 
            'motif': self.motif_a
        }

        second = {
            'begin': self.motif_b_begin, 
            'end': self.motif_b_end, 
            'motif': self.motif_b
        }

        return first, second
