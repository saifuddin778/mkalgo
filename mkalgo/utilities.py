import random
import numpy

class utils(object):
    def __init__(self):
        pass


    @staticmethod
    def chunker(ts, l):
        """
        chunks a given ts into n-chunks
        """
        k = len(ts)
        for i in range(0, k, l):
            yield ts[i:i + l]


    @staticmethod
    def random_sample(m):
        rsample = random.sample(range(m), 1)
        for rs in rsample:
            return rs


    @staticmethod
    def std(item):
        s = sum(item)
        mean = s * 1.0/len(item)
        variance = sum([(a-mean)**2 for a in item])/len(item)
        return variance * 0.5



class funcs(object):
    def __init__(self):
        pass


    @staticmethod
    def euclidean(x1, x2):
        d = 0
        for a, b in zip(x1, x2):
            d += (a - b) ** 2
        if isinstance(d, type(numpy.array)) or isinstance(d, type(list)):
            d = sum(d)
        return d ** 0.5


    @staticmethod
    def dtw(x1, x2):
        mapp = {};
        for i, e in enumerate(x1):
            mapp[(i, -1)] = float('inf')
        for i, e in enumerate(x2):
            mapp[(-1, i)] = float('inf')
        mapp[(-1, -1)] = 0

        for i, e in enumerate(x1):
            for j, f in enumerate(x2):
                mapp[(i, j)] = (x1[i] - x2[j])**2  + \
                            min(mapp[(i-1, j)], mapp[(i, j-1)], mapp[(i-1, j-1)])
        return mapp[(len(x1)-1, len(x2)-1)] ** 0.5

