from mkalgo.mk import mk_eab
from mkalgo.data import hospital
from pprint import pprint

if __name__=='__main__':
  obj = mk_eab(l=5, metric='euclidean')
  x = hospital()
  motif_a, motif_b = obj.search(x)
  pprint(motif_a)
  pprint(motif_b)
