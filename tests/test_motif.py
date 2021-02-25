from mkalgo.mk import mk_eab
from mkalgo.data import hospital

def test_motif():
  obj = mk_eab(l=5, metric='euclidean')
  x = hospital()
  motif_a, motif_b = obj.search(x)
  
