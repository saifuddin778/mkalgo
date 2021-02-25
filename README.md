## mkalgo (MK Algorithm)
This package provides Pythonic (>=2.6) implementation of Mueen-Keogh (MK) algorithm for finding motifs in the time series.

### Citation
_Mueen, A., Keogh, E., Zhu, Q., Cash, S. and Westover, B. Exact Discovery of Time Series Motif. SDM 2009_ ([link](http://alumni.cs.ucr.edu/~mueen/pdf/EM.pdf)).

![]ttps://i.imgur.com/MfbxULV.png)

### Explanation
This algorithm aims to find significant pairs of motifs (identical subsequences in a time-series which repeat more than once) in a given time-series. The usual bruteforce methods of finding these motifis are usually O(n^2) in terms of their complexity, while MK Algorithm aims to find them in O(n * log(n)) order.

### Usage
This implementation contains both MK-Early Abandoning and standard MK algorithm implementations:
```python
>>> from mkalgo.mk import mk_eab, mk
```

##### MK-EAB

```python
>>> obj = mk_eab(l=5, metric='euclidean')
>>> motif_a, motif_b = obj.search(timeseries)
```

##### MK
```python
>>> obj = mk_eab(l=5, metric='euclidean', r=10)
>>> motif_a, motif_b = obj.search(timeseries)
```

* `l` refers to the desired length of motif pairs to be found
* `metric` refers to the distance metric to be utilized, which can be `euclidean` or `dtw` (dynamic time warping).


The returned object is a dictionary, containing the respective motif and its beginning and ending indices in the given timeseries.

Both the implementations can be used in the similar fashion, with `mk` taking an additional parameter of `r`, which refers to the number of rounds made over `n/l` chunks of the time-series to find the significant reference point.
