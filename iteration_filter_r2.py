import sys
import scipy
import itertools
import numpy as np
from scipy import stats

#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Disjoint set data structure <http://code.activestate.com/recipes/387776/>
Author: Michael Droettboom
"""

class Grouper(object):
    """
    This class provides a lightweight way to group arbitrary objects
    together into disjoint sets when a full-blown graph data structure
    would be overkill. 
    
    Objects can be joined using .join(), tested for connectedness
    using .joined(), and all disjoint sets can be retrieved using list(g)
    The objects being joined must be hashable.
    >>> g = Grouper()
    >>> g.join('a', 'b')
    >>> g.join('b', 'c')
    >>> g.join('d', 'e')
    >>> list(g)
    [['a', 'b', 'c'], ['d', 'e']]
    >>> g.joined('a', 'b')
    True
    >>> g.joined('a', 'c')
    True
    >>> 'f' in g
    False
    >>> g.joined('a', 'd')
    False
    """   
    def __init__(self, init=[]):
        mapping = self._mapping = {}
        for x in init:
            mapping[x] = [x]
        
    def join(self, a, *args):
        """
        Join given arguments into the same set. Accepts one or more arguments.
        """
        mapping = self._mapping
        set_a = mapping.setdefault(a, [a])

        for arg in args:
            set_b = mapping.get(arg)
            if set_b is None:
                set_a.append(arg)
                mapping[arg] = set_a
            elif set_b is not set_a:
                if len(set_b) > len(set_a):
                    set_a, set_b = set_b, set_a
                set_a.extend(set_b)
                for elem in set_b:
                    mapping[elem] = set_a


    def joined(self, a, b):
        """
        Returns True if a and b are members of the same set.
        """
        mapping = self._mapping
        try:
            return mapping[a] is mapping[b]
        except KeyError:
            return False


    def __iter__(self):
        """
        Returns an iterator returning each of the disjoint sets as a list.
        """
        seen = set()
        for elem, group in self._mapping.iteritems():
            if elem not in seen:
                yield group
                seen.update(group)


    def __getitem__(self, key):
        """
        Returns the set that a certain key belongs.
        """
        return tuple(self._mapping[key])


    def __contains__(self, key):
        return key in self._mapping


    def __len__(self):
        group = set()
        for v in self._mapping.values():
           group.update([tuple(v)])
        return len(group)


if __name__ == '__main__':

    import doctest
    doctest.testmod()
#python iteration_filter.py NIC_allgenotypes_hybrid_beagle_filtermaf_zlname_classicalgene.vcf NIC_allgenotypes_hybrid_beagle_filte
#rmaf_zlname_classicalgene_rmSimilarity.vcf
fh = open(sys.argv[1],'r')
head = fh.readline()
mdict = {}
for line in fh:
	new = line.strip().split('\t')
	if new[0] not in mdict:
		mdict[new[0]] = {}
	if new[1] not in mdict[new[0]]:
		mdict[new[0]][new[1]] = map(float,new[2:])
fh.close()

out = open(sys.argv[2],'w')
out.write(head)
for i in mdict:
	if len(mdict[i]) == 1:
		for m in mdict[i]:
#			print i,m
			final = map(str,mdict[i][m])
			out.write(i+'\t'+m+'\t'+'\t'.join(final)+'\n')
	else:
#		print i
		tf = set([])
		#print list(mdict[i])
		g = Grouper()
		for j in itertools.combinations(list(mdict[i]),2):
			perc = 0
			tt = len(mdict[i][j[0]])
			nlist = []
			k = scipy.stats.pearsonr(mdict[i][j[0]],mdict[i][j[1]])
			r = k[0]*k[0]
			if r > 0.8:
				g.join(j[0],j[1])
			else:
				g.join(j[0])
				g.join(j[1])
		#print sorted(list(g))
		for m in sorted(list(g)):
			snp = m[0]
			final = map(str,mdict[i][snp])
			out.write(i+'\t'+snp+'\t'+'\t'.join(final)+'\n')
out.close()	
