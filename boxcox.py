from scipy import stats
import matplotlib.pyplot as plt
import sys

fh = open(sys.argv[1], 'r')
fh.readline()
flist = []
slist = []
tlist = []
for line in fh:
    new = line.strip().split('\t')
    flist.append(float(new[1]))
    slist.append(float(new[2]))
    tlist.append(float(new[3]))

stlist = []
for i in tlist:
    n = i + abs(min(tlist)) + 0.1
    stlist.append(n)

fig = plt.figure()
ax1 = fig.add_subplot(321)
prob = stats.probplot(flist, dist=stats.norm, plot=ax1)
ax1.set_title('Probplot against normal distribution')
ax2 = fig.add_subplot(322)
xt, l = stats.boxcox(flist)
print l
stats.probplot(xt, dist=stats.norm, plot=ax2)

ax1 = fig.add_subplot(323)
prob = stats.probplot(slist, dist=stats.norm, plot=ax1)
ax1.set_title('Probplot against normal distribution')
ax2 = fig.add_subplot(324)
xt, l = stats.boxcox(slist)
print l
stats.probplot(xt, dist=stats.norm, plot=ax2)

ax1 = fig.add_subplot(325)
prob = stats.probplot(tlist, dist=stats.norm, plot=ax1)
ax1.set_title('Probplot against normal distribution')
ax2 = fig.add_subplot(326)
xt, l = stats.boxcox(stlist)
print l
stats.probplot(xt, dist=stats.norm, plot=ax2)
plt.show()
