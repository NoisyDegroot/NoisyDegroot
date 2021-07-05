import sys
from scipy import stats
import numpy
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
numpy.random.seed(sum(map(ord, "aesthetics")))

legend_properties = {'size' : 13, 'weight':'bold'}
plt.rcParams['font.family'] = "Times New Roman"


def get_dev(r, mi = 0):
    a = r.values()
    return stats.variation(a)


def main():

    rst = []

    rst.append(('Windsurfers', [0.9375, 0.0625, 0.0, 0.0, 0.0]))
    rst.append(('CA-HepPh', [0.976449, 0.0171794, 0.00560919, 0.000762225, 0.0]))
    rst.append(('Facebook(NIPS)', [0.990175, 0.00442576, 0.00106218, 0.000663864, 0.00367338]))
    rst.append(('Hamster full', [0.870484, 0.0829852, 0.0198247, 0.0105677, 0.0161383]))
    rst.append(('Jazz musicians', [0.971239, 0.0195428, 0.00626844, 0.00221239, 0.000737463]))
    rst.append(('Reactome', [0.990658, 0.00665735, 0.00145202, 0.000301362, 0.000931481]))
    rst.append(('Hamsterster Friendships', [0.832379, 0.106338, 0.0259641, 0.0124093, 0.0229095]))
    rst.append(('Pretty Good Privacy', [0.880963, 0.0435071, 0.0192138, 0.00817138, 0.0481449]))

    width = 1
    size = 5
    x = [y * width * 10 for y in numpy.arange(size)]
    cur_pal = sns.color_palette("deep")

    xticks = [
        "(0,0.1]",
        "(0.1,0.2]",
        "(0.2,0.3]",
        "(0.3,0.4]",
        " > 0.4",
    ]

    for i in range(len(rst)):
        #print(map(lambda x: x + i * width, x))
        plt.bar([x + i * width for x in x], [x * 100 for x in rst[i][1]], width=width, label=rst[i][0], color=cur_pal[i])

    plt.axis()
    plt.xlim((- width * 2.5, size * width * 5.5))
    plt.ylim((0, 100))

    ax = plt.subplot(1,1,1)
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['top'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)

    '''
    plt.xticks([x + 2 * width for x in x], xticks)#, fontsize = 16, fontweight='bold')
    #plt.yticks(fontsize=16, fontweight='bold')
    '''
    plt.xticks([x + 2 * width for x in x], xticks, fontsize=14, fontweight='bold', rotation=25)
    plt.yticks(fontsize=16, fontweight='bold')

    plt.ylabel('Percentage of Edges', fontsize = 18, fontweight='bold')
    plt.xlabel('Relative Error', fontsize = 18, fontweight='bold')

    plt.legend(loc="upper right")#, prop=legend_properties)

    #plt.draw()
    #plt.savefig('deviation.eps', format='eps')

    plt.show()


if __name__ == "__main__":
    main()
