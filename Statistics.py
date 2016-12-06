"""
Put Relevant Plots and Graphs Here:
Basic:
Min/Max/Mean per Generation
Histogram of Score per Generation
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from GeneticAlgorithm import runit
from Data import Data

def sampler(test):
    _mean = []
    _max = []
    _min = []
    _best = []
    _times = []
    _scores = []
    run_times, score_watch, best_perc = runit(test)
    for run_time, scores in zip(run_times,score_watch):
        _times.append(run_time)
        _mean.append(scores[2])
        _max.append(scores[3])
        _min.append(scores[1])
        _best.append(scores[0])
        _scores.append(scores[4])
    return _times, _mean, _max, _min, _best,_scores, best_perc

def scatter(test):
    _times, _mean, _max, _min, _best = sampler(test)
    plt.figure()
    plt.plot(range(len(_times)),_mean)
    plt.plot(range(len(_times)),_max)
    plt.plot(range(len(_times)),_min)
    plt.xlim([0,len(_times)])
    plt.ylim([min(_min),max(_max)])
    plt.legend(["Mean","Max","Min"])
    plt.title("Mean, Max, and Min Scores ({0})".format(test.upper()))
    plt.xlabel("Generations")
    plt.ylabel("Scores")
    plt.grid(True)
    plt.savefig("scatter_{0}.png".format(test))
    plt.close()

def histogram(test):
    _times, _mean, _max, _min, _best, _scores, best_perc = sampler(test)
    f,a = plt.subplots(2,2)
    a = a.ravel()
    data = []
    titles = []
    for i in range(0,len(_times),len(_times)/4):
        data.append(_scores[i])
        titles.append("Generation {0}".format(i))
    for idx,ax in enumerate(a):
        ax.hist(data[idx])
        ax.set_xlim([min(_min),max(_max)])
        ax.set_title(titles[idx])
        ax.set_xlabel("Density")
        ax.set_ylabel("Scores")
    plt.tight_layout()
    plt.savefig("histogram_{0}.png".format(test))

def vary_mutation():
    return

def vary_crossover():
    return

def percentage_bar(tests,runs=5,width=.35):
    string_size = []
    avg_perc = []
    plt.figure()
    for test in tests:
        percs = []
        string_size.append(len(Data("{0}.txt".format(test.lower())).getSequence()))
        for run in xrange(runs):
            _times, _mean, _max, _min, _best, _scores, perc = sampler(test)
            percs.append(object)
        avg_perc.append(mean(percs))
        errors.append(numpy.std(percs))
    plt.bar(range(len(string_size)),avg_perc, width, color='r', yerr=menStd)
    plt.title("Accuracy (Percentage) vs. RNA Length")
    plt.xticks(ind + width)
    plt.xticklabels(string_size)
    plt.ylabel("Percentage")
    plt.xlabel("RNA Length")

def time_bar(tests,runs=5,width=.35):
    string_size = []
    avg_perc = []
    plt.figure()
    for test in tests:
        percs = []
        string_size.append(len(Data("{0}.txt".format(test.lower())).getSequence()))
        for run in xrange(runs):
            _times, _mean, _max, _min, _best, _scores, perc = sampler(test)
            percs.append(object)
        avg_perc.append(mean(percs))
        errors.append(numpy.std(percs))
    plt.bar(range(len(string_size)),avg_perc, width, color='r', yerr=menStd)
    plt.title("Accuracy (Percentage) vs. RNA Length")
    plt.xticks(ind + width)
    plt.xticklabels(string_size)
    plt.ylabel("Percentage")
    plt.xlabel("RNA Length")


def all_stats():
    tests = ["test6","test7","test4","test5"]
    for test in tests:
        #scatter(test)
        #histogram(test)
    percentage_bar(test)
    time_bar(test)


if __name__ == '__main__':
    all_stats()