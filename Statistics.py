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
import numpy
import pylab
import matplotlib.cm as cm
from PIL import Image

def sampler(test,mutation=0.1,crossover=0.7):
    _mean = []
    _max = []
    _min = []
    _best = []
    _times = []
    _scores = []
    run_times, score_watch, best_perc = runit(test,mutation=mutation,crossover=crossover)
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
        ax.set_ylabel("Density")
        ax.set_xlabel("Scores")
    plt.tight_layout()
    plt.savefig("histogram_{0}.png".format(test))

def vary_mutation(test,mutations,runs=3):
    means = []
    errors = []
    fig, ax = plt.subplots()
    for mutation in mutations:
        avg_perc = []
        for run in range(runs):
            _times, _mean, _max, _min, _best, _scores, best_perc = sampler(test,mutation=mutation)
            avg_perc.append(best_perc)
        means.append(numpy.mean(avg_perc))
        errors.append(numpy.std(avg_perc))
    plt.bar(range(len(mutations)),means, width, color='r', yerr=errors)
    ax.set_title("Accuracy (Percentage) vs. Mutation Rate ({0})".format(test.upper()))
    ax.set_xticks(range(len(mutations)))
    ax.xticklabels(mutations)
    ax.set_ylabel("Percentage")
    ax.set_xlabel("Mutation Rate")
    plt.savefig("mutation_{0}.png".format(test.lower()))

def vary_crossover(test,crossovers,runs=3):
    means = []
    errors = []
    fig, ax = plt.subplots()
    for crossover in crossovers:
        avg_perc = []
        for run in range(runs):
            _times, _mean, _max, _min, _best, _scores, best_perc = sampler(test,crossover=crossover)
            avg_perc.append(best_perc)
        means.append(numpy.mean(avg_perc))
        errors.append(numpy.std(avg_perc))
    plt.bar(range(len(crossovers)),means, width, color='r', yerr=errors)
    ax.set_title("Accuracy (Percentage) vs. Crossover Rate ({0})".format(test.upper()))
    ax.set_xticks(range(len(crossovers)))
    ax.set_xticklabels(crossovers)
    ax.set_ylabel("Percentage")
    ax.set_xlabel("Crossover Rate")
    plt.savefig("crossover_{0}.png".format(test.lower()))

def percentage_bar(tests,runs=3,width=.35):
    string_size = []
    avg_perc = []
    errors = []
    fig, ax = plt.subplots()
    for test in tests:
        percs = []
        string_size.append(len(Data("{0}.txt".format(test.lower())).getSequence()))
        for run in xrange(runs):
            _times, _mean, _max, _min, _best, _scores, perc = sampler(test)
            percs.append(perc)
        avg_perc.append(numpy.mean(percs))
        errors.append(numpy.std(percs))
    ax.bar(range(len(string_size)),avg_perc, width, color='r', yerr=errors)
    ax.set_title("Accuracy (Percentage) vs. RNA Length")
    ax.set_xticks(range(len(string_size)))
    ax.set_xticklabels(string_size)
    ax.set_ylabel("Percentage")
    ax.set_xlabel("RNA Length")
    plt.savefig("percentage.png")

def time_bar(tests,runs=3,width=.35):
    string_size = []
    avg_times = []
    errors = []
    fig, ax = plt.subplots()
    for test in tests:
        times = []
        string_size.append(len(Data("{0}.txt".format(test.lower())).getSequence()))
        for run in xrange(runs):
            _times, _mean, _max, _min, _best, _scores, perc = sampler(test)
            times += _times
        avg_times.append(numpy.mean(times))
        errors.append(numpy.std(times))
    ax.bar(range(len(string_size)),avg_times, width, color='r', yerr=errors)
    ax.set_title("Average Time per Generation vs. RNA Length")
    ax.set_xticks(range(len(string_size)))
    ax.set_xticklabels(string_size)
    ax.set_ylabel("Seconds")
    ax.set_xlabel("RNA Length")
    plt.savefig("time.png")


def iter_bar(tests,runs=3,width=.35):
    string_size = []
    avg_times = []
    errors = []
    fig, ax = plt.subplots()
    for test in tests:
        iters = []
        string_size.append(len(Data("{0}.txt".format(test.lower())).getSequence()))
        for run in xrange(runs):
            _times, _mean, _max, _min, _best, _scores, perc = sampler(test)
            iters.append(len(_times))
        avg_times.append(numpy.mean(iters))
        errors.append(numpy.std(iters))
    ax.bar(range(len(string_size)),avg_times, width, color='r', yerr=errors)
    ax.set_title("Average Iterations (Convergence) vs. RNA Length")
    ax.set_xticks(range(len(string_size)))
    ax.set_xticklabels(string_size)
    ax.set_ylabel("Average Iterations")
    ax.set_xlabel("RNA Length")
    plt.savefig("time.png")

def compare_diagrams(test_image,true_image):
    f = plt.figure(figsize=(4,8))
    perc = 40.506
    titles=["Predicted (Accuracy={0}%)".format(perc),"True Structure"]
    for n, fname in enumerate((test_image, true_image)):
        image=Image.open(fname).convert("L")
        arr=numpy.asarray(image)
        f.add_subplot(2, 1, n+1)
        plt.axis('off')
        plt.tick_params(axis='both', left='off', top='off', right='off', bottom='off', labelleft='off', labeltop='off', labelright='off', labelbottom='off')
        plt.imshow(arr,cmap=cm.Greys_r)
        plt.title(titles[n])
    plt.tight_layout()
    plt.savefig("comparison.png")

def all_stats():
    tests = ["test6","test7","test4","test5"]
    # for test in tests:
    #     #scatter(test)
    #     histogram(test)
    percentage_bar(tests)
    #time_bar(tests)
    #iter_bar(tests)
    #vary_mutation("test4", [0.01,0.05,0.1,0.15,0.2,0.3])
    #vary_crossover("test4", [0.2,0.3,0.4,0.5,0.6,0.7,0.8])
    # compare_diagrams("Test5Output_ss.ps","test5_true.png")

if __name__ == '__main__':
    all_stats()