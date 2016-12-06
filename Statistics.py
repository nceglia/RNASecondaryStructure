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

def sampler(test,runs=10):
    _mean = []
    _max = []
    _min = []
    _best = []
    _times = []
    run_times, score_watch, best_perc = runit(test)
    for run_time, scores in zip(run_times,score_watch):
        _times.append(run_time)
        _mean.append(scores[2])
        _max.append(scores[3])
        _min.append(scores[1])
    return _times, _mean, _max, _min, _best

def scatter(test,runs,prefix):
    _times, _mean, _max, _min, _best = sampler(test,runs)
    plt.plot(range(len(_times)),_mean)
    plt.plot(range(len(_times)),_max)
    plt.plot(range(len(_times)),_min)
    plt.plot(range(len(_times)),_best)
    plt.legend(["Mean","Max","Min","Best"])
    plt.title("Mean, Max, Min, and Best Scores")
    plt.xlabel("Generations")
    plt.ylabel("Scores")
    plt.savefig("scatter.png")

def score_histogram():
    return

def vary_mutation():
    return

def vary_crossover():
    return

def percentage_per_size():
    return

def time_with_size():
    return

def all_stats():
    scatter("test6", 10, "test6")

if __name__ == '__main__':
    all_stats()