import petls
import numpy as np
import scipy
import tadasets
import timeit 
import pandas as pd
import statistics

def harmonic_stats(nonzero_eigs):
    if len(nonzero_eigs) == 0:
        return [0, 0, 0, 0, 0]
    min_eig = nonzero_eigs[0]
    max_eig = nonzero_eigs[-1]
    avg_eig = sum(nonzero_eigs) / len(nonzero_eigs)
    med_eig = statistics.median(nonzero_eigs)
    harm_eig = statistics.harmonic_mean(nonzero_eigs)
    return [min_eig, max_eig, avg_eig, med_eig, harm_eig]

def spectra_request():
    # every 0.025, from 0 to 5
    step = 0.25
    max_filtration = 5
    samples = [step*i for i in range(int(max_filtration/step)+1)]
    dims = [0, 1, 2, 3]
    requests = []
    for i in range(len(samples)-1):
        for dim in dims:
            request = (dim, samples[i], samples[i+1])    
            requests.append(request)
    return requests

def t_spec(n):
    torus = tadasets.torus(n=n, c=3, a=2, noise = 0.0)
    complex_torus = petls.Alpha(torus,3)


    all_spectra = complex_torus.spectra(spectra_request())
    all_stats = []
    for spectrum in all_spectra:
        nonzero_eigs = [x for x in spectrum[3] if x > 0]
        # print(nonzero_eigs)
        spectrum_stats = harmonic_stats(nonzero_eigs)
        betti = len([x for x in spectrum if x == 0])
        spectrum_stats.append(betti)
        all_stats.append([spectrum[0], spectrum[1], spectrum[2]] + spectrum_stats)
        # print(all_stats[-1])
    df = pd.DataFrame(all_stats,columns=["dim","a","b","min","max","avg","med","harmonic_mean","betti"])
    return df

import matplotlib.pyplot as plt
def plot_full_df(df, n):
    dims = df["dim"].unique()
    for dim in dims:
        df_dim = df.loc[df["dim"] == dim]
        df_dim.to_csv(f"subsamples_df_{dim}_{n}.csv", index = False)
        # plt.scatter(x=df_dim["a"], y=df_dim["min"])
        # plt.show()



execution_times = {}
# point_samples = [1200,1000,800,600,500,400,300]
# point_samples = [1200]
# point_samples = [300]
# point_samples = [20, 30]
# point_samples = 
# point_samples = [i for i in range(5,51,1)] + [i for i in range(55,101,5)]
# point_samples = [125, 150, 175, 200, 250] # already done 300 
# point_samples = [350, 400, 500]
point_samples = [200]
def parallel(point_samples):
    from multiprocessing import Pool
    with Pool() as p:
        dfs = p.map(t_spec,point_samples)
        for i in range(len(point_samples)):
            plot_full_df(dfs[i],point_samples[i])

def timed(point_samples):
    import time

    for n in point_samples:
        start_time = time.time()
        df = t_spec(n)
        plot_full_df(df, n)
        end_time = time.time()
        execution_times[n] = end_time - start_time

    for key, value in execution_times.items():
        print(f"n={key}, time={value}")

    series = pd.Series(execution_times)
    print(series)
    series.to_excel("subsample_execution_times_2.xlsx")

# parallel(point_samples)