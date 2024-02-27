import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
import networkx as nx
from collections import defaultdict
from scipy import stats
import multiprocessing
import datetime
import graph_tool.all as gt
import scipy as sp
from itertools import product
import os

def STREAMLINE_multicore(inputSettings, mode, ncores = 12, verbose=True, directed=False, saveraw=True):
    global opts
    opts = dotdict({
        'datadir': str(inputSettings.datadir),
        'outDir': "outputs/"+str(inputSettings.datadir).split("inputs/")[1]+ '/',
        'verbose': verbose,
        'mode': mode,
        'directed': directed,
        'saveraw': saveraw,
        'str_directed': ("directed" if directed else "undirected"),
        'datasets': [dataset["name"] for dataset in inputSettings.datasets],
        'algos': [algo[0] for algo in inputSettings.algorithms if algo[1]['should_run'] == True]
    })
    dataDict = defaultdict(lambda: defaultdict(lambda: defaultdict(np_nan)))

    multiprocessing.set_start_method("fork", force=True)
    ncores = min(ncores, len(opts.algos) * len(opts.datasets))
    if opts.verbose: print("\n## Running on " + str(ncores) + " cores", datetime.datetime.now())

    if opts.mode == "local":
        if opts.directed:
            opts.metrics = ['PageRank', 'Out Centrality', 'Betweenness', 'Clustering Coefficient', 'Node In Degrees']
        else:
            opts.metrics = ['PageRank', 'Centrality', 'Betweenness', 'Radiality', 'Clustering Coefficient', 'Local Efficiency', 'Node Degrees']
    elif opts.mode == "global":
        if opts.directed:
            opts.metrics = ['Clustering Coefficient', 'Assortativity', 'In Centralization']
        else:
            opts.metrics = ['Clustering Coefficient', 'Global Efficiency', 'Local Efficiency', 'Shortest Path Length', 'Assortativity', 'Centralization']

    runners = list(product(opts.algos + ['ground truth'], inputSettings.datasets))

    ## computation of metrics
    pool = multiprocessing.Pool(ncores)
    results = pool.map(func=computeMetrics, iterable=runners, chunksize=1)
    pool.close()
    pool.join()
    for idx, runner in enumerate(runners):
        dataDict[runner[1]["name"]][runner[0]] = results[idx]
        
    ## computation of summary statistics
    if opts.mode == "global":
        evals = ["mse", "mse_rel"]
        if opts.saveraw:
            save_csv(dataDict, opts, "raw")
    elif opts.mode == "local":
        evals = ["jaccard_ratio", "jaccard", "correlation"] 

    for eval in evals:
        resDict = compute_summary(dataDict, mode=eval)
        save_csv(resDict, opts, eval)

    return

def computeMetrics(runner):
    metricDict = defaultdict(np_nan)
    algo, dataset = runner

    trueEdgesDF = load_data(path=(opts.datadir +'/'+ dataset["name"] + '/' + dataset["trueEdges"]), sep=",", opts=opts)
    refGraph = nx_largest_subgraph(trueEdgesDF, opts)

    predPath = opts.outDir + dataset["name"] + '/' + algo+'/rankedEdges.csv'

    if algo == 'ground truth':
        predGraph = refGraph
    elif Path(predPath).exists():

        predDF = load_data(path=predPath, sep = '\t', opts=opts)
        predDF.EdgeWeight = predDF.EdgeWeight.abs()
        #predDF = predDF[predDF['EdgeWeight'] > 0] #remove edges with EdgeWeight 0
        predDF = predDF[predDF['Gene1'].isin(refGraph.nodes) & predDF['Gene2'].isin(refGraph.nodes)]
        predDF.reset_index(drop=True, inplace=True)

        if predDF.shape[0] > 0:
        
            numEdges = refGraph.number_of_edges()
            if predDF.shape[0] < numEdges:
                numEdges = predDF.shape[0]
                print("Warning - Less edges available from", algo, "than in ground truth", dataset["name"])
            newDF = predDF.loc[0:numEdges-1]
            predGraph = nx_largest_subgraph(newDF, opts)

            if opts.verbose: print("Finished loading predicted network " + dataset["name"], datetime.datetime.now())
            
        else:
            print("Warning - No edges predicted for:", algo, "on", dataset["name"])
            return {}
    else:
        print(predPath, ' does not exist. Skipping...')
        return {}
        
    ## calculate actual metrics
    for metric in opts.metrics:
        if opts.mode == "local":
            metricDict[metric] = get_local_metric(metric, predGraph)
        elif opts.mode == "global":
            metricDict[metric] = get_global_metric(metric, predGraph)
        
    rawValues = pd.DataFrame(metricDict) if opts.mode == "local" else pd.DataFrame([metricDict])
        
    ## save raw data if needed
    if opts.saveraw:
        savepath = opts.outDir + dataset["name"] + '/' + algo
        if os.path.isdir(savepath)==False:
            os.mkdir(savepath)
        rawValues.to_csv(savepath + f"/STREAMLINE_{opts.mode}_{opts.str_directed}_raw.csv", index=(opts.mode == "local"))

    if opts.verbose: print("Finished", algo, "on", dataset["name"], datetime.datetime.now())

    return metricDict

def load_data(path, sep, opts):
    edgesDF = pd.read_csv(path, sep = sep, header = 0, index_col = None)
    
    edgesDF = edgesDF[edgesDF['Gene1'] != edgesDF['Gene2']] #remove self-loops
    if not opts.directed:
        edgesDF.loc[edgesDF['Gene1'] > edgesDF['Gene2'], edgesDF.columns[[0,1]]] = edgesDF.loc[edgesDF['Gene1'] > edgesDF['Gene2'], edgesDF.columns[[1,0]]].values
    edgesDF.drop_duplicates(subset=['Gene1', 'Gene2'], keep = 'first', inplace=True)
    edgesDF.reset_index(drop=True, inplace=True)
    
    return edgesDF

def nx_largest_subgraph(edgeDF, opts):
    if opts.directed:
        Graph = nx.from_pandas_edgelist(edgeDF, source='Gene1', target='Gene2', create_using=nx.DiGraph())
        Graph = max((Graph.subgraph(c) for c in nx.weakly_connected_components(Graph)), key=len)
    else:
        Graph = nx.from_pandas_edgelist(edgeDF, source='Gene1', target='Gene2', create_using=nx.Graph())
        Graph = max((Graph.subgraph(c) for c in nx.connected_components(Graph)), key=len)
    return Graph

def get_local_metric(metric, Graph):

    if metric=='PageRank':
        scores = pagerank(Graph)
    elif metric=='Centrality':
        scores = centrality(Graph)
    elif metric=='Out Centrality':
        scores = out_centrality(Graph)
    elif metric=='Betweenness':
        scores = betweenness(Graph) 
    elif metric=='Radiality':
        scores = radiality(Graph)
    elif metric=='Clustering Coefficient':
        scores = node_clustering_coefficient(Graph)
    elif metric=='Local Efficiency':
        scores = node_local_efficiency(Graph)
    elif metric=='Node Degrees':
        scores = node_degrees(Graph)
    elif metric=='Node In Degrees': #directed only
        scores = node_in_degrees(Graph)
    
    return scores

def get_global_metric(metric, Graph):

    if metric=='Clustering Coefficient':
        score = clustering_coefficient(Graph)
    elif metric=='Global Efficiency':
        score = global_efficiency(Graph)
    elif metric=='Local Efficiency':
        score = local_efficiency(Graph)
    elif metric=='Shortest Path Length':
        score = average_shortest_path(Graph)
    elif metric=='Assortativity':
        score = assortativity(Graph)
    elif metric=='Centralization':
        score = centralization(Graph)
    elif metric=='In Centralization': #directed only
        score = in_centralization(Graph)

    return score

### Hub identification metrics

def radiality(G):
    if nx.is_directed(G) or not(nx.is_connected(G)):
        return np.nan
    
    G_gt = nx_to_gt(G)

    all_sp = gt.shortest_distance(G_gt)
    sum_sp =  [sum(i) for i in all_sp]
    c = 1 + gt.pseudo_diameter(G_gt)[0]
    denom = G_gt.num_vertices()-1
    rad = [(c - i/denom) for i in sum_sp]
    rad_dict = {list(G.nodes)[i]: rad[i] for i in range(G_gt.num_vertices())}
    return rad_dict

def centrality(G):
    scores = nx.degree_centrality(G)
    return scores

def out_centrality(G):
    scores = nx.out_degree_centrality(G)
    return scores

def betweenness(G):
    #better for directed, but defined for undirected
    G_gt = nx_to_gt(G)

    vp, _ = gt.betweenness(G_gt, norm=True)
    vp_dict = {list(G.nodes)[i]: vp[i] for i in range(G_gt.num_vertices())}
    return(vp_dict)

def pagerank(G):
    #undirected wil be converted to directed automatically
    scores = nx.pagerank(G, max_iter=500) 
    return scores

def node_clustering_coefficient(G):
    scores = nx.clustering(G)
    return scores

def node_local_efficiency(G):
    scores = dict((node, global_efficiency(G.subgraph(G[node]))) for node in G)
    return scores

def node_degrees(G):
    scores = dict(G.degree)
    return scores

def node_in_degrees(G):
    scores = dict(G.in_degree)
    return scores

## Hub topology metrics

def centralization(G):
    N = G.order()
    degrees = [x[1] for x in list(G.degree)]
    score = float((N*max(degrees) - sum(degrees)))/((N-1)*(N-2))
    return score

def in_centralization(G):
    N = G.order()
    degrees = [x[1] for x in list(G.in_degree)]
    score = float((N*max(degrees) - sum(degrees)))/((N-1)*(N-1))
    return score

def clustering_coefficient(G):
    score = nx.average_clustering(G)
    return score

def assortativity(G):
    if nx.is_directed(G):
        score = nx.degree_assortativity_coefficient(G, x='in', y='in')
    else:
        score = nx.degree_assortativity_coefficient(G)
    return score

## Information exchange metrics

def global_efficiency(G):
    n = len(G)
    denom = n*(n-1)

    if denom>0:
        G_gt = nx_to_gt(G)
        counts, bins = gt.distance_histogram(G_gt)
        efficiency = sum(counts[1:]/bins[1:-1]) / denom
        return efficiency
    else:
        return 0

def local_efficiency(G):
    efficiency_list = (global_efficiency(G.subgraph(G[v])) for v in G)
    return sum(efficiency_list) / len(G)

def average_shortest_path(G):
    if not(nx.is_connected(G)):
        return np.nan
    
    G_gt = nx_to_gt(G)

    all_sp = gt.shortest_distance(G_gt)
    avg_path = sum([sum(i) for i in all_sp])/(G_gt.num_vertices()**2-G_gt.num_vertices())

    return avg_path

## Utility functions

class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
    
def np_nan(): 
    return np.nan

def nx_to_gt(G):
    G_gt = gt.Graph(directed=G.is_directed())

    vertices = {}
    for node in G.nodes:
        v = G_gt.add_vertex()
        vertices[node] = v

    for src, dst in G.edges():
        G_gt.add_edge(vertices[src], vertices[dst])

    return(G_gt)

def jaccard(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    return float(intersection) / union

def save_csv(resDict, opts, eval):
    datasets = list(resDict.keys())
    for dataset in datasets:
        resDF = pd.DataFrame(resDict[dataset]).T
        resDF.to_csv(opts.outDir + dataset + f"/STREAMLINE_{opts.mode}_{opts.str_directed}_{eval}.csv")


## Evaluation summaries
def compute_summary(dataDict, mode):
    resDict = defaultdict(lambda: defaultdict(lambda: defaultdict(np_nan)))
    for dataset in opts.datasets:
        for algo in opts.algos:
            for metric in opts.metrics:
                if dataset in dataDict.keys() and algo in dataDict[dataset].keys() and 'ground truth' in dataDict[dataset].keys() and metric in dataDict[dataset][algo].keys() and metric in dataDict[dataset]['ground truth'].keys():
                    if(opts.mode=="global"):
                        score_pred = dataDict[dataset][algo][metric]
                        score_ref = dataDict[dataset]['ground truth'][metric]
                        if mode == "mse":
                            resDict[dataset][algo][metric] = get_mse(score_pred, score_ref, rel=False)
                        if mode == "mse_rel":
                            resDict[dataset][algo][metric] = get_mse(score_pred, score_ref, rel=True)
                    elif(opts.mode=="local"):
                        scores_ref = dataDict[dataset]['ground truth'][metric]
                        scores_pred = dataDict[dataset][algo][metric]
                    
                        if isinstance(scores_ref, dict) and isinstance(scores_ref, dict):
                            scores_ref = dict((node, scores_ref[node]) for node in scores_pred) #predicted network may have missing nodes
                            
                            if mode == "jaccard_ratio":
                                resDict[dataset][algo][metric] = get_jaccard(scores_pred, scores_ref, norm_to_random = True)
                            elif mode == "jaccard":
                                resDict[dataset][algo][metric] = get_jaccard(scores_pred, scores_ref, norm_to_random = False)
                            elif mode == "correlation":
                                resDict[dataset][algo][metric] = get_correlation(scores_pred, scores_ref)
                        else:
                            resDict[dataset][algo][metric] = np.nan
                else:
                    resDict[dataset][algo][metric] = np.nan
        
    return resDict

def get_mse(score_pred, score_ref, rel=False):
    mse = score_pred - score_ref
    if rel: 
        mse = mse / abs(score_ref)
    return mse
                
def get_jaccard(scores_pred, scores_ref, norm_to_random = True):
    N = len(scores_ref)
    Nselect = int(np.ceil(0.1*N))
    
    #np.random.seed(0)
    scores_pred = sorted(scores_pred.items(), key=lambda kv: (kv[1], np.random.rand()), reverse=True) #randomize choice between nodes with similar score
    scores_ref = sorted(scores_ref.items(), key=lambda kv: kv[1],reverse=True)
    hubs_pred = [node[0] for node in scores_pred[0:Nselect]]
    hubs_ref = [node[0] for node in scores_ref[0:Nselect]]

    jac_values = jaccard(hubs_pred, hubs_ref)
    scores = jac_values

    if norm_to_random:
        #N = min(len(scores_ref), 20000) skip larger N to reduce computation time (minimal deviation only)
        #Nselect = int(np.ceil(0.1*N))

        def prob (x ,N , Nselect):
            p =  (sp.special.comb(Nselect, x, exact=True)*sp.special.comb((N-Nselect), (Nselect-x), exact=True))/(sp.special.comb(N,Nselect, exact=True))
            return(p)

        random_jac = sum([prob(x, N, Nselect)*x/(2*Nselect- x) for x in range(Nselect+1)])
        scores = jac_values/random_jac
    
    return scores              

def get_correlation(scores_pred, scores_ref):
    if np.unique(list(scores_pred.values())).size == 1 or np.unique(list(scores_ref.values())).size == 1:
        score = np.nan
    else:
        score = stats.spearmanr(list(scores_pred.values()), list(scores_ref.values())).correlation
    return score
