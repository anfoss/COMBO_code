from itertools import combinations
import time
import dask
from sklearn import preprocessing
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LassoCV
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, confusion_matrix, precision_score, recall_score, matthews_corrcoef, f1_score
from sklearn import metrics
from itertools import combinations
from dask.diagnostics import ProgressBar
from more_itertools import chunked


from sklearn.model_selection import StratifiedKFold


def parse_input():
    dd = pd.ExcelFile('meta/COMBO PRO Banked Samples Aliquot Clinical Data.xlsx')
    kk = []
    for x in dd.sheet_names:
        tmp = pd.read_excel('meta/COMBO PRO Banked Samples Aliquot Clinical Data.xlsx', sheet_name=x)
        tmp['loc'] = x
        kk.append(tmp)
    tot = pd.concat(kk)
    return tot


def preprocess_metadata(flt):
    # metadata
    cols = "#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2"
    cols = [x.lower() for x in list(cols)]

    meta = parse_input()
    dd = list(set(meta['COMBO Plasma Box Number']))
    meta = meta[meta['COMBO ID'].isin(flt)]
    enc = dict(zip(dd, cols))
    col_colors = meta['COMBO Plasma Box Number'].map(enc)

    col_colors.index = meta['COMBO ID']
    meta['COMBO Plasma Box Number'] = pd.Categorical(
        meta['COMBO Plasma Box Number'])

    meta['TB code'] = pd.Categorical(
        meta['TB Classification'])

    meta['code'] = meta['COMBO Plasma Box Number'].cat.codes
    meta['TB code'] = meta['TB code'].cat.codes
    meta.dropna(subset=['TB Classification'], inplace=True)
    return meta, col_colors


def get_prot_cl(std=True, nmr=False):
    X = pd.read_csv('data/protein_level_melted_merged.csv')
    X['Pr_gn'] = X['Protein.Group'] + '@' + X['Genes']
    X = pd.pivot_table(X, columns='COMBO ID', index='Pr_gn', values='value')
    meta, col_colors = preprocess_metadata(list(X))
    cls = meta[meta['TB Classification'].isin(['Confirmed TB', 'Unlikely TB'])]
    prot_cl = X[list(cls['COMBO ID'])].T
    if std:
        prot_cl = prot_cl.apply(lambda x: stats.zscore(x), axis=1)
    elif nmr: 
        prot_cl = pd.DataFrame(preprocessing.minmax_scale(prot_cl, feature_range=(0, 1), axis=1, copy=True), index=prot_cl.index, columns=prot_cl.columns)
    else:
        pass
    cls = dict(zip(cls['COMBO ID'], cls['TB Classification']))
    y = [cls[x] for x in list(prot_cl.index)]
    #print(len([x for x in y if x=='Confirmed TB']), len(y))
    y = [1 if x == 'Confirmed TB' else 0 for x in y]
    return prot_cl, y


def get_coef_lasso(prot_cl, y):
    reg = LassoCV(cv=StratifiedKFold(20), max_iter=10000, tol=0.0001, random_state=1).fit(prot_cl, y)
    # reg = LassoCV(cv=20, max_iter=10000, tol=0.001, random_state=1).fit(prot_cl, y)

    print("Best alpha using built-in LassoCV: %f" % reg.alpha_)
    print("Best score using built-in LassoCV: %f" % reg.score(prot_cl, y))
    coef = pd.Series(reg.coef_, index=prot_cl.columns)
    metr = []
    imp_coef = coef.sort_values()
    imp_coef = imp_coef[imp_coef!=0]
    imp_coef = imp_coef.to_frame()
    imp_coef['abs'] = imp_coef[0].abs()
    imp_coef = imp_coef.sort_values(by='abs', ascending=False)
    imp_coef.to_csv('output/imp_coef_lasso.csv')
    imp_coef.drop(0, axis=1, inplace=True)
    imp_coef['rank_lasso'] = imp_coef['abs'].rank(ascending=False)
    return imp_coef
  

def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if "log_time" in kw:
            name = kw.get("log_name", method.__name__.upper())
            kw["log_time"][name] = int((te - ts) * 1000)
        else:
            print("%r  %2.2f ms" % (method.__name__, (te - ts) * 1000))
        return result
    return timed


def evaluate_combination(X_train, X_test, y_train, y_test, feature_names, cmb, best_threshold):
    X_train_sub = X_train[:, cmb]
    X_test_sub = X_test[:, cmb]

    clf = LogisticRegression(penalty='l2', solver='liblinear')
    clf.fit(X_train_sub, y_train)

    y_prob = clf.predict_proba(X_test_sub)[:, 1]
    fpr, tpr, _ = metrics.roc_curve(y_test, y_prob, pos_label=1)
    prec70sens = tpr[np.abs(fpr - 0.3).argmin()]
    auc = metrics.auc(fpr, tpr)

    if prec70sens > best_threshold:
        combo_names = ','.join(feature_names[list(cmb)])
        return [combo_names, prec70sens, auc]
    return None


def evaluate_combinations_batch(X_train, X_test, y_train, y_test, feature_names, cmb_list, best_threshold):
    results = []
    for cmb in cmb_list:
        res = evaluate_combination(X_train, X_test, y_train, y_test, feature_names, cmb, best_threshold)
        if res is not None:
            results.append(res)
    return results


def get_combi_parallel(X_train, X_test, y_train, y_test, feature_names, n=3, best_threshold=0.9, chunk_size=10000):
    all_combs = list(combinations(range(X_train.shape[1]), n))
    tasks = [dask.delayed(evaluate_combinations_batch)(
        X_train, X_test, y_train, y_test, feature_names, list(chunk), best_threshold)
        for chunk in chunked(all_combs, chunk_size)]
    batched_results = dask.compute(*tasks, threads_per_worker=1, processes=True, n_workers=61)
    results = [r for batch in batched_results for r in batch]
    df = pd.DataFrame(results, columns=['comb', 'prec70sens', 'auc'])
    df.to_csv(f'output/allAUC_lasso_{n}.csv', sep='\t', index=False)
   


def main():
    ProgressBar().register()
    prot_cl, y = get_prot_cl(True, False)
    lassocoef = pd.read_csv('meta\imp_coef_lasso.csv')
    coefs = lassocoef[lassocoef['abs']!=0]
    selected_features = prot_cl[coefs['Pr_gn']]
    X = selected_features.to_numpy()
    y_array = np.asarray(y)
    feature_names = selected_features.columns.to_numpy()
    X_train, X_test, y_train, y_test = train_test_split(
        X, y_array, test_size=0.25, stratify=y_array, random_state=1
    )
    print(X_train.shape)
    #all_results = []
    for i in range(1, 7):
        get_combi_parallel(X_train, X_test, y_train, y_test, feature_names, n=i)
        #all_results.append(df)

    # final_df = pd.concat(all_results, ignore_index=True)
    #get_combi_parallel(X_train, X_test, y_train, y_test, feature_names, n=6)


if __name__ == "__main__":
    main()
