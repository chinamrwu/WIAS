#!/home/guomics/anaconda3/bin/python
#coding=utf-8
import pandas as pd
import os, sys
from pandas import Series, DataFrame
from numpy import NaN
import numpy as np
import datetime
import time
from sklearn.linear_model import LinearRegression
from quantile_norm import *


"""
1.Prepare the peptide matrix from TRIC, where each row represents a peptide and each column a sample. 55.97% values at 
peptide are NA.
2.Log2 transform the peptide matrix. NA values remain as NA.
3.Quantile normalization based on samples. (pppb_pep2.txt)
4.Technical imputation (Qing please describe) Missing values?  (pppb_pep3.txt, missing rate is 42.24%)

5.Replace all NA values by 0.
6.Batch correction based on the six MS instrument using the Python package Combat. 
7.Replace all values, which are original NA values from Step 4 by NA. (180519pppb_pep4.txt, missing rate is 42.24%)
8.Calculate the mean expression of each row to rank the peptide precursor intensity.
9.Order peptide precursors in protein group, first by number of NAs (Intensity equal to 0, ascending) and then by the
 mean expression (descending)
10.Keep no more than top 3 peptides (# of NAs ascending and order by mean expression descending) for each protein group.
 (2018-05-21pppb_pep5.2.txt, missing rate is 34.14%).
11.To impute some missing values at protein level from multiple peptides quantified, we built linear model for 
step 10 matrix, always use the top protein groups as dependant variable (y) and chose values greater than 0 as multiple 
independent variables (X) by decreasing order of priority. Moreover, impute y using a linear combination of X. 
Only the regression coefficient with P value<=0.05 and R2 > 0.36 were accepted and rounded up to two decimal places. 
Then we used the intensity value from the top1 peptide precursor to represent the protein intensity 
12.Keep only prototypic peptides and their corresponding proteins. (2018-05-21pppb_prot.txt) missing rate is 33.03% 
(including all 1566 samples)
"""


def open_txt(filename, sep):
    print("start reading data from %s" % filename)
    try:
        data_frame = pd.read_table(filename, sep=sep)
    except:
        print("Error occured when reading data!")
        sys.exit()
    data_frame.replace("NA", np.NaN, inplace=True)
    data_frame.replace("Na", np.NaN, inplace=True)
    data_frame.replace("na", np.NaN, inplace=True)

    na_num = data_frame.isnull().sum().sum()
    rows, cols = data_frame.shape
    total = rows * (cols - 2)
    print('NA ratio = %f'%(na_num / float(total)))
    columns = list(data_frame.columns)
    if len(columns) < 3:
        print(columns[0])
        print("error in data's columns, please check data column delimiter or parameter")
        sys.exit()
    columns[1] = 'prot'
    data_frame.columns = columns
    print("reset column_2 name = prot")
    return data_frame, na_num


def log_2(data_frame):
    print("start taking log2 on input data")
    data_frame.replace("NA", np.NaN, inplace=True)
    columns = data_frame.columns.size
    data_frame.iloc[:, 2: columns] = pd.DataFrame(data_frame.iloc[:, 2: columns], dtype=np.float64)
    data_frame.replace(0.0, NaN, inplace=True)
    format = lambda x: np.log2((x))
    data_frame.iloc[:, 2: columns] = data_frame.iloc[:, 2: columns].applymap(format)
    return data_frame


def quantile_normalization(data_frame):
    print("#" * 50)
    print("start quantile normalization")
    data_arr = data_frame.iloc[:, 2:data_frame.columns.size].copy()
    normed = quantile_normalize(data_arr)
    data_frame.iloc[:, 2:data_frame.columns.size] = normed[:, :]
    return data_frame


def read_replicate(sample_doc, names, rep_type=1, rep_header=True, sep='\t'):
    pairs = []
    with open(sample_doc, 'r') as f:
        if rep_header:
            _ = f.readline()
        for line in f:
            cols = line.strip().split(sep)
            if len(cols) != 2:
                print('- '*30)
                print("Warning: Check sample-id pair in %s, item=%d"%(sample_doc, len(cols)))
                continue
            else:
                pairs.append(cols)

    if len(names) != len(pairs):
        print('#' * 30)
        print("Warning: sample numbers of peptides not equal with that in %s" % sample_doc)

    if int(rep_type) == 1:
        samples_wrt_id = {pair[0]: pair[1] for pair in pairs}
        serial_no_sample = {p[1]:[] for p in pairs}
        for k, v in samples_wrt_id.items():
            serial_no_sample[v].append(k)
    else:
        samples_wrt_id = {pair[0]: pair[1][:-1] for pair in pairs}
        serial_no_sample = {p[1][:-1] :[] for p in pairs}
        for k, v in samples_wrt_id.items():
            serial_no_sample[v[:-1]].append(k)

    return samples_wrt_id, serial_no_sample


def technical_imputation(normed_frame, sample_doc, rep_type=1, header=True, sep='\t'):
    print("starting to do technical imputation")
    names = list(normed_frame.columns[2:])
    sample_wrt_int, serial_no_sample = read_replicate(sample_doc, names, rep_type, header, sep)
    while len(names) > 0:
        try:
            name = sample_wrt_int[names[0]]
            merge_name = serial_no_sample[name]
        except:
            print("ERROR !!!!!! one sample name not in %s !" % sample_doc)
            # merge_name = names[0]
            merge_name = []
            # raise NameError
        if len(merge_name) == 1:
            names.remove(names[0])
            continue
        
        merge_arr = np.array(normed_frame[merge_name], dtype=np.float64)
        row_mean = np.nanmean(merge_arr, axis=1)
        inds = np.where(np.isnan(merge_arr))
        merge_arr[inds] = np.take(row_mean, inds[0])
        normed_frame[merge_name] = merge_arr

        # delete names from names_list
        for col in merge_name:
            names.remove(col)

    print("imputation done")
    return normed_frame


# ranking
def ranking_the_protein(new_pep_expr, group_num=3):
    print("#" * 50)
    print("start ranking by prot, NA number, mean expression")
    new_pep_expr.replace(0.0, NaN, inplace=True)
    pep_2_frame = new_pep_expr.iloc[:, 0: 3].copy()
    nrows, ncols = new_pep_expr.shape

    ## build df as ranking parameters. make na to negative for co_sort with mean expression.
    for i in range(0, nrows):
        NA_number = ncols - new_pep_expr.iloc[i, :].count()
        pep_1 = 0.0 - (new_pep_expr.iloc[i, 2: ncols].sum()) / (ncols - NA_number)
        pep_2 = (new_pep_expr.iloc[i, 2: ncols].sum()) / ncols
        pep = np.array([NA_number, pep_1, pep_2])
        pep_2_frame.iloc[i] = pep

    ##concat dataframe
    print("computing NAs and average expression done")
    pep_2_frame.columns = ["NAs", "ave_expr1", "ave_expr2"]
    pep_expr = pd.concat([new_pep_expr, pep_2_frame], axis=1)
    pep_expr = pep_expr.drop_duplicates()

    ##ranking, and drop auxiliary columns
    sorted_expr = pep_expr.sort_values(["prot", "NAs", "ave_expr1"], axis=0, ascending=True)
    ranked_expr = sorted_expr.drop(["NAs", "ave_expr1", "ave_expr2"], axis=1)

    # only choose protein name start with "1/" and top3
    wanted_rows = ranked_expr['prot'].str.startswith("1/")
    rows_selected = ranked_expr[wanted_rows]
    prot_top3 = rows_selected.groupby('prot').head(group_num).reset_index(drop=True)

    print("Ranking is done !")
    return prot_top3


def doing_lr(row1_series, row2_series):
    row1, row2 = np.array(row1_series[2:], dtype=np.float64), np.array(row2_series[2:], dtype=np.float64)
    # print(row1)
    row1_notna, row2_notna = ~np.isnan(row1), ~np.isnan(row2)
    #print(len(row1_notna), len(row2_notna))
    both = row1_notna & row2_notna
    #print(len(both))

    if np.sum(both) > 0:
        # reshape to required shape in sklearn
        row1_both = np.reshape(row1[both], [len(row1[both]), 1])
        row2_both = np.reshape(row2[both], [len(row2[both]), 1])

        model = LinearRegression()
        model.fit(row2_both, row1_both)
        score = model.score(row2_both, row1_both)

        # change score to 0.64, 2018-11-06
        if score < 0.64:
            return row1_series
        else:
            # print(np.sum(np.isnan(row1)))
            for i in range(len(row1)):
                # after human_curated, there're some nan should be replaced.
                if np.isnan(row1[i]) and not np.isnan(row2[i]):
                    row1[i] = model.predict(row2[i])

                    # added two lines, 2018-11-06
                    if row1[i] < 0.0 or row1[i] > 25:
                        row1[i] = np.NaN
            row1_series.iloc[2:] = row1
            # print(np.sum(pd.isnull(row1_series.iloc[2:])))
            return row1_series
    else:
        return row1_series


def replace_nas_lr(ranked_protein_df):
    """
    :param ranked_protein_df: must have 'prot' and 'tg' columns at begin.
        value should start from the third columns.
    :return:
    """
    rows, cols = ranked_protein_df.shape
    protein_hash, prot_names = {}, ranked_protein_df['prot']

    # doing lr
    index = []
    for i in range(rows):
        prot = prot_names.iloc[i]
        if prot not in protein_hash:
            protein_hash[prot] = i
            index.append(True)
        else:
            # print(i, prot_names.iloc[i])
            row1 = ranked_protein_df.iloc[protein_hash[prot]]
            row2 = ranked_protein_df.iloc[i]
            row1, row2 = row1.copy(), row2.copy()
            new_row1 = doing_lr(row1, row2)
            ranked_protein_df.iloc[protein_hash[prot]] = new_row1
            index.append(False)
        if i % 200 == 0:
            print("linear regression: total rows=%d, %d rows"%(rows, i+1))
    unique_prot = ranked_protein_df[index]
    return unique_prot


def save_result(df, output_path):
    if output_path and '.txt' in output_path:
        df.to_csv(output_path, header=True, index=False, sep='\t')
    elif output_path and '.csv' in output_path:
        df.to_csv(output_path, index=False)
    else:
        print('saving to present dir, file_name=dia_preprocess_output.csv')
        df.to_csv('dia_preprocess_output.csv', index=False)


def total_interface(data_doc, sep='\t', 
                    tech_rep_file=None, rep_type=1, rep_sep='\t', rep_header=True,
                    output_path='dia_preprocess_output.csv', log2=True, quantile_norm=True,
                    combat=False, ranking=True, group_number=3, linear_reg=True):
    """
    User interface for dia_preprocessing. Default setting is log2, quantile_norm,
    then get the top3 rows in every prot group(ranking by Na, mean expression),
    lastly, using linear regression to fill in NA
    :param data_doc: the absolute path or relative path. data should organized as follow:
    tg  prot    sample1 sample2
    xxx yyyy    100.0   3000.0
    
    :sep: columns' delimiter of data_doc

    :tech_rep_file：sample replication name and its group_id file. Must have same sample name with data_doc's header.
        otherwise, error may occur.
    :rep_type: organization of technical_replicate file.
            if type=1 means sample_id is "10", otherwise means "A10a"
    :rep_sep: separation symbol of technical_replicate file. Default setting is '\t'.
    :rep_header: whether there is header in technical_replicate file.

    :param output_path: strongly recommend to explicitly type in output file_name in present working directory.
        defaultly, result will save into "dia_pre_output.csv" in present working directory

    :param log2:
    :param quantile_norm:
    :param combat:
    :param group_number:    rank结果后选择某个蛋白的前几个peptide.
    :return:  The result will write into output path, with same organization as input.
        Typically result will have fewer lines because we only select most wanted row for every prot.
    """
    df, ori_na_num = open_txt(data_doc, sep)
    if log2:
        df = log_2(df)
    #df.to_csv('log2.csv', index=False)

    if quantile_norm:
        df = quantile_normalization(df)
    #df.to_csv("quantile.csv", index=False)

    if combat:
        df.replace(NaN, 0.0, inplace=True)
        df = combat(df)

    if tech_rep_file:
        df = technical_imputation(df, tech_rep_file, rep_type, rep_header, rep_sep)
    #df.to_csv("tech_impu.csv", index=False)

    if ranking:
        df = ranking_the_protein(df, group_number)
    #df.to_csv("rank.csv", index=False)
    #df = pd.read_csv('ranking1.csv')
    if linear_reg:
        df = replace_nas_lr(df)

        new_na_num = df.isnull().sum().sum()
        rows, cols = df.shape
        total = rows * (cols - 2)
        perct = new_na_num / float(total)
        print('result NA ration=%f'%perct)

        filled = ori_na_num - new_na_num
        print('Totally, %d NAs have been filled using linear regression'%(filled))
    save_result(df, output_path)

    print(' *' * 50)
    print("All preprocess is done!!!!")

# about speed: For a input file with 43000 rows and 1000 columns,
#               may spend 40 mins to 50 mins.
if __name__ == '__main__':
    # total_interface('./peptides.txt', '\t',
    #             tech_rep_file='technical_rep.txt', rep_type=1, rep_sep='\t',rep_header=False,
    #             output_path='output.csv')
    total_interface("/work/output600/peptides_20181115.txt",'\t',tech_rep_file=None, rep_type=int(1),rep_sep='\t',rep_header=False,output_path='TPD_Win600_protMatrix_20181115.csv')
