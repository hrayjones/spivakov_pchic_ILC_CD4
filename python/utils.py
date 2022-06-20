import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.metrics import r2_score
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from decimal import Decimal
import matplotlib.ticker as ticker

def get_spearman(matrix, x, y):
    tmp_df = matrix[[x,y]].copy()
    
    tmp_df.dropna(inplace=True)
    
    corr, pval = spearmanr(tmp_df[x], tmp_df[y])
    
    return corr, pval

class pchic:
    def __init__(self, input_pchic, dropna=True, drop_off_target=True, drop_p2p=True, drop_trans_chrom=True):        
        pchic_df =  pd.read_csv(input_pchic, sep="\t", header=0)

        pchic_df["baitChr"] = "chr" + pchic_df["baitChr"].apply(str)
        pchic_df["oeChr"] = "chr" + pchic_df["oeChr"].apply(str)

        if dropna:
            pchic_df.dropna(subset=["dist"], inplace=True)

        if drop_off_target:
            pchic_df[pchic_df["baitName"] != "off_target"]

        if drop_p2p:
            pchic_df = pchic_df[pchic_df["oeName"] == "."]

        if drop_trans_chrom:
            pchic_df = pchic_df[pchic_df["baitChr"] == pchic_df["oeChr"]]

        pchic_df["OE_width"] = pchic_df["oeEnd"] - pchic_df["oeStart"]
        
        pchic_df["ID"] = pchic_df["baitChr"] + ":" + \
                                     pchic_df["baitStart"].apply(str) + "-" + \
                                     pchic_df["baitEnd"].apply(str) + "_" + \
                                     pchic_df["oeChr"] + ":" + \
                                     pchic_df["oeStart"].apply(str) + "-" + \
                                     pchic_df["oeEnd"].apply(str)

        self.pchic_df = pchic_df

        self.OE_DF = self.pchic_df[["oeChr", "oeStart", "oeEnd", "OE_width"]].drop_duplicates(subset=["oeChr", "oeStart", "oeEnd"], keep="first")

def map_counts(input_filename):
    counts_df = pd.read_csv(input_filename, sep="\t", header=None, names=["chr", "start", "stop", "ID", "counts"])
            
    counts_dict =  pd.Series(counts_df["counts"].values,index=counts_df["ID"]).to_dict()

    return counts_dict

def expression_cutoff_analysis(df_expression, 
                               feature_col="", 
                               mean_col="", 
                               title="Mean Expression of Genes with ATAC-seq overlapping PIR - Fragments", 
                               xlabel="Number of Bait-PIR Interactions", 
                               log=True, 
                               plot="",
                               plot_filename=""):

    # Get the expression data for all data that is above 0 counts
    #df_expression = df_expression[df_expression[mean_col] > 0]
    
    # Calculate the spearman correlation
    corr, pval = get_spearman(df_expression, feature_col, mean_col)

    # Create the groupby dataframe
    thresholded_expression_df = pd.DataFrame(df_expression.groupby([feature_col])[mean_col].apply(list)).reset_index().explode(mean_col)
    
    thresholded_expression_df[feature_col] = thresholded_expression_df[feature_col].apply(int)

    if plot:
        plot_expression_analysis(thresholded_expression_df, feature_col, mean_col, corr, pval, plot_filename=plot_filename, xlabel=xlabel, title=title)
    else:
        return thresholded_expression_df

def plot_expression_analysis(thresholded_expression_df, feature_col, mean_col, corr, pval, plot_filename="", xlabel="Number of Bait-PIR Interactions", title="Mean Expression of Genes with ATAC-seq overlapping PIR - Fragments"):
    with plt.style.context("default"):
        fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True, dpi=200)
        
        plt.suptitle(title, fontsize=16)
        
        max_xval = len(thresholded_expression_df[feature_col].unique())

        plt.yticks(fontsize=14)
        
        sns.boxplot(data=thresholded_expression_df, x=feature_col, y=mean_col, color="#d9d9d9", ax=axes[0])
        sns.regplot(data=thresholded_expression_df, x=feature_col, y=mean_col, color="#fdae6b", ax=axes[1])

        if max_xval >= 30:
            temp = axes[0].xaxis.get_ticklabels()
            temp = list(set(temp) - set(temp[::5]))
            for label in temp:
                label.set_visible(False)
        
        elif max_xval > 15 & max_xval <= 30:
            temp = axes[0].xaxis.get_ticklabels()
            temp = list(set(temp) - set(temp[::2]))
            for label in temp:
                label.set_visible(False)
        else:
            pass
               
        plt.yscale("log")
        plt.ylim(0,10**6)

        for axis in axes:
            axis.grid(zorder=0)
            axis.minorticks_on()
            axis.grid(axis="y", which='minor', linestyle=':', linewidth='0.5', color='black', alpha=.2)

            textstr = ' | '.join((
                'Spearman=%.3f' % corr,
                'Pval=%.2E' % pval))

            props = dict(boxstyle='round', facecolor='white', alpha=0.5)

            axis.text(2, (10**6.4), textstr, fontsize=12,
            verticalalignment='top', bbox=props)

        for ax in axes.flat:
            ax.set(xlabel=xlabel, ylabel="Gene Expression (Length Scaled TPM)")

        # Hide x labels and tick labels for top plots and y ticks for right plots.
        for ax in axes.flat:
            ax.label_outer()

        if plot_filename:
            fig.savefig(plot_filename + ".png", dpi=150)
            thresholded_expression_df.to_csv(plot_filename + ".tsv", sep="\t", index=False)
        else:
            pass

def binary_expression_cutoff_analysis(df_expression, feature_col="", mean_col=""):    
    one_df = df_expression[df_expression[feature_col] >= 1][mean_col].to_list()
    
    zero_df = df_expression[df_expression[feature_col] == 0][mean_col].to_list()

    return one_df, zero_df

def map_ID(input_filename):
    counts_df = pd.read_csv(input_filename, sep="\t", header=0, names=["Gene stable ID", "Transcript stable ID", "Gene name", "Gene stable ID version"])
            
    counts_dict =  pd.Series(counts_df["Gene name"].values,index=counts_df["Gene stable ID version"]).to_dict()

    return counts_dict

def chicago_to_PIR_BED(input_pchic, drop_off_target=True, drop_trans_chrom=True, score_col=""):        
        pchic_df =  pd.read_csv(input_pchic, sep="\t", header=0, low_memory=False)

        pchic_df["baitChr"] = "chr" + pchic_df["baitChr"].apply(str)
        pchic_df["oeChr"] = "chr" + pchic_df["oeChr"].apply(str)

        if drop_off_target:
            pchic_df[pchic_df["baitName"] != "off_target"]

        if drop_trans_chrom:
            pchic_df = pchic_df[pchic_df["baitChr"] == pchic_df["oeChr"]]

        if score_col:
            pchic_df = pchic_df[pchic_df[score_col] >= 5]
            
        pchic_df["OE_width"] = pchic_df["oeEnd"] - pchic_df["oeStart"]
        
        pchic_df["ID"] = pchic_df["baitChr"] + ":" + \
                                     pchic_df["baitStart"].apply(str) + "-" + \
                                     pchic_df["baitEnd"].apply(str) + "_" + \
                                     pchic_df["oeChr"] + ":" + \
                                     pchic_df["oeStart"].apply(str) + "-" + \
                                     pchic_df["oeEnd"].apply(str)

        OE_DF = pchic_df[["oeChr", "oeStart", "oeEnd", "ID"]].drop_duplicates(subset=["oeChr", "oeStart", "oeEnd"], keep="first")

        return pchic_df, OE_DF
