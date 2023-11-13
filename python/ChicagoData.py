"""ChicagoData.py
Import and filter Chicago results
"""

from unittest import skip
import pandas as pd
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import pybedtools
import os
from os import makedirs, error

class ChicagoData(object):
    """Import CHiCAGO data
    """
    def __init__(self,
                 filename: str,
                 drop_off_target_bait: bool = True,
                 drop_off_target_oe: bool = True,
                 drop_trans_chrom: bool = True,
                 remove_p2p: bool= True,
                 features_to_count: dict = {},
                 gene_expression: str = "",
                 nonzero_expression: bool = True,
                 dropna_expression: bool = True,
                 output_dir: str = "",
                 output_basename: str = "",
                 chain_file: str = ""
                 ):
        """Initialize the object

        This object will perform the following methods in order upon initialization:
        
        # Read file into DF
        self._read_file_()
        
        # Format the DF
        self._format_file_()
        
        # Filter the formatted DF
        self._filter_file_()
        
        # Get the PIR df
        self._get_PIR_df_()
        
        # Get the bait df
        self._get_bait_df_()
        
        # Get the combined df
        self._get_combined_df_()  
        
        # Get the feature counts per PIR and map to hg19
        self._get_feature_counts_()
        
        # Import the gene expression matrix
        self._import_gene_counts_()
        
        # Map the feature counts to genes in the expression matrix
        self._map_feature_counts_to_genes_()
        
        # Map the ABC count or CHiCAGO count as a column in the gene expression matrix
        self._map_ABC_counts_to_genes_()
        
        # Filter the gene expression matrix file
        self._filter_expression_()
        
        # Get the PIR count v mean gene expression
        self._get_PIR_count_v_mean_()
        
        # Write the newly formatted CHiCAGO data to a file
        self._write_new_chicago_data_()
        
        Args:
            filename (str): Input filename
            drop_off_target_bait (bool, optional): Drop off target baits. Defaults to True.
            drop_off_target_oe (bool, optional): Drop off target OE. Defaults to True.
            drop_trans_chrom (bool, optional): Drop trans chromosomal interactions. Defaults to True.
            remove_p2p (bool, optional): _description_. Defaults to True.
            features_to_count (dict, optional): _description_. Defaults to {}.
            gene_expression (str, optional): _description_. Defaults to "".
            nonzero_expression (bool, optional): _description_. Defaults to True.
            dropna_expression (bool, optional): _description_. Defaults to True.
            output_dir (str, optional): _description_. Defaults to "".
            output_basename (str, optional): _description_. Defaults to "".
        """
        # Set filename to the input filename
        self.filename = filename
        # Set whether to drop off target baits
        self.drop_off_target_bait = drop_off_target_bait
        # Set whether to drop off target oe
        self.drop_off_target_oe = drop_off_target_oe
        # Set whether to drop transchromosomal interactions
        self.drop_trans_chrom = drop_trans_chrom
        # Remove promoter to promoter interactions
        self.remove_p2p = remove_p2p
        # Map feature counts to PIR if provided
        self.features_to_count = features_to_count
        # Import the gene expression matrix
        self.gene_expression = gene_expression
        # Only keep non-zero expression
        self.nonzero_expression = nonzero_expression
        # Drop na Values
        self.dropna_expression = dropna_expression
        # Set the output directory
        self.output_dir = output_dir
        # Set the output basename
        self.basename = output_basename
        # Set the chain file
        self.chainfile = chain_file
        
        # Read file into DF
        self._read_file_()
        
        # Format the DF
        self._format_file_()
        
        # Filter the formatted DF
        self._filter_file_()
        
        # Get the PIR df
        self._get_PIR_df_()
        
        # Get the bait df
        self._get_bait_df_()
        
        # Get the combined df
        self._get_combined_df_()
        
        # Get the feature counts per PIR        
        self._get_feature_counts_()
        
        # Import the gene expression matrix
        self._import_gene_counts_()
        
        # Map the feature counts to genes in the expression matrix
        self._map_feature_counts_to_genes_()
        
        # Filter the gene expression matrix file
        self._filter_expression_()
        
        # Get the PIR count v mean gene expression
        self._get_PIR_count_v_mean_()
        
        # Write the newly formatted CHiCAGO data to a file        
        self._write_new_chicago_data_()        
                    
    def _read_file_(self):
        """Read in original file
        """
        # Read in original file and save
        self.input_df =  pd.read_csv(self.filename, sep="\t", header=0, low_memory=False)
    
    def _write_new_chicago_data_(self):
        output_dir = os.path.join(self.output_dir, "modified_chicago")
        
        output_filename = os.path.join(output_dir, f"{self.basename}_modified.tsv")
        
        self._get_dir_(output_dir)

        self.df.to_csv(output_filename, sep="\t", index=False)
        
    def _format_file_(self):
        """Format CHICAGO file
        """
        # Create a copy of the raw input to be manipulated
        df = self.input_df.copy()
        
        # Format the chromosome names
        df["baitChr"] = "chr" + df["baitChr"].apply(str)
        df["oeChr"] = "chr" + df["oeChr"].apply(str)
        
        # Create an ID column for the OE
        df["oe_interval_ID"] = df["oeChr"] + ":" + \
                   df["oeStart"].apply(str) + "-" + \
                   df["oeEnd"].apply(str)

        # Create a bait ID column
        df["bait_interval_ID"] = df["baitChr"] + ":" + \
                   df["baitStart"].apply(str) + "-" + \
                   df["baitEnd"].apply(str)

        # Create an interaction ID from the bait interval ID and the OE ID
        df["interaction_ID"] = df["bait_interval_ID"] + "_" + df["oe_interval_ID"].apply(str)
        
        # Find the bait ID from CHiCAGO
        self.bait_ID = df["baitID"].unique()
        
        # Find the unique OE ID
        self.oe_ID = df["oeID"].unique()

        # Find the unique bait interval IDs
        self.bait_interval_ID = df["bait_interval_ID"].unique()
        
        # Find the unique OE ID
        self.oe_interval_ID = df["oe_interval_ID"].unique()
        
        self.df = df
        
    def _filter_file_(self):
        """Filter the formatted CHICAGO results
        """
        # Drop the off target baits
        if self.drop_off_target_bait:
            self.df[self.df["baitName"] != "off_target"]

        # Drop the off target OE names
        if self.drop_off_target_oe:
            self.df[self.df["oeName"] != "off_target"]

        # Drop the trans chromosomal interactions
        if self.drop_trans_chrom:
            self.df = self.df[self.df["baitChr"] == self.df["oeChr"]]

        # Drop promoter to promoter interactions
        if self.remove_p2p:            
            self.df = self.df[~self.df.oe_interval_ID.isin(self.bait_interval_ID)]
        
    def _get_PIR_df_(self):
        """Get a DF of all PIR interactions
        """
        self.pir_df = self.df[["oeChr", "oeStart", "oeEnd", "oe_interval_ID"]].drop_duplicates(subset=["oeChr", "oeStart", "oeEnd"], keep="first")
        
        self.PIR_bt = pybedtools.BedTool.from_dataframe(self.pir_df)

    def _get_bait_df_(self):
        """Get a DF of baits
        """
        self.bait_df = self.df[["baitChr", "baitStart", "baitEnd", "bait_interval_ID"]].drop_duplicates(subset=["baitChr", "baitStart", "baitEnd"], keep="first")
        
    def _get_combined_df_(self):
        """Get a comined DF
        """
        tmp_df = self.pir_df.copy()
        tmp_df2 = self.bait_df.copy()

        tmp_df.columns = ["Chr", "Start", "Stop", "ID"]
        tmp_df2.columns = ["Chr", "Start", "Stop", "ID"]

        self.unique_features = pd.concat([tmp_df, tmp_df2])
    
    def _get_feature_counts_(self):
        """Get the counts of features that overlap PIRs and map them back to the pcHiC interaction
        """
        for file, tag in self.features_to_count.items():
            print(f"Importing {file} : Column will be saved as {tag}")
            output_intersection_dir = os.path.join(self.output_dir, "PIR_intersection")

            output_intersection_fname = os.path.join(output_intersection_dir, f"{self.basename}_PIR_intersect_{tag}.bed")
            
            # Import and convert the featuers to a bedtools
            feature_bt = pybedtools.BedTool(file)
            
            # Intersect the features to get counts
            feature_counts = self.PIR_bt.intersect(feature_bt, c=True)
            
            # Intersect the features to get overlaps (true intersections)
            feature_intersection = self.PIR_bt.intersect(feature_bt)
            
            feature_intersection_sort = feature_intersection.sort()
            
            # Convert bedtools to pandas
            feature_counts_df = feature_counts.to_dataframe()
                        
            # Convert bedtools to dataframe
            feature_intersection_df = feature_intersection_sort.to_dataframe()
                        
            # Create a dictionary of counts 
            counts_dict =  pd.Series(feature_counts_df["score"].values,index=feature_counts_df["name"]).to_dict()

            # Map the counts back to the CHICAGO dataframe
            self.df[tag] = self.df["oe_interval_ID"].map(counts_dict)
            
            self._get_dir_(output_intersection_dir)
            
            feature_intersection_df[["chrom", "start", "end"]].to_csv(output_intersection_fname, sep="\t", index=False, header=False)
            
            # Liftover the PIR files (here from hg38 to hg19)
            if self.chainfile:
                print(f"Lifting over {tag} PIR intersections using {self.chainfile}") 
                
                liftover_dir = os.path.join(self.output_dir, "PIR_intersection_liftover")
                
                # Perform liftOver on the bedtools object 
                lifted = feature_intersection_sort.liftover(self.chainfile, unmapped=None) 
                
                # Converted lifted coordinates to data frame
                lifted_feature_intersection_df = lifted.to_dataframe()
                
                # Save the results
                liftover_intersection_fname = os.path.join(liftover_dir, f"{self.basename}_PIR_intersect_{tag}_hg19.bed")
                lifted_feature_intersection_df[["chrom", "start", "end"]].to_csv(liftover_intersection_fname, sep="\t", index=False, header=False)

    def _import_gene_counts_(self):
        self.gene_counts = pd.read_csv(self.gene_expression, sep="\t", header=0, names=["GeneName", "Expression"])
        
    def _map_feature_counts_to_genes_(self):
        self.gene_counts["enhancer_count"] = self.gene_counts["GeneName"].map(self.df.groupby(["baitName"]).count()["baitChr"])
        
        for _, tag in self.features_to_count.items():
            self.gene_counts[f"{tag}_count"] = self.gene_counts["GeneName"].map(self.df.groupby(["baitName"]).sum()[tag])

    def _get_dir_(self, dir: str, permissions=0o0775, exist_ok : bool=True):
        """Makes a directory at the given location
        Args:
            dir (str): Path of the directory
            permissions ([type], optional): Permissions of directory. Defaults to 0o0775.
            exist_ok (bool, optional): If True, program will continue if directory exists. Defaults to True.
        Returns:
            str: Absolute path to the created directory
            
        Example:
        
        >>> output_dir = get_dir("./output/")
        """
        try:
            makedirs(dir, mode=permissions)
        except error:
            if not exist_ok:
                raise

    def _filter_expression_(self):
        output_analysis_dir = os.path.join(self.output_dir, "expression_matrix")

        output_analysis_fname_unfiltered = os.path.join(output_analysis_dir, f"{self.basename}_unfiltered_expression_matrix.tsv")
        output_analysis_fname_filtered = os.path.join(output_analysis_dir, f"{self.basename}_filtered_expression_matrix.tsv")

        self._get_dir_(output_analysis_dir)

        self.gene_counts.to_csv(output_analysis_fname_unfiltered, sep="\t", index=False)

        if self.nonzero_expression:
            self.gene_counts = self.gene_counts[self.gene_counts["Expression"] > 0]

        if self.dropna_expression:
            self.gene_counts = self.gene_counts.dropna()

        self.gene_counts["GeneName_MeanExpression"] = self.gene_counts["GeneName"] \
            + " " + self.gene_counts["Expression"].apply(str)
            
        self.gene_counts.to_csv(output_analysis_fname_filtered, sep="\t", index=False)
        
    def _calculate_spearman_(self):
        self.corr = []
        output_analysis_dir = os.path.join(self.output_dir, "correlation_analysis")

        output_analysis_fname = os.path.join(output_analysis_dir, f"{self.basename}_correlation_stats.tsv")
             
        for _, tag in self.features_to_count.items():
            tmp_df = self.gene_counts[[f"{tag}_count", "Expression"]].copy()
            
            s_corr, s_pval = spearmanr(tmp_df[f"{tag}_count"], tmp_df["Expression"])
            
            p_corr, p_pval = pearsonr(tmp_df[f"{tag}_count"], tmp_df["Expression"])

            self.corr.append([s_corr, s_pval, p_corr, p_pval, tag])

        self.corr_df = pd.DataFrame(self.corr, columns=["Spearman_corr", "Spearman_pval", "Pearson_corr", "Pearson_pval", "Feature"])
            
        self._get_dir_(output_analysis_dir)

        self.corr_df.to_csv(output_analysis_fname, sep="\t", index=False)

    def _get_PIR_count_v_mean_(self):
        """Create the dataframe of the number of features overlapping PIRs by mean gene expression
        
        First you groupby the feature counts column. Then you find all of the mean gene expression values
        associated with the number of features overlapping a PIR. The output is a dataframe that can be plotted
        """
    
        for _, tag in self.features_to_count.items():
            col_name = f"{tag}_count"
            output_analysis_dir = os.path.join(self.output_dir, "expression_analysis")

            output_analysis_fname = os.path.join(output_analysis_dir, f"{self.basename}_{tag}.tsv")
            
            # Create the dataframe
            self.pir_count_v_mean = pd.DataFrame(self.gene_counts.groupby([col_name])["GeneName_MeanExpression"].apply(list)).reset_index().explode("GeneName_MeanExpression")
        
            self.pir_count_v_mean[col_name] = self.pir_count_v_mean[col_name].apply(int)
            
            self.pir_count_v_mean[["Gene_Name", "Mean_Gene_Expression"]] = self.pir_count_v_mean["GeneName_MeanExpression"].str.split(" ", n=2, expand=True)
            
            self.pir_count_v_mean.drop("GeneName_MeanExpression",axis=1, inplace=True)
            
            self.pir_count_v_mean["Mean_Gene_Expression"] = self.pir_count_v_mean["Mean_Gene_Expression"].apply(float)
            
            self._get_dir_(output_analysis_dir)

            self._calculate_spearman_()
            
            self.pir_count_v_mean.to_csv(output_analysis_fname, sep="\t", index=False)