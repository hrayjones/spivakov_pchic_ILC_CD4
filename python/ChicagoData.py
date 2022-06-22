"""ChicagoData.py
Import and filter Chicago results
"""

import pandas as pd

class ChicagoData(object):
    """Import CHICACO data
    """
    def __init__(self,
                 filename: str,
                 dropna: bool = True,
                 drop_off_target_bait: bool = True,
                 drop_off_target_oe: bool = True,
                 drop_trans_chrom: bool = True,
                 score_col: str = None,
                 score_val: int = 5,
                 remove_p2p: bool= True
                 ):
        """Initialize the object

        Args:
            filename (str): Input CHICAGO txt file
            dropna (bool, optional): Drop the interactions with NA. Defaults to =True.
            drop_off_target (bool, optional): Drop off target interactions. Defaults to =True.
            drop_trans_chrom (bool, optional): Drop transchromosomal interactions. Defaults to =True.
        """
        # Set filename to the input filename
        self.filename = filename
        # Set whether to drop na interactions
        self.dropna = dropna
        # Set whether to drop off target baits
        self.drop_off_target_bait = drop_off_target_bait
        # Set whether to drop off target oe
        self.drop_off_target_oe = drop_off_target_oe
        # Set whether to drop transchromosomal interactions
        self.drop_trans_chrom = drop_trans_chrom
        # Score column name for filtering
        self.score_col = score_col
        # Score value threshold
        self.score_val = score_val
        # Remove promoter to promoter interactions
        self.remove_p2p = remove_p2p
        
        # Read file into DF
        self._read_file_()
        
        # Format the DF
        self._format_file_()
        
        # Filter the formatted DF
        self._filter_file_()

        # Get the PIR df
        self._get_PIR_df_()
        
    def _read_file_(self):
        """Read in original file
        """
        # Read in original file and save
        self.input_df =  pd.read_csv(self.filename, sep="\t", header=0, low_memory=False)
    
    def _format_file_(self):
        """Format CHICAGO file
        """
        # Create a copy of the raw input to be manipulated
        df = self.input_df.copy()
        
        # Format the chromosome names
        df["baitChr"] = "chr" + df["baitChr"].apply(str)
        df["oeChr"] = "chr" + df["oeChr"].apply(str)
        
        # Create an ID column that can be used to track the intervals
        df["ID"] = df["baitChr"] + ":" + \
                   df["baitStart"].apply(str) + "-" + \
                   df["baitEnd"].apply(str) + "_" + \
                   df["oeChr"] + ":" + \
                   df["oeStart"].apply(str) + "-" + \
                   df["oeEnd"].apply(str)

        # Set the variable to the formatted df
        self.df = df
    
    def _filter_file_(self):
        """Filter the formatted CHICAGO results
        """
        # There should not be any NA values in this column, but we will double check
        if self.dropna:
            self.df.dropna(subset=["dist"], inplace=True)

        # Drop the off target baits
        if self.drop_off_target_bait:
            self.df[self.df["baitName"] != "off_target"]

        # Drop the off target OE names
        if self.drop_off_target_oe:
            self.df[self.df["oeName"] != "off_target"]

        # Drop the trans chromosomal interactions
        if self.drop_trans_chrom:
            self.df = self.df[self.df["baitChr"] == self.df["oeChr"]]

        # Filter the specific score column by a specific value
        if self.score_col:
            self.df = self.df[self.df[self.score_col] >= self.score_val]
            
        # Drop promoter to promoter interactions
        if self.remove_p2p:
            self.df = self.df[self.df.oeName == "."]
        
    def _get_PIR_df_(self):
        """Get a DF of all PIR interactions
        """
        self.pir_df = self.df[["oeChr", "oeStart", "oeEnd", "OE_width"]].drop_duplicates(subset=["oeChr", "oeStart", "oeEnd"], keep="first")

    def write_PIR_bed(self, output_filename): 
        """Write PIRs from the filtered CHICAGO results to a bed file
        """       
        self.pir_df.to_csv(output_filename, sep="\t", header=False, index=False)
