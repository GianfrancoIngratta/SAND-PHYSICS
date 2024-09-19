import uproot4 as upr
import numpy as np
import pandas as pd
import ROOT

class Converter:
    def __init__(self, file_name: str, tree_name: str):
        """
        Initializes the ROOT2PandasConverter class.
        
        Parameters:
        - file_name (str): input file root.
        - tree_name (str): tree name in the root input file.
        """
        if not file_name.endswith('.root'):
            raise ValueError(f"Invalid file extension: {file_name}. The file must have a .root extension.")
        
        self.file_name = file_name
        self.tree = upr.open(file_name)[tree_name]
    
    def keys(self):
        """
        Print tree columns
        """
        print(self.tree.keys())

    def CreatePandas(self, columns: list, rename: bool = False, indices: list = None):
        """
        Create a Pandas DataFrame from the ROOT tree.
        - rename(bool): simplify columns'names
        - indices: list of colName to set as indices 
        - columns (str): name of the root branches for the output file
        Returns:
        - pd.DataFrame: A DataFrame containing the requested columns from the ROOT tree.
        
        Raises:
        - TypeError: If the resulting object is not a Pandas DataFrame.
        """
        df = self.tree.arrays(columns, library='pd')
        # Check if the result is indeed a DataFrame
        
        if not isinstance(df, pd.DataFrame):
            print("Result is not a DataFrame. Here's the tuple returned:", df)
            raise TypeError("Expected a Pandas DataFrame, but got a tuple instead.")

        if(rename):
            new_columns = {}
            for column in df.columns:
                new_columns[column] = self.GenerateNewColumn(column)

            df = df.rename(columns = new_columns)

        if(indices):
            df = df.set_index(indices)    
        
        return df
    
    def GenerateNewColumn(self, input_column: tuple):
        """
        parse input column (tuple) and return a short name
        Example:
        ('PrimariesFirstHitECAL', 'fP', 'fX') -> PrimariesFirstHitECAL_x
        ('PrimariesFirstHitECAL', 'fE', '')   -> PrimariesFirstHitECAL_t,
        """
        colName, suffix1, suffix2 = input_column
        if(suffix1=='fP'):
            if(suffix2=='fX'):
                return colName+"_x"
            elif(suffix2=='fY'):
                return colName+"_y"
            elif(suffix2=='fZ'):
                return colName+"_z"
            else:
                return colName
        elif(suffix1=='fE'):
            return colName+"_t"
        else:
            return colName

