import uproot4 as upr
import numpy as np
import pandas as pd
import ROOT

# class Converter:
#     def __init__(self, file_name: str, tree_name: str):
#         """
#         Initializes the ROOT2PandasConverter class.
        
#         Parameters:
#         - file_name (str): input file root.
#         - tree_name (str): tree name in the root input file.
#         """
#         if not file_name.endswith('.root'):
#             raise ValueError(f"Invalid file extension: {file_name}. The file must have a .root extension.")
        
#         self.file_name = file_name
#         self.tree = upr.open(file_name)[tree_name]
    
#     def keys(self):
#         """
#         Print tree columns
#         """
#         print(self.tree.keys())

#     def CreatePandas(self, columns: list, rename: bool = False, indices: list = None):
#         """
#         Create a Pandas DataFrame from the ROOT tree.
#         - rename(bool): simplify columns'names
#         - indices: list of colName to set as indices 
#         - columns (str): name of the root branches for the output file
#         Returns:
#         - pd.DataFrame: A DataFrame containing the requested columns from the ROOT tree.
        
#         Raises:
#         - TypeError: If the resulting object is not a Pandas DataFrame.
#         """
#         df = self.tree.arrays(columns, library='pd')
#         # Check if the result is indeed a DataFrame
        
#         if not isinstance(df, pd.DataFrame):
#             print("Result is not a DataFrame. Here's the tuple returned:", df)
#             raise TypeError("Expected a Pandas DataFrame, but got a tuple instead.")

#         if(rename):
#             new_columns = {}
#             for column in df.columns:
#                 new_columns[column] = self.Tuple2String(column)

#             # check if df has multi index structure
#             if isinstance(df.index, pd.MultiIndex):
#                 print("DataFrame has MultiIndex structure, keeping entry and subentry as columns")
#                 df = df.reset_index()
#                 # new_columns['entry'] = 'entry'
#                 new_columns['subentry'] = 'subentry'
#                 if(indices): indices.append("subentry") # safe: ensure unique index for each row
#             df = df.rename(columns = new_columns)

#         if(indices):
#             df = df.set_index(indices)
        
#         return df
    
#     def Tuple2String(self, input_column: tuple) -> str:
#         """
#         parse input column (tuple) and return a short name
#         Example:
#         ('PrimariesFirstHitECAL', 'fP', 'fX') -> PrimariesFirstHitECAL_x
#         ('PrimariesFirstHitECAL', 'fE', '')   -> PrimariesFirstHitECAL_t,
#         """
#         colName, suffix1, suffix2 = input_column
#         if(suffix1=='fP'):
#             if(suffix2=='fX'):
#                 return colName+"_x"
#             elif(suffix2=='fY'):
#                 return colName+"_y"
#             elif(suffix2=='fZ'):
#                 return colName+"_z"
#             else:
#                 return colName
#         elif(suffix1=='fE'):
#             return colName+"_t"
#         else:
#             return colName

class Converter:
    def __init__(self, file_names: list, tree_name: str):
        """
        Initializes the Converter class.
        
        Parameters:
        - file_names (list): List of input ROOT files.
        - tree_name (str): Tree name in the ROOT input files.
        """
        for file_name in file_names:
            if not file_name.endswith('.root'):
                raise ValueError(f"Invalid file extension: {file_name}. The file must have a .root extension.")
        
        self.file_names = file_names
        self.tree_name = tree_name
    
    def keys(self, file_index: int = 0):
        """
        Print tree columns for the specified file.
        
        Parameters:
        - file_index (int): Index of the file in the list of file names. Default is 0 (first file).
        """
        if file_index < 0 or file_index >= len(self.file_names):
            raise IndexError("File index out of range.")
        
        tree = upr.open(self.file_names[file_index])[self.tree_name]
        print(tree.keys())

    def CreatePandas(self, columns: list, rename: bool = False, indices: list = None):
        """
        Create a Pandas DataFrame by concatenating the data from all the ROOT files.
        
        Parameters:
        - columns (list): List of column names (branches) to extract from the ROOT files.
        - rename (bool): If True, rename the columns using a simplified naming scheme.
        - indices (list): List of column names to set as the DataFrame index.
        
        Returns:
        - pd.DataFrame: A concatenated DataFrame containing the requested columns from all ROOT files.
        
        Raises:
        - TypeError: If the resulting object is not a Pandas DataFrame.
        """
        df_list = []
        
        for file_name in self.file_names:
            try:
                tree = upr.open(file_name)[self.tree_name]
                df = tree.arrays(columns, library='pd')
            except KeyError:
                print(f"\033[91m {self.tree_name} not found in {file_name}. Skipping file \033[0m")
                continue
            
            if not isinstance(df, pd.DataFrame):
                print(f"Result from file {file_name} is not a DataFrame. Here's the tuple returned:", df)
                raise TypeError(f"Expected a Pandas DataFrame from file {file_name}, but got a tuple instead.")
            
            if rename:
                new_columns = {col: self.Tuple2String(col) for col in df.columns}
                
                if isinstance(df.index, pd.MultiIndex):
                    print(f"DataFrame from file {file_name} has MultiIndex structure, keeping entry and subentry as columns")
                    df = df.reset_index()
                    new_columns['subentry'] = 'subentry'
                    if indices:
                        indices.append("subentry")  # Ensure unique index for each row
            
                df = df.rename(columns=new_columns)
            
            if indices:
                df = df.set_index(indices)
            
            df_list.append(df)
        
        # Concatenate all DataFrames
        concatenated_df = pd.concat(df_list, ignore_index=False)
        
        return concatenated_df
    
    def Tuple2String(self, input_column: tuple) -> str:
        """
        Parse input column (tuple) and return a short name.
        
        Example:
        ('PrimariesFirstHitECAL', 'fP', 'fX') -> PrimariesFirstHitECAL_x
        ('PrimariesFirstHitECAL', 'fE', '')   -> PrimariesFirstHitECAL_t
        """
        colName, suffix1, suffix2 = input_column
        if suffix1 == 'fP':
            if suffix2 == 'fX':
                return colName + "_x"
            elif suffix2 == 'fY':
                return colName + "_y"
            elif suffix2 == 'fZ':
                return colName + "_z"
            else:
                return colName
        elif suffix1 == 'fE':
            return colName + "_t"
        else:
            return colName
