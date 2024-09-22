import pandas as pd

class Sample:
    def __init__(self, name: str, condition: str, dataframe: pd.DataFrame, reference_index: str = None):
        """
        Initialize a Sample object.
        
        Parameters:
        - name (str): The name of the sample.
        - condition (str): The condition string used to filter the DataFrame.
        - dataframe (pd.DataFrame): The filtered DataFrame that matches the condition.
        - reference_index (str): Name of the DataFrame index to use to perform operations with different Samples.
        """
        self.name = name
        self.condition = condition
        self.dataframe = dataframe
        self.reference_index = reference_index
    
    def SampleIndex(self) -> pd.Index:
        """
        Retrieve the index values for the sample based on the reference index.
        
        Returns:
        - pd.Index: Index values based on the reference index if specified, otherwise the default index.
        """
        if self.reference_index is not None:
            return self.dataframe.index.get_level_values(self.reference_index)
        else:
            return self.dataframe.index
    
    def ExpandIndex(self, list_index: list) -> pd.Series:
        """
        Get a pd.Series of [True, False] of length len(self.SampleIndex)
        from a list_index.
        
        Parameters:
        - list_index (list): List of indices to check against the DataFrame index.
        
        Returns:
        - pd.Series: Boolean series of the same length as self.SampleIndex.
        """
        return self.SampleIndex().isin(list_index)
    
    def loc(self, list_rows: pd.Series) -> pd.DataFrame:
        """
        Locate and return rows of the DataFrame based on a boolean series.
        
        Parameters:
        - list_rows (pd.Series): A boolean series indicating which rows to keep.
        
        Returns:
        - pd.DataFrame: The filtered DataFrame.
        """
        return self.dataframe.loc[list_rows]
    
    def size(self):
        """
        Returns: sample size
        """
        return len(self.dataframe)
    
    def AddInfo(self, other_sample:"Sample", variable:str) -> pd.DataFrame:
        if(other_sample.reference_index != self.reference_index):
            raise ValueError(f"cannot find info in sample with diffrent reference index")
        return self.dataframe.merge(
            other_sample.dataframe[variable],
            on = self.reference_index,
            how ='left'
        )[variable]

class Manager:
    def __init__(self, dataframe: pd.DataFrame, manager_name: str, reference_index: str = None):
        """
        Initialize a Manager object.
        
        Parameters:
        - dataframe (pd.DataFrame): The DataFrame managed by this Manager.
        - manager_name (str): The name of the manager.
        - reference_index (str): Name of the DataFrame index used for referencing.
        """
        self.df = dataframe
        self.manager_name = manager_name
        self.reference_index = reference_index
        self.samples = {}
        
        if self.reference_index is not None:
            self.df_index = self.df.index.get_level_values(self.reference_index)
        else:
            self.df_index = self.df.index
    
    def DefineSample(self, sample_name: str, condition: str):
        """
        Define a sample using a condition string and add it to the manager's samples.
        
        Parameters:
        - sample_name (str): The name of the sample.
        - condition (str): The condition string used to filter the DataFrame.
        """
        self.samples[sample_name] = Sample(sample_name, condition, self.df.query(condition), self.reference_index)
    
    def GetSample(self, sample_name: str) -> Sample:
        """
        Retrieve a Sample object by its name.
        
        Parameters:
        - sample_name (str): The name of the sample to retrieve.
        
        Returns:
        - Sample: The Sample object associated with the given name.
        
        Raises:
        - ValueError: If the sample name is not defined.
        """
        if sample_name not in self.samples:
            raise ValueError(f"Sample '{sample_name}' not defined in manager '{self.manager_name}'")
        return self.samples[sample_name]
    
    def CombineSamples(self, this_sample_name: str, other_manager: "Manager", other_sample_name: str, new_sample_name: str):
        """
        Combine two samples from different managers based on the intersection of their indices
        and create a new sample in this manager.
        
        Parameters:
        - this_sample_name (str): The name of the sample in this manager.
        - other_manager (Manager): The other manager containing the sample to combine.
        - other_sample_name (str): The name of the sample in the other manager.
        - new_sample_name (str): The name of the new sample to create in this manager.
        
        Raises:
        - ValueError: If the samples do not share a common reference index.
        """
        this_sample = self.GetSample(this_sample_name)
        other_sample = other_manager.GetSample(other_sample_name)

        if this_sample.reference_index != other_sample.reference_index:
            raise ValueError(f"Samples '{this_sample_name}' and '{other_sample_name}' do not share a common index.")

        this_sample_index = this_sample.SampleIndex()
        other_sample_index = other_sample.SampleIndex()

        index_intersection = this_sample_index.intersection(other_sample_index)
        common_index_expanded = this_sample.ExpandIndex(index_intersection)

        combined_condition = f"({this_sample.condition}) & ({other_sample.condition})"

        self.samples[new_sample_name] = Sample(new_sample_name, combined_condition, this_sample.loc(common_index_expanded), this_sample.reference_index)