import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

class ExcelReader:
    def __init__(self, filename, sheet_name, start_row, end_row, columns_to_read):
        self.filename = filename
        self.sheet_name = sheet_name
        self.start_row = start_row
        self.end_row = end_row
        self.columns_to_read = columns_to_read

    def read_data(self):
        try:
            df = pd.read_excel(self.filename, self.sheet_name, skiprows=self.start_row, nrows=self.end_row - self.start_row + 1, usecols=self.columns_to_read)
            return df
        except FileNotFoundError:
            print(f"Error: File '{self.filename}' not found.")
            return None
        except Exception as e:
            print(f"Error: {e}")
            return None

diesel = ExcelReader('pruebas_lab.xlsx', 'CH4_1', 42, 49, 'C:O').read_data()

print(diesel)